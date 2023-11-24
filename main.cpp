#include <chrono>
#include <cstring>
#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/MatrixMarketGraphReader.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>
#include <networkit/matching/DynamicBSuitorMatcher.hpp>

using std::chrono::high_resolution_clock;
using std::chrono::seconds;
typedef std::chrono::duration<double> dur;

using namespace NetworKit;

std::string operation;
uint32_t batch_size;
std::optional<count> num_b;
std::optional<std::vector<count>> vec_b;
template <typename T> std::vector<T> edges;
Graph G;

std::vector<edgeweight> dyn_ws;
edgeweight stat_w;
std::vector<count> dyn_num_affected;

std::vector<dur> dyn_rt;
std::vector<dur> stat_rt;
// double dyn_mean = 0.0;
// double stat_mean = 0.0;
// double dyn_var;
// double stat_var;

int num_runs = 50;
count n;
// limits for the generation of random b values
[[maybe_unused]] count b_min = 1;
[[maybe_unused]] count b_max = 10;

// void calculateMean() {
//   for (int i = 0; i < dyn_rt.size(); i++) {
//     dyn_mean += dyn_rt.at(i).count();
//     stat_mean += stat_rt.at(i).count();
//   }
// }

// void calculateVariance() {
//   double dyn_sum_diff_sq = 0.0;
//   double stat_sum_diff_sq = 0.0;

//   for (int i = 0; i < dyn_rt.size(); i++) {
//     double dyn_diff = dyn_rt.at(i).count() - dyn_mean;
//     double stat_diff = stat_rt.at(i).count() - stat_mean;

//     dyn_sum_diff_sq += dyn_diff * dyn_diff;
//     stat_sum_diff_sq += stat_diff * stat_diff;
//   }

//   dyn_var = dyn_sum_diff_sq / dyn_rt.size();
//   stat_var = stat_sum_diff_sq / stat_rt.size();
// }

std::string pluralS(int num) {
  return (num == 1) ? "" : "s";
}
std::string fromOrInto(std::string &s) {
  return (s == "insert") ? "into" : "from";
}

std::string getFileFormat(const std::string &file) {
  auto pos = file.find_last_of(".");
  if (pos != std::string::npos) {
    return file.substr(pos + 1);
  }
  return "";
}

std::vector<count> generateRandomVector() {
  std::vector<count> rand_v;
  rand_v.reserve(n);
  for (int i = 0; i < n; i++) {
    rand_v.emplace_back(Aux::Random::integer(b_min, b_max));
  }
  return rand_v;
}

void printUse() {
  std::cerr << "Usage: ./dyn-b-suitor [0] [1] [2] [3]\n"
            << "[0]: path to graph file in MatrixMarket, METIS or graph-tool "
               "binary format\n"
            << "[1]: operation 'insert' or 'remove'\n"
            << "[2]: batch size\n"
            << "[3]: constant b value or 'random'\n";
}

bool parseInput(std::vector<std::string> args) {
  const auto file = args.at(1);
  const auto format = getFileFormat(file);
  if (format.empty()) {
    printUse();
    std::cerr << "first argument path must a path to a graph file" << std::endl;
    return false;
  } else if (format == "mtx") {
    G = MatrixMarketGraphReader{}.read(file);
    G.removeSelfLoops();
    G.removeMultiEdges();
  } else if (format == "gt") {
    G = GraphToolBinaryReader{}.read(
        file); // doesn't work with weighted graphs (GraphToolBinaryReader line
               // 23 weighted=false) therefore adding random weights below
    G.removeSelfLoops();
    G.removeMultiEdges();
    GraphTools::randomizeWeights(G);
  } else if (format == "graph") {
    G = METISGraphReader{}.read(file);
    G.removeSelfLoops();
    G.removeMultiEdges();
  } else {
    printUse();
    std::cerr << "invalid graph format" << std::endl;
    return false;
  }
  n = G.numberOfNodes();

  if (args.at(2) != "insert" && args.at(2) != "remove") {
    printUse();
    std::cerr << "second argument operation must be 'insert' or 'remove'"
              << std::endl;
    return false;
  } else {
    operation = args.at(2);
  }

  batch_size = std::stoi(args.at(3));

  std::istringstream iss(args.at(4));
  [[maybe_unused]] int num;
  if (iss >> num) {
    num_b = num;
  } else {
    vec_b = generateRandomVector();
  }

  return true;
}

template <typename EdgeType>
dur edgeInsertion(Graph &G, DynamicBSuitorMatcher &dbsm) {
  for (auto &edge : edges<WeightedEdge>) {
    G.addEdge(edge.u, edge.v, edge.weight);
  }

  const auto t1 = high_resolution_clock::now();
  dbsm.addEdges(edges<WeightedEdge>);
  dbsm.buildBMatching();
  const auto t2 = high_resolution_clock::now();
  return t2 - t1;
}

template <typename EdgeType>
dur edgeRemoval(Graph &G, DynamicBSuitorMatcher &dbsm) {
  for (auto e : edges<Edge>)
    G.removeEdge(e.u, e.v);
  const auto t1 = high_resolution_clock::now();
  dbsm.removeEdges(edges<Edge>);
  dbsm.buildBMatching();
  const auto t2 = high_resolution_clock::now();
  return t2 - t1;
}

template <typename BType> void runDynamicBSuitor(Graph &G, BType &b) {
  if (operation == "insert") {
    // Select batch_size edges of the graph, remove them but put them into edges
    // for later insertion. This will make sure that the graph is valid and th
    // after insertion.
    for (auto j = 0; j < batch_size; j++) {
      const auto [u, v] = GraphTools::randomEdge(G);
      assert(G.hasEdge(u, v));
      edges<WeightedEdge>.emplace_back(u, v, G.weight(u, v));
      G.removeEdge(u, v);
    }
  } else {
    // Select batch_size edges to be added to the graph. Put them into edges for
    // later removal. This will make sure that the graph is valid and th after
    // removal.

    for (auto j = 0; j < batch_size; j++) {
      node u, v;
      do {
        u = GraphTools::randomNode(G);
        v = GraphTools::randomNode(G);
      } while (u == v || G.hasEdge(u, v));
      edgeweight w = 1.0;
      edges<Edge>.emplace_back(u, v);
      G.addEdge(u, v, w);
    }
  }

  DynamicBSuitorMatcher dbsm(G, b);

  dbsm.run();

  dur dyn_t = (operation == "insert") ? edgeInsertion<WeightedEdge>(G, dbsm)
                                      : edgeRemoval<Edge>(G, dbsm);
  dyn_num_affected.emplace_back(dbsm.getNumberOfAffected());

  const auto dm = dbsm.getBMatching();

  dyn_ws.emplace_back(dm.weight(G));
  dyn_rt.emplace_back(dyn_t);
}

template <typename BType> void runStaticBSuitor(Graph &G, BType &b) {
  BSuitorMatcher bsm(G, b);

  const auto t3 = high_resolution_clock::now();
  bsm.run();
  bsm.buildBMatching();
  const auto t4 = high_resolution_clock::now();
  dur stat_t = t4 - t3;

  const auto sm = bsm.getBMatching();
  stat_w = sm.weight(G);

  stat_rt.emplace_back(stat_t);
}

void printResults() {
  std::cout << "(Dynamic) affected nodes per run:\n";
  for (auto v : dyn_num_affected) {
    std::cout << v << std::endl;
  }
  std::cout << std::endl;
  // std::cout << "(Dynamic) Adding " << batch_size << " edges took " <<
  // dyn_mean << "s on average.\n";
  std::cout << "(Dynamic) runtimes [s] per run:\n";
  for (auto r : dyn_rt) {
    std::cout << r.count() << std::endl;
  }
  std::cout << std::endl;
  // std::cout << "(Static) Running b-suitor took " << stat_mean << "s on
  // average.\n";
  std::cout << "(Static) runtimes [s] per run:\n";
  for (auto r : stat_rt) {
    std::cout << r.count() << std::endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    printUse();
    return 1;
  }
  Aux::Random::setSeed(0, true);

  if (!parseInput(std::vector<std::string>(argv, argv + argc))) {
    return 1;
  }

  std::cout << "Start comparison for the dynamic " << argv[4]
            << "-matching: " << operation << " " << batch_size << " edge"
            << pluralS(batch_size) << " " << fromOrInto(operation)
            << " the graph located at " << argv[1] << " with "
            << G.numberOfNodes() << " nodes and " << G.numberOfEdges()
            << " edges." << std::endl;
  std::cout << std::endl;

  for (int i = 0; i < num_runs; i++) {
    Aux::Random::setSeed(i, true);
    (operation == "insert") ? edges<WeightedEdge>.clear() : edges<Edge>.clear();
    num_b.has_value() ? runDynamicBSuitor(G, num_b.value())
                      : runDynamicBSuitor(G, vec_b.value());
  }
  num_b.has_value() ? runStaticBSuitor(G, num_b.value())
                    : runStaticBSuitor(G, vec_b.value());

  // calculateMean();
  // calculateVariance();

  printResults();

  for (auto w : dyn_ws)
    assert(std::abs(w - stat_w) < FLOAT_EPSILON);

  return 0;
}
