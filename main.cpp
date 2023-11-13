#include <chrono>
#include <cstring>
#include <filesystem>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/MatrixMarketGraphReader.hpp>
#include <networkit/matching/BSuitorMatcher.hpp>
#include <networkit/matching/DynamicBSuitorMatcher.hpp>

using std::chrono::high_resolution_clock;
using std::chrono::seconds;
typedef std::chrono::duration<double> dur;

using namespace NetworKit;

std::string operation;
uint32_t batch_size;
std::optional<int> num_b;
std::optional<std::string> path_b;
template <typename T> std::vector<T> edges;
Graph G;

[[maybe_unused]] count n = 100;
[[maybe_unused]] double p = 0.5;

std::string pluralS(int num) {
  return (num == 1) ? "" : "s";
}
std::string fromOrInto(const char *s) {
  return (std::strcmp(s, "insert") == 0) ? "into" : "from";
}

void printUse() {
  std::cerr << "Usage: ./dyn-b-suitor [0] [1] [2] [3]\n"
            << "[0]: path to mtx graph file or 'random' for using a random "
               "generated graph\n"
            << "[1]: operation 'insert' or 'remove'\n"
            << "[2]: batch size\n"
            << "[3]: constant b value or path to file with n (number of "
               "nodes in the graph) b-values (only if graph provided)"
            << std::endl;
}

bool parseInput(std::vector<std::string> args) {
  if (std::strcmp(args.at(1).c_str(), "random") == 0) {
    // Generate a random undirected Graph with n nodes, a chance of p for an
    // edge  between a pair of nodes and without self-loops. Set a random weight
    // between 1 and 100 for each edge.
    G = ErdosRenyiGenerator(n, p, false, false).generate();
    G = GraphTools::toWeighted(G);
    G.forEdges([&](node u, node v) {
      G.setWeight(u, v, Aux::Random::integer(1, 100));
    });
  } else {
    G = MatrixMarketGraphReader{}.read(args.at(1));
    G.removeSelfLoops();
    G.removeMultiEdges();
  }

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
    path_b = iss.str();
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
  const auto t1 = high_resolution_clock::now();
  dbsm.removeEdges(edges<Edge>);
  dbsm.buildBMatching();
  const auto t2 = high_resolution_clock::now();
  return t2 - t1;
}

template <typename EdgeType> void runComparison(Graph &G, EdgeType b) {
  // Select m edges of the graph, remove them but put them into edges for
  // later insertion. This will make sure that the graph is valid.
  for (auto j = 0; j < batch_size; j++) {
    const auto [u, v] = GraphTools::randomEdge(G);
    assert(G.hasEdge(u, v));
    (std::strcmp(operation.c_str(), "insert") == 0)
        ? edges<WeightedEdge>.emplace_back(u, v, G.weight(u, v))
        : edges<Edge>.emplace_back(u, v);
    G.removeEdge(u, v);
  }

  DynamicBSuitorMatcher dbsm(G, b.value());
  dbsm.run();

  dur dyn_t = (std::strcmp(operation.c_str(), "insert") == 0)
                  ? edgeInsertion<WeightedEdge>(G, dbsm)
                  : edgeRemoval<Edge>(G, dbsm);

  const auto dm = dbsm.getBMatching();
  const auto dwm = dm.weight(G);

  BSuitorMatcher bsm(G, b.value());

  const auto t3 = high_resolution_clock::now();
  bsm.run();
  bsm.buildBMatching();
  const auto t4 = high_resolution_clock::now();
  dur stat_t = t4 - t3;

  const auto sm = bsm.getBMatching();
  const auto wm = sm.weight(G);

  assert(dwm == wm);

  std::cout << "(Dynamic) Adding " << batch_size << " edge"
            << pluralS(batch_size) << " took\n"
            << dyn_t.count() << " s.\n";
  std::cout << "(Static) Running b-Suitor took\n" << stat_t.count() << " s.\n";
}

int main(int argc, char *argv[]) {
  std::cout << "Start runtime comparison for the dynamic " << argv[4]
            << "-matching: " << argv[2] << " " << argv[3] << " edge"
            << pluralS(std::atoi(argv[3])) << " " << fromOrInto(argv[2])
            << " graph " << std::filesystem::path(argv[1]).filename()
            << std::endl;
  if (argc != 5) {
    std::cout << argc << std::endl;
    printUse();
    return 1;
  }
  Aux::Random::setSeed(0, true);
  // Graph G;

  if (!parseInput(std::vector<std::string>(argv, argv + argc))) {
    return 1;
  }

  num_b.has_value() ? runComparison(G, num_b) : runComparison(G, path_b);

  return 0;
}
