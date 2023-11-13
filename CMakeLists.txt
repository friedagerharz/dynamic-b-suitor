cmake_minimum_required (VERSION 3.17)
project(dyn-b-suitor-comparison CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(PkgConfig)
pkg_check_modules(NETWORKIT REQUIRED IMPORTED_TARGET networkit)

add_executable(dyn-b-suitor-comparison main.cpp)
target_include_directories(dyn-b-suitor-comparison PUBLIC "${PROJECT_SOURCE_DIR}/." ${NETWORKIT_INCLUDE_DIRS})
target_link_libraries(dyn-b-suitor-comparison PUBLIC PkgConfig::NETWORKIT)
install(TARGETS dyn-b-suitor-comparison DESTINATION ${CMAKE_INSTALL_PREFIX})

# link_directories(/usr/local/lib)

# # Link against the installed Networkit library
# target_link_libraries(dyn-b-suitor-comparison
#     PRIVATE networkit
# )

