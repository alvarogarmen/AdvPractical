#include "barycenter_algorithm.h"

#include <gtest/gtest.h>

#include "ds/bipartite_graph.h"
#include "gmock/gmock-matchers.h"

namespace {
using ::testing::ElementsAre;
}

TEST(BarycenterTest, TrivialGraph) {
  BipartiteGraph myGraph = BipartiteGraph(4, 4, 0);

  myGraph.addEdge(4, 1);
  myGraph.addEdge(5, 2);
  myGraph.addEdge(6, 3);
  myGraph.addEdge(7, 4);

  barycenterAlgorithm(myGraph);
  for (auto i = 0; i < myGraph.getFreeNodesSize(); i++) {
    EXPECT_EQ(myGraph.getFreeNodes()[i], i + 1);
  }
}
