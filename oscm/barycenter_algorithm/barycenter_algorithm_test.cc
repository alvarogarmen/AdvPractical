#include "barycenter_algorithm.h"

#include <gtest/gtest.h>

#include "ds/bipartite_graph.h"
#include "ds/reduction_graph/reduction_graph.h"
#include "gmock/gmock-matchers.h"

namespace {
using ::testing::ElementsAre;
}

TEST(BarycenterTest, TrivialGraph) {
  BipartiteGraph myGraph = BipartiteGraph(4, 4, 0);
  myGraph.addEdge(0, 0);
  myGraph.addEdge(1, 1);
  myGraph.addEdge(2, 2);
  myGraph.addEdge(3, 3);
  barycenter_algorithm::barycenterAlgorithm(myGraph);

  for (auto i = 0; i < myGraph.getFreeNodesSize(); i++) {
    EXPECT_EQ(myGraph.getFreeNodes()[i], i);
  }
}

TEST(BarycenterTest, TrivialGraphOutdegreeTwo) {
  BipartiteGraph myGraph = BipartiteGraph(8, 4, 0);

  myGraph.addEdge(0, 0);
  myGraph.addEdge(0, 1);
  myGraph.addEdge(1, 2);
  myGraph.addEdge(1, 3);
  myGraph.addEdge(2, 4);
  myGraph.addEdge(2, 5);
  myGraph.addEdge(3, 6);
  myGraph.addEdge(3, 7);

  barycenter_algorithm::barycenterAlgorithm(myGraph);
  for (auto i = 0; i < myGraph.getFreeNodesSize(); i++) {
    EXPECT_EQ(myGraph.getFreeNodes()[i], i);
  }
}

TEST(BarycenterTest, ReverseOrderGraph) {
  BipartiteGraph myGraph = BipartiteGraph(8, 4, 0);

  myGraph.addEdge(3, 0);
  myGraph.addEdge(3, 1);
  myGraph.addEdge(2, 2);
  myGraph.addEdge(2, 3);
  myGraph.addEdge(1, 4);
  myGraph.addEdge(1, 5);
  myGraph.addEdge(0, 6);
  myGraph.addEdge(0, 7);
  std::cout << myGraph.getEdges()[0][0] << std::endl;
  barycenter_algorithm::barycenterAlgorithm(myGraph);

  for (auto i = 0; i < myGraph.getFreeNodesSize(); i++) {
    EXPECT_EQ(myGraph.getFreeNodes()[i], myGraph.getFreeNodesSize() - 1 - i);
  }
}
