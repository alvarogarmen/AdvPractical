#include "ds/bipartite_graph.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(BipartiteTest, SizeTest) {
  int freeNodes = 4;
  int fixedNodes = 3;
  BipartiteGraph myGraph = BipartiteGraph<int>(fixedNodes, freeNodes, 4);

  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph.getFreeNodesSize(), 4);
}

TEST(BipartiteTest, PositionTest) {
  int freeNodes = 4;
  int fixedNodes = 3;
  BipartiteGraph myGraph = BipartiteGraph<int>(fixedNodes, freeNodes, 4);
  std::vector<int> storedFreeNodes = myGraph.getFreeNodes();
  std::vector<int> storedFixedNodes = myGraph.getFixedNodes();

  for (size_t i = 0; i < storedFreeNodes.size(); i++) {
    EXPECT_EQ(storedFreeNodes[i], i);
  }
  for (size_t i = 0; i < storedFixedNodes.size(); i++) {
    EXPECT_EQ(storedFixedNodes[i], i);
  }
}
TEST(BipartiteTest, AddEdgeTest) {
  int freeNodes = 4;
  int fixedNodes = 3;
  BipartiteGraph myGraph = BipartiteGraph<int>(fixedNodes, freeNodes, 4);

  myGraph.addEdge(3, 0);
  myGraph.addEdge(4, 2);

  std::vector<std::vector<int>> edges = myGraph.getEdges();
  // Reminder that nodeIDs of freeNodes start at fixedNodes.size(), hence 3 corresponds here to the
  // index 0
  EXPECT_EQ(edges[0][0], 0);
  EXPECT_EQ(edges[1][0], 2);
}

TEST(BipartiteTest, SwitchNodesTest) {
  int freeNodes = 4;
  int fixedNodes = 3;
  BipartiteGraph myGraph = BipartiteGraph<int>(fixedNodes, freeNodes, 4);

  myGraph.switchNodes(3, 4);
  std::vector<int> gottenFreeNodes = myGraph.getFreeNodes();
  // Reminder that nodeIDs of freeNodes start at fixedNodes.size(), hence 3 corresponds here to the
  // index 0
  EXPECT_EQ(gottenFreeNodes[0], 1);
  EXPECT_EQ(gottenFreeNodes[1], 0);
}