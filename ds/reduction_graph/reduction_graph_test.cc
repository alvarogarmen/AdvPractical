#include "reduction_graph.h"

#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include "undo_algorithm_step.h"

TEST(GraphTest, SimpleTest) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  EXPECT_EQ(myGraph.getFreeNodesSize(), 3);
  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph.getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph.getFreeNodeNeighbours(0)[1], 1);

  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);

  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(1).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 0);
}