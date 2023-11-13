#include "h_graph.h"

#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraphTest, SimpleTest) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  HGraph myGraph = HGraph<int, int>(freeNodes, fixedNodes);

  EXPECT_EQ(myGraph.getFreeNodesSize(), 3);
  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph.getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph.getFreeNodeNeighbours(0)[1], 1);

  EXPECT_EQ(myGraph.getLeftCrossings(0), 0);
  EXPECT_EQ(myGraph.getRightCrossings(0), 2);
  EXPECT_EQ(myGraph.getLeftCrossings(1), 1);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 1);
  EXPECT_EQ(myGraph.getRightCrossings(2), 0);
  myGraph.switchNeighbours(0, 1);
  EXPECT_EQ(myGraph.getLeftCrossings(1), 0);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);
  EXPECT_EQ(myGraph.getPermutation()[0], 1);
  EXPECT_EQ(myGraph.getPermutation()[1], 0);
  EXPECT_EQ(myGraph.getPermutation()[2], 2);

  myGraph.switchNeighbours(0, 2);
  EXPECT_EQ(myGraph.getLeftCrossings(0), 3);
  EXPECT_EQ(myGraph.getRightCrossings(0), 0);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 0);
  EXPECT_EQ(myGraph.getRightCrossings(2), 3);
}