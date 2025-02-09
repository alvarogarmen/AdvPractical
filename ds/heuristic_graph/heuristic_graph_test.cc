#include "heuristic_graph.h"

#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraphTest, SimpleTest) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  HeuristicGraph myGraph = HeuristicGraph<int, int>(freeNodes, fixedNodes);

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
  myGraph.switchNeighbours(0, 1, false);
  EXPECT_EQ(myGraph.getLeftCrossings(1), 0);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);
  EXPECT_EQ(myGraph.getPermutation()[0], 1);
  EXPECT_EQ(myGraph.getPermutation()[1], 0);
  EXPECT_EQ(myGraph.getPermutation()[2], 2);

  myGraph.switchNeighbours(0, 2, false);
  EXPECT_EQ(myGraph.getLeftCrossings(0), 3);
  EXPECT_EQ(myGraph.getRightCrossings(0), 0);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 0);
  EXPECT_EQ(myGraph.getRightCrossings(2), 3);
}

TEST(GraphTest, addEdge) {
  HeuristicGraph myGraph = HeuristicGraph<int, int>(3, 3, 0);
  myGraph.addEdge(0, 0);
  myGraph.addEdge(0, 1);
  myGraph.addEdge(1, 0);
  myGraph.addEdge(2, 0);
  myGraph.addEdge(2, 1);
  myGraph.addEdge(2, 2);
  EXPECT_EQ(myGraph.getFreeNodesSize(), 3);
  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph.getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph.getFreeNodeNeighbours(0)[1], 1);
}

TEST(GraphTest, setFreeNodes) {
  HeuristicGraph myGraph = HeuristicGraph<int, int>(3, 3, 0);
  myGraph.addEdge(0, 0);
  myGraph.addEdge(0, 1);
  myGraph.addEdge(1, 0);
  myGraph.addEdge(2, 0);
  myGraph.addEdge(2, 1);
  myGraph.addEdge(2, 2);

  std::vector<int> permutation = {2, 1, 0};
  myGraph.setFreeNodes(permutation);
  EXPECT_EQ(myGraph.getFreeNodesPosition().size(), 3);
  EXPECT_EQ(myGraph.getFreeNodesPosition()[0], 2);
  EXPECT_EQ(myGraph.getFreeNodesPosition()[1], 1);
  EXPECT_EQ(myGraph.getFreeNodesPosition()[2], 0);
  EXPECT_EQ(myGraph.getPermutation().size(), 3);
  EXPECT_EQ(myGraph.getPermutation()[0], 2);
  EXPECT_EQ(myGraph.getPermutation()[1], 1);
  EXPECT_EQ(myGraph.getPermutation()[2], 0);
}

TEST(GraphTest, getCrossings) {
  std::vector<std::vector<int>> freeNodes = {{5},    {6}, {7}, {8}, {0, 9},
                                             {0, 9}, {1}, {2}, {3}, {4}};
  std::vector<std::vector<int>> fixedNodes = {{4, 5}, {6}, {7}, {8}, {9},
                                              {0},    {1}, {2}, {3}, {4, 5}};
  HeuristicGraph myGraph = HeuristicGraph<int, int>(freeNodes, fixedNodes);
  EXPECT_EQ(myGraph.getCrossings(), 33);
}

