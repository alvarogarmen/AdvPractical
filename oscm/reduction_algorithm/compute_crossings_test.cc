#include <map>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include "reduction_algorithm.h"

TEST(ComputeCrossings, computeUVcrossing) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();
  EXPECT_EQ((computeUVcrossing<ReductionGraph<int, int>>(myGraph, 0, 2)), 1);
  EXPECT_EQ((computeUVcrossing<ReductionGraph<int, int>>(myGraph, 1, 0)), 0);
  EXPECT_EQ((computeUVcrossing<ReductionGraph<int, int>>(myGraph, 2, 0)), 3);
}

TEST(ComputeCrossings, computeCrossingSum) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  EXPECT_EQ(myGraph.getFreeNodesSize(), 3);
  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph.getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph.getFreeNodeNeighbours(0)[1], 1);
  EXPECT_EQ(myGraph.getNodeCrossing(0).at(2), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(2).at(0), 3);
  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(1).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 1);

  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 2);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);
  EXPECT_EQ(*(myGraph.getLeftNodes(0).find(1)), 1);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(0)), 0);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(2)), 2);

  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(1).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 1);
}