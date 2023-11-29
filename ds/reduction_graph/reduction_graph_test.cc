#include "reduction_graph.h"

#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include "undo.h"

TEST(GraphTest, SimpleTest) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
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

  myGraph.parameterAccounting(2, 0, currentSolution);
  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 0);

  std::vector<std::vector<int>> freeNodes1 = {{0, 1}, {0}, {0, 1, 2}, {0}};
  std::vector<std::vector<int>> fixedNodes1 = {{0, 1, 2, 3}, {0, 2}, {2}};
  ReductionGraph myGraph1 = ReductionGraph<int, int>(freeNodes1, fixedNodes1);
  std::set<int> aLeftSet = {2};
  std::set<int> aRightSet = {};
  std::set<int> bLeftSet = {};
  std::set<int> bRightSet = {3};
  std::set<int> cLeftSet = {};
  std::set<int> cRightSet = {0};
  std::set<int> dLeftSet = {1};
  std::set<int> dRightSet = {};

  myGraph1.setLeftNodes(0, aLeftSet);
  myGraph1.setRightNodes(0, aRightSet);
  myGraph1.setLeftNodes(1, bLeftSet);
  myGraph1.setRightNodes(1, bRightSet);
  myGraph1.setLeftNodes(2, cLeftSet);
  myGraph1.setRightNodes(2, cRightSet);
  myGraph1.setLeftNodes(3, dLeftSet);
  myGraph1.setRightNodes(3, dRightSet);

  EXPECT_EQ(myGraph1.getLeftNodes(0).size(), 1);

  EXPECT_EQ(*(myGraph1.getLeftNodes(0).find(2)), 2);
  EXPECT_EQ(*(myGraph1.getRightNodes(1).find(3)), 3);
  EXPECT_EQ(*(myGraph1.getRightNodes(2).find(0)), 0);
  EXPECT_EQ(*(myGraph1.getLeftNodes(3).find(1)), 1);
  myGraph1.parameterAccounting(0, 1, currentSolution);

  EXPECT_EQ(myGraph1.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph1.getRightNodes(0).size(), 2);
  EXPECT_EQ(*(myGraph1.getLeftNodes(0).find(2)), 2);
  EXPECT_EQ(*(myGraph1.getRightNodes(0).find(1)), 1);
  EXPECT_EQ(*(myGraph1.getRightNodes(0).find(3)), 3);

  EXPECT_EQ(myGraph1.getLeftNodes(1).size(), 2);
  EXPECT_EQ(myGraph1.getRightNodes(1).size(), 1);
  EXPECT_EQ(*(myGraph1.getLeftNodes(1).find(2)), 2);
  EXPECT_EQ(*(myGraph1.getLeftNodes(1).find(0)), 0);
  EXPECT_EQ(*(myGraph1.getRightNodes(1).find(3)), 3);

  EXPECT_EQ(myGraph1.getLeftNodes(2).size(), 0);
  EXPECT_EQ(myGraph1.getRightNodes(2).size(), 3);
  EXPECT_EQ(*(myGraph1.getRightNodes(2).find(0)), 0);
  EXPECT_EQ(*(myGraph1.getRightNodes(2).find(1)), 1);
  EXPECT_EQ(*(myGraph1.getRightNodes(2).find(3)), 3);

  EXPECT_EQ(myGraph1.getLeftNodes(3).size(), 3);
  EXPECT_EQ(myGraph1.getRightNodes(3).size(), 0);

  EXPECT_EQ(myGraph1.getRightNodes(3).size(), 0);
  EXPECT_EQ(*(myGraph1.getLeftNodes(3).find(0)), 0);
  EXPECT_EQ(*(myGraph1.getLeftNodes(3).find(1)), 1);
  EXPECT_EQ(*(myGraph1.getLeftNodes(3).find(2)), 2);
}
TEST(GraphTest, undo) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  Undo undo = Undo<int, int>();
  int currentSolution = 0;

  myGraph.parameterAccounting(0, 2, currentSolution, &undo);
  EXPECT_EQ(undo.getParameterAccountingUndo().size(), 1);
  EXPECT_EQ(undo.getParameterAccountingUndo()[0].leftNode, 0);
  EXPECT_EQ(undo.getParameterAccountingUndo()[0].rightNode, 2);
  EXPECT_EQ(undo.getParameterAccountingUndo()[0].leftRightCrossing, 1);
  EXPECT_EQ(undo.getParameterAccountingUndo()[0].rightLeftCrossing, 3);
  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(1).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 2);
  myGraph.doUndo(undo);
  EXPECT_EQ(myGraph.getNodeCrossing(0).size(), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(1).size(), 0);
  EXPECT_EQ(myGraph.getNodeCrossing(2).size(), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(0).at(2), 1);
  EXPECT_EQ(myGraph.getNodeCrossing(2).at(0), 3);

  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 2);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);
}