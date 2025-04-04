#include "reduction_algorithm.h"

#include <map>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
namespace {
using namespace reductionalgorithms;
TEST(AlgorithmTest, ParameterAccounting) {
  std::vector<std::vector<int>> freeNodes = {{0, 2}, {0, 3}, {0, 1}, {1, 3}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {2, 3}, {0}, {1, 3}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();
  int currentSolution = 0;

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  parameterAccounting<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, 0, 1,
                                                                             currentSolution);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 1);

  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 0);

  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 1);

  EXPECT_EQ(myGraph.getLeftNodes(3).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(3).size(), 0);
  EXPECT_EQ(*(myGraph.getRightNodes(0).find(1)), 1);
  EXPECT_EQ(*(myGraph.getLeftNodes(1).find(0)), 0);
  EXPECT_EQ(currentSolution, 1);
  // now check transitvity
  parameterAccounting<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, 1, 2,
                                                                             currentSolution);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 3);

  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 2);

  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 2);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 1);

  EXPECT_EQ(myGraph.getLeftNodes(3).size(), 3);
  EXPECT_EQ(myGraph.getRightNodes(3).size(), 0);
  EXPECT_EQ(*(myGraph.getRightNodes(0).find(1)), 1);
  EXPECT_EQ(*(myGraph.getRightNodes(0).find(2)), 2);
  EXPECT_EQ(*(myGraph.getRightNodes(0).find(3)), 3);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(2)), 2);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(3)), 3);
  EXPECT_EQ(*(myGraph.getLeftNodes(1).find(0)), 0);
  EXPECT_EQ(*(myGraph.getRightNodes(2).find(3)), 3);
  EXPECT_EQ(*(myGraph.getLeftNodes(2).find(0)), 0);
  EXPECT_EQ(*(myGraph.getLeftNodes(2).find(1)), 1);
  EXPECT_EQ(*(myGraph.getLeftNodes(3).find(0)), 0);
  EXPECT_EQ(*(myGraph.getLeftNodes(3).find(1)), 1);
  EXPECT_EQ(*(myGraph.getLeftNodes(3).find(2)), 2);
  EXPECT_EQ(currentSolution, 7);
}

TEST(AlgorithmTest, rrlo1Rrlo2) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  rrlo1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &undo);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
  EXPECT_EQ(myGraph.getFixedPosition()[1], 0);
  EXPECT_EQ(myGraph.getFixedPosition()[2], 0);
  rrlo2<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution, &undo);
  rrlo1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &undo);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
  EXPECT_EQ(myGraph.getFixedPosition()[1], 0);
  EXPECT_EQ(myGraph.getFixedPosition()[2], 2);
}

TEST(AlgorithmTest, rr3) {
  std::vector<std::vector<int>> freeNodes = {{0, 2}, {0, 1}, {0, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {1}, {0, 2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();
  int currentSolution = 0;
  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  rr3<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution, &undo);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 0);
  EXPECT_EQ(*(myGraph.getLeftNodes(0).find(1)), 1);

  EXPECT_EQ(myGraph.getLeftNodes(1).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(1).size(), 2);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(0)), 0);
  EXPECT_EQ(*(myGraph.getRightNodes(1).find(2)), 2);

  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);
  EXPECT_EQ(*(myGraph.getLeftNodes(2).find(1)), 1);
  rrlo1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &undo);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
}

TEST(AlgorithmTest, rr2) {
  std::vector<std::vector<int>> freeNodes = {{0, 2}, {0, 1}, {0, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {1}, {0, 2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  rr2<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 0);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 1);

  EXPECT_EQ(*(myGraph.getRightNodes(0).find(2)), 2);

  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);

  EXPECT_EQ(*(myGraph.getLeftNodes(2).find(0)), 0);
}

TEST(AlgorithmTest, rrLarge) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  rrLarge<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, 2, currentSolution);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 2);
}

TEST(AlgorithmTest, IJBiggerThenFour) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  auto [best, u, v] = IJBiggerThenFour<ReductionGraph<int, int>>(myGraph);
  EXPECT_EQ(u, 0);
  EXPECT_EQ(v, 2);
}

TEST(AlgorithmTest, IJEqualToThree) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  auto [success, u, v] = IJEqualToThree<ReductionGraph<int, int>>(myGraph);
  EXPECT_EQ(u, 2);
  EXPECT_EQ(v, 0);
}

TEST(AlgorithmTest, IJEqualToTwo) {
  std::vector<std::vector<int>> freeNodes = {{1, 2}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{1, 2}, {0, 2}, {0, 2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);
  auto [success, u, v] = IJEqualToTwo<ReductionGraph<int, int>>(myGraph);
  EXPECT_EQ(u, 0);
  EXPECT_EQ(v, 2);
}

TEST(AlgorithmTest, algorithm) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  auto [crossingSum, orderVector] =
      algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(crossingSum, 1);
  EXPECT_EQ(orderVector[0], 1);
  EXPECT_EQ(orderVector[1], 0);
  EXPECT_EQ(orderVector[2], 2);
}

TEST(AlgorithmTest, algorithmV1) {
  std::vector<std::vector<int>> freeNodes = {{5},    {6}, {7}, {8}, {0, 9},
                                             {0, 9}, {1}, {2}, {3}, {4}};
  std::vector<std::vector<int>> fixedNodes = {{4, 5}, {6}, {7}, {8}, {9},
                                              {0},    {1}, {2}, {3}, {4, 5}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  auto [crossingSum, orderVector] =
      algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(crossingSum, 17);
}

TEST(GraphTest, UndoAlgorithmStep) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();
  int currentSolution = 0;

  computeCrossingSums<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  rr1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, currentSolution);

  parameterAccounting<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(
      myGraph, 0, 2, currentSolution, &undo);
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

TEST(AlgorithmTest, getCrossings) {
  std::vector<std::vector<int>> freeNodes = {{5},    {6}, {7}, {8}, {0, 9},
                                             {0, 9}, {1}, {2}, {3}, {4}};
  std::vector<std::vector<int>> fixedNodes = {{4, 5}, {6}, {7}, {8}, {9},
                                              {0},    {1}, {2}, {3}, {4, 5}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  auto [crossingSum, orderVector] =
      algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(crossingSum, 17);
  EXPECT_EQ(myGraph.getCrossings(), 17);

}

}