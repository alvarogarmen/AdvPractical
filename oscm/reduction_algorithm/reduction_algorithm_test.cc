#include "reduction_algorithm.h"

#include <map>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(AlgorithmTest, rrlo1Rrlo2) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int currentSolution = 0;
  UndoAlgorithmStep undo = UndoAlgorithmStep<int, int>();

  rrlo1<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &undo);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
  EXPECT_EQ(myGraph.getFixedPosition()[1], 0);
  EXPECT_EQ(myGraph.getFixedPosition()[2], 0);
  rrlo2<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &currentSolution, &undo);
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

  rr3<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, &currentSolution, &undo);
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

  rr2<ReductionGraph<int, int>>(myGraph, &currentSolution);
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
  rrLarge<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph, 2, &currentSolution);
  EXPECT_EQ(myGraph.getLeftNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(0).size(), 1);
  EXPECT_EQ(myGraph.getRightNodes(2).size(), 0);
  EXPECT_EQ(myGraph.getLeftNodes(2).size(), 2);
}

TEST(AlgorithmTest, IJBiggerThenFour) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int u = 0;
  int v = 0;
  IJBiggerThenFour<ReductionGraph<int, int>>(myGraph, &u, &v);
  EXPECT_EQ(u, 0);
  EXPECT_EQ(v, 2);
}

TEST(AlgorithmTest, IJEqualToThree) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int u = 0;
  int v = 0;
  IJEqualToThree<ReductionGraph<int, int>>(myGraph, &u, &v);
  EXPECT_EQ(u, 2);
  EXPECT_EQ(v, 0);
}

TEST(AlgorithmTest, IJEqualToTwo) {
  std::vector<std::vector<int>> freeNodes = {{1, 2}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{1, 2}, {0, 2}, {0, 2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  int u = 0;
  int v = 0;
  IJEqualToTwo<ReductionGraph<int, int>>(myGraph, &u, &v);
  EXPECT_EQ(u, 0);
  EXPECT_EQ(v, 2);
}

TEST(AlgorithmTest, algorithm) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(myGraph.getBestSolution(), 1);
  EXPECT_EQ(myGraph.getBestOrder()[0], 1);
  EXPECT_EQ(myGraph.getBestOrder()[1], 0);
  EXPECT_EQ(myGraph.getBestOrder()[2], 2);
}