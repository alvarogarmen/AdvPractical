#include "r_algorithm.h"

#include <map>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(AlgorithmTest, rrlo1Rrlo2) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  RGraph myGraph = RGraph<int, int>(freeNodes, fixedNodes);
  ReductionAlgorithm<int, int>::rrlo1(myGraph);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
  EXPECT_EQ(myGraph.getFixedPosition()[1], 0);
  EXPECT_EQ(myGraph.getFixedPosition()[2], 0);
  ReductionAlgorithm<int, int>::rrlo2(myGraph);
  ReductionAlgorithm<int, int>::rrlo1(myGraph);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
  EXPECT_EQ(myGraph.getFixedPosition()[1], 0);
  EXPECT_EQ(myGraph.getFixedPosition()[2], 2);
}
TEST(AlgorithmTest, rr3) {
  std::vector<std::vector<int>> freeNodes = {{0, 2}, {0, 1}, {0, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {1}, {0, 2}};
  RGraph myGraph = RGraph<int, int>(freeNodes, fixedNodes);
  ReductionAlgorithm<int, int>::rr3(myGraph);
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
  ReductionAlgorithm<int, int>::rrlo1(myGraph);
  EXPECT_EQ(myGraph.getFixedPosition()[0], 1);
}