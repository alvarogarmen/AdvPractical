#include "heuristic_algorithm.h"

#include <map>
#include <vector>

#include "ds/heuristic_graph/heuristic_graph.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include "oscm/heuristic_algorithm/heuristic_algorithm.h"
using namespace heuristic_algorithm;
TEST(HeuristicTest, HeuristicWithR1) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  HeuristicGraph myGraph = HeuristicGraph<int, int>(freeNodes, fixedNodes);

  bool didSwitch = heuristicAlgorithm<HeuristicGraph<int, int>>(myGraph, true, false, false);
  EXPECT_EQ(didSwitch, true);
  EXPECT_EQ(myGraph.getPermutation()[0], 1);
  EXPECT_EQ(myGraph.getPermutation()[1], 0);
  EXPECT_EQ(myGraph.getPermutation()[2], 2);

  EXPECT_EQ(myGraph.getLeftCrossings(1), 0);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);

  EXPECT_EQ(myGraph.getLeftCrossings(0), 0);
  EXPECT_EQ(myGraph.getRightCrossings(0), 1);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 1);
  EXPECT_EQ(myGraph.getRightCrossings(2), 0);
}

TEST(HeuristicTest, HeuristicWithR2) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  HeuristicGraph myGraph = HeuristicGraph<int, int>(freeNodes, fixedNodes);

  bool didSwitch = heuristicAlgorithm<HeuristicGraph<int, int>>(myGraph, false, true, false);
  EXPECT_EQ(didSwitch, true);
  EXPECT_EQ(myGraph.getPermutation()[0], 1);
  EXPECT_EQ(myGraph.getPermutation()[1], 0);
  EXPECT_EQ(myGraph.getPermutation()[2], 2);

  EXPECT_EQ(myGraph.getLeftCrossings(1), 0);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);

  EXPECT_EQ(myGraph.getLeftCrossings(0), 0);
  EXPECT_EQ(myGraph.getRightCrossings(0), 1);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 1);
  EXPECT_EQ(myGraph.getRightCrossings(2), 0);
}

TEST(HeuristicTest, HeuristicWithR3) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  HeuristicGraph myGraph = HeuristicGraph<int, int>(freeNodes, fixedNodes);

  bool didSwitch = heuristicAlgorithm<HeuristicGraph<int, int>>(myGraph, false, false, true);
  EXPECT_EQ(didSwitch, true);
  EXPECT_EQ(myGraph.getPermutation()[0], 1);
  EXPECT_EQ(myGraph.getPermutation()[1], 0);
  EXPECT_EQ(myGraph.getPermutation()[2], 2);

  EXPECT_EQ(myGraph.getLeftCrossings(1), 0);
  EXPECT_EQ(myGraph.getRightCrossings(1), 0);

  EXPECT_EQ(myGraph.getLeftCrossings(0), 0);
  EXPECT_EQ(myGraph.getRightCrossings(0), 1);
  EXPECT_EQ(myGraph.getLeftCrossings(2), 1);
  EXPECT_EQ(myGraph.getRightCrossings(2), 0);
}