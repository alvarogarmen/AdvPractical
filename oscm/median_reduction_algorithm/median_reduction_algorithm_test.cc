#include "median_reduction_algorithm.h"

#include <map>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
namespace {

TEST(MedianReductionTest, algorithm) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  auto [crossingSum, orderVector] =
      medianReduction_algorithm::medianReductionAlgorithm<ReductionGraph<int, int>,
                                                          UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(crossingSum, 1);
  EXPECT_EQ(orderVector[0], 1);
  EXPECT_EQ(orderVector[1], 0);
  EXPECT_EQ(orderVector[2], 2);
}

TEST(medianReductionTest, algorithmV1) {
  std::vector<std::vector<int>> freeNodes = {{5},    {6}, {7}, {8}, {0, 9},
                                             {0, 9}, {1}, {2}, {3}, {4}};
  std::vector<std::vector<int>> fixedNodes = {{4, 5}, {6}, {7}, {8}, {9},
                                              {0},    {1}, {2}, {3}, {4, 5}};
  ReductionGraph myGraph = ReductionGraph<int, int>(freeNodes, fixedNodes);
  auto [crossingSum, orderVector] =
      medianReduction_algorithm::medianReductionAlgorithm<ReductionGraph<int, int>,
                                                          UndoAlgorithmStep<int, int>>(myGraph);
  EXPECT_EQ(crossingSum, 17);
}
}  // namespace