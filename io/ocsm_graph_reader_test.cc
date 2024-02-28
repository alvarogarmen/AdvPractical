#include "ocsm_graph_reader.h"

#include <gtest/gtest.h>

#include "ds/bipartite_graph.h"
#include "ds/heuristic_graph/heuristic_graph.h"
#include "ds/reduction_graph/reduction_graph.h"
#include "gmock/gmock-matchers.h"
#include "oscm/reduction_algorithm/reduction_algorithm.h"

namespace {
using ::testing::ElementsAre;
}

TEST(ReadGraphTest, ValidGraph) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 4 5\n"
      "1 5\n"
      "2 6\n"
      "3 4\n"
      "2 7\n"
      "3 5\n");
  auto result = readGraph<ReductionGraph<int, int>>(testData);
  EXPECT_EQ(result.status(), absl::OkStatus());
  auto graph = std::move(result.value());

  // Check n0, n1, and m
  EXPECT_EQ(graph->getEdge(1, 0),
            0);  // Target of first edge of node 1 (5) is 0 (we use 0-indexation)
  EXPECT_EQ(graph->getEdge(2, 0),
            1);  // Target of first edge of node 2 (6) is 1 (we use 0-indexation)
  EXPECT_EQ(graph->getEdge(0, 0),
            2);  // Target of first edge of node 0 (4) is 2 (we use 0-indexation)

  ASSERT_EQ(graph->getFixedNodesSize(), 3);  // We have 3 fixed Nodes

  ASSERT_EQ(graph->getFreeNodesSize(), 4);  // We have 4 free Nodes
}

// Test case for invalid header format
TEST(ReadGraphTest, InvalidHeaderFormat) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "invalid header\n");
  // Invalid header format: missing 'p' or incorrect format

  auto result = readGraph<ReductionGraph<int, int>>(testData);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for nodes out of bounds
TEST(ReadGraphTest, NodesOutOfBounds) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 4 5\n"
      "4 2\n");
  // Nodes out of bounds: source node greater than n0
  // Source node is out of bounds

  auto result = readGraph<ReductionGraph<int, int>>(testData);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for negative number of nodes or edges
TEST(ReadGraphTest, NegativeNumbers) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr -3 4 5\n");
  // Negative number of nodes
  // Negative n0

  auto result = readGraph<ReductionGraph<int, int>>(testData);
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for missing 'p' line
TEST(ReadGraphTest, MissingPLine) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "1 2\n");
  // Missing 'p' line
  // Edge line without 'p' line

  auto result = readGraph<ReductionGraph<int, int>>(testData);
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for file with invalid content
TEST(ReadGraphTest, InvalidFileContent) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr abc 4 5\n");
  // Invalid content: non-integer characters instead of number of nodes/edges
  // Invalid n0

  auto result = readGraph<ReductionGraph<int, int>>(testData);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case to check if the function handles errors correctly.
TEST(ReadGraphTest, MissingInput) {
  // Test with a non-existing file/empty stringstream
  std::stringstream nonExistentFileName;
  auto result = readGraph<ReductionGraph<int, int>>(nonExistentFileName);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kNotFound);
}

TEST(ReadGraphTest, heuristicGraph) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 3 6\n"
      "1 4\n"
      "1 5\n"
      "1 6\n"
      "2 4\n"
      "2 6\n"
      "3 6\n");
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};

  auto result = readGraph<HeuristicGraph<int, int>>(testData);

  EXPECT_EQ(result.status(), absl::OkStatus());

  auto myGraph = std::move(result.value());
  myGraph->computeCrossingSums();
  EXPECT_EQ(myGraph->getFreeNodesSize(), 3);
  EXPECT_EQ(myGraph->getFixedNodesSize(), 3);
  EXPECT_EQ(myGraph->getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(0)[0], 0);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(0)[1], 1);
  EXPECT_EQ(myGraph->getFreeNodeNeighboursSize(1), 1);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(1)[0], 0);
  EXPECT_EQ(myGraph->getFreeNodeNeighboursSize(2), 3);

  EXPECT_EQ(myGraph->getFreeNodeNeighbours(2)[0], 0);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(2)[1], 1);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(2)[2], 2);

  EXPECT_EQ(myGraph->getFreeNodeNeighboursSize(0), 2);
  EXPECT_EQ(myGraph->getFreeNodeNeighbours(0)[1], 1);

  EXPECT_EQ(myGraph->getLeftCrossings(0), 0);
  EXPECT_EQ(myGraph->getRightCrossings(0), 2);
  EXPECT_EQ(myGraph->getLeftCrossings(1), 1);
  EXPECT_EQ(myGraph->getRightCrossings(1), 0);
  EXPECT_EQ(myGraph->getLeftCrossings(2), 1);
  EXPECT_EQ(myGraph->getRightCrossings(2), 0);
}

TEST(ReadGraphTest, reductionGraph) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 3 5\n"
      "1 4\n"
      "1 5\n"
      "1 6\n"
      "2 4\n"
      "2 6\n"
      "3 6\n");
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};

  auto result = readGraph<ReductionGraph<int, int>>(testData);

  EXPECT_EQ(result.status(), absl::OkStatus());

  auto graph = std::move(result.value());

  auto [crossingSum, orderVector] =
      reductionalgorithms::algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(*graph);

  std::ofstream outputFile("results.gr");

  // Check if the file is open
  if (!outputFile.is_open()) {
    std::cerr << "Error opening file: results.gr" << std::endl;
  }

  graph->writeResultsToFile(outputFile, orderVector);
  EXPECT_EQ(crossingSum, 1);
  EXPECT_EQ(orderVector[0], 1);
  EXPECT_EQ(orderVector[1], 0);
  EXPECT_EQ(orderVector[2], 2);

  // Check n0, n1, and m
  std::cout << "the number of FreeNodes is " << graph->getFreeNodesSize() << std::endl;
  std::cout << "the number of FixedNodes is " << graph->getFixedNodesSize() << std::endl;
  // ASSERT_EQ(graph->getFixedNodesSize(), 3);  // We have 3 fixed Nodes
  // ASSERT_EQ(graph->getFreeNodesSize(), 4);   // We have 4 free Nodes
}