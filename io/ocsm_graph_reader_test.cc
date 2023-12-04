#include "ocsm_graph_reader.h"

#include <gtest/gtest.h>

#include "ds/bipartite_graph.h"
#include "gmock/gmock-matchers.h"

namespace {
using ::testing::ElementsAre;
}

TEST(ReadGraphTest, ValidGraph) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 4 5\n"
      "1 5\n"
      "2 6\n"
      "3 7\n"
      "4 7\n");

  auto result = readGraph<BipartiteGraph<int>>(testData);
  EXPECT_EQ(result.status(), absl::OkStatus());

  auto graph = std::move(result.value());

  // Check n0, n1, and m
  EXPECT_EQ(graph->edges[0].source, 0);  // 0 instead of 1 because internally we use 0-indexation
  EXPECT_EQ(graph->edges[0].target, 4);  // instead of 1-indexation
  ASSERT_EQ(graph->getFixedNodesSize(), 3);
  ASSERT_EQ(graph->getFreeNodesSize(), 4);
  ASSERT_EQ(graph->getEdgesSize(), 4);
}

// Test case for invalid header format
TEST(ReadGraphTest, InvalidHeaderFormat) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "invalid header\n");
  // Invalid header format: missing 'p' or incorrect format

  auto result = readGraph<BipartiteGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for invalid line format
// Removed because we won´t get such files
/*TEST(ReadGraphTest, InvalidLineFormat) {
  std::stringstream testData;
  // Invalid line format: characters instead of integers
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr 3 4 5" << std::endl;
  testData << "invalid line" << std::endl;

  auto result = readGraph<BipartiteGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}*/

// Test case for nodes out of bounds
TEST(ReadGraphTest, NodesOutOfBounds) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr 3 4 5\n"
      "4 2\n");
  // Nodes out of bounds: source node greater than n0
  // Source node is out of bounds

  auto result = readGraph<BipartiteGraph<int>>(testData);
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

  auto result = readGraph<BipartiteGraph<int>>(testData);
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for missing 'p' line
TEST(ReadGraphTest, MissingPLine) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "1 2\n");
  // Missing 'p' line
  // Edge line without 'p' line

  auto result = readGraph<BipartiteGraph<int>>(testData);
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for file with invalid content
TEST(ReadGraphTest, InvalidFileContent) {
  std::stringstream testData(
      "c this is a comment and should not be read\n"
      "p ocr abc 4 5\n");
  // Invalid content: non-integer characters instead of number of nodes/edges
  // Invalid n0

  auto result = readGraph<BipartiteGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case to check if the function handles errors correctly.
TEST(ReadGraphTest, MissingInput) {
  // Test with a non-existing file/empty stringstream
  std::stringstream nonExistentFileName;
  auto result = readGraph<BipartiteGraph<int>>(nonExistentFileName);
  ASSERT_FALSE(result.ok());
  EXPECT_EQ(result.status().code(), absl::StatusCode::kNotFound);
}