#include "ocsm_graph_reader.h"

#include <gtest/gtest.h>

#include "ds/fpt/fpt_graph.h"

TEST(ReadGraphTest, ValidGraph) {
  std::stringstream testData;
  // Create a test graph with stringsteam
  // Edges are not sorted here. If the actual inputs are sorted, then we have to change this but we
  // don´t know yer
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr 3 4 5" << std::endl;
  testData << "1 5" << std::endl;
  testData << "2 6" << std::endl;
  testData << "3 7" << std::endl;
  testData << "4 7" << std::endl;

  auto result = readGraph<FptGraph<int>>(testData);
  EXPECT_EQ(result.status(), absl::OkStatus());
  std::cout << " Pass " << std::endl;

  auto graph = std::move(result.value());

  // Check n0, n1, and m
  ASSERT_EQ(graph->getFixedNodesSize(), 3);
  ASSERT_EQ(graph->getFreeNodesSize(), 4);
  // ASSERT_EQ(graph->getEdgesSize(), 5);   No getEdges function in fptGraph
}

// Test case for invalid header format
TEST(ReadGraphTest, InvalidHeaderFormat) {
  std::stringstream testData;
  // Invalid header format: missing 'p' or incorrect format
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "invalid_header" << std::endl;

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for invalid line format
// Removed because we won´t get such files
/*TEST(ReadGraphTest, InvalidLineFormat) {
  std::stringstream testData;
  // Invalid line format: characters instead of integers
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr 3 4 5" << std::endl;
  testData << "invalid line" << std::endl;

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}
*/
// Test case for nodes out of bounds
TEST(ReadGraphTest, NodesOutOfBounds) {
  std::stringstream testData;
  // Nodes out of bounds: source node greater than n0
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr 3 4 5" << std::endl;
  testData << "4 2" << std::endl;  // Source node is out of bounds

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for negative number of nodes or edges
TEST(ReadGraphTest, NegativeNumbers) {
  std::stringstream testData;
  // Negative number of nodes
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr -3 4 5" << std::endl;  // Negative n0

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for missing 'p' line
TEST(ReadGraphTest, MissingPLine) {
  std::stringstream testData;
  // Missing 'p' line
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "1 2" << std::endl;  // Edge line without 'p' line

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case for file with invalid content
TEST(ReadGraphTest, InvalidFileContent) {
  std::stringstream testData;
  // Invalid content: non-integer characters instead of number of nodes/edges
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr abc 4 5" << std::endl;  // Invalid n0

  auto result = readGraph<FptGraph<int>>(testData);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kInvalidArgument);
}

// Test case to check if the function handles errors correctly.
TEST(ReadGraphTest, MissingInput) {
  // Test with a non-existing file/empty stringstream
  std::stringstream nonExistentFileName;
  auto result = readGraph<FptGraph<int>>(nonExistentFileName);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kNotFound);
}