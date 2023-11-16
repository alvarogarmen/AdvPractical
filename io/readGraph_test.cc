#include "readGraph.h"

#include <gtest/gtest.h>

#include "ds/fpt/fpt_graph.h"
TEST(ReadGraphTest, ValidGraph) {
  std::stringstream testData;
  // Create a test graph with stringsteam
  // Edges are not sorted here. If the actual inputs are sorted, then we have to change this but we
  // don´t know yer
  testData << "c this is a comment and should not be read" << std::endl;
  testData << "p ocr 3 4 5" << std::endl;
  testData << "1 4" << std::endl;
  testData << "2 3" << std::endl;
  testData << "2 4" << std::endl;
  testData << "3 3" << std::endl;
  testData << "1 2" << std::endl;

  auto result = readGraph<BipartiteGraph>(testData);
  ASSERT_TRUE(result.ok());
  auto graph = std::move(result.value());

  // Check n0, n1, and m
  ASSERT_EQ(graph->getNumNodes0(), 3);
  ASSERT_EQ(graph->getNumNodes1(), 4);
  ASSERT_EQ(graph->getNumEdges(), 5);
}

// Test case to check if the function handles errors correctly.
TEST(ReadGraphTest, ErrorHandling) {
  // Test with a non-existing file. Please, never write a non_existent.txt file
  const std::string nonExistentFileName = "non_existent.txt";
  auto result = readGraph<BipartiteGraph>(nonExistentFileName);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kNotFound);
}