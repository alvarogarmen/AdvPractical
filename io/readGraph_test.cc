#include "readGraph.h"

#include <gtest/gtest.h>

#include "ds/fpt/fpt_graph.h"

// Test case to check if the function reads a valid graph from a file.
TEST(ReadGraphTest, ValidGraph) {
  const std::string testFileName = "test_valid_graph.txt";
  // Create a test graph in file "test_valid_graph.txt"
  // Edges are not sorted here. If the actual inputs are sorted, then we have to change this but we
  // don´t know yer
  std::ofstream testFile(testFileName);
  testFile << "c this is a comment and should not be read" << std::endl;
  testFile << "p ocr 3 4 5" << std::endl;
  testFile << "1 4" << std::endl;
  testFile << "2 3" << std::endl;
  testFile << "2 4" << std::endl;
  testFile << "3 3" << std::endl;
  testFile << "1 2" << std::endl;
  testFile.close();

  auto result = readGraph<BipartiteGraph>(testFileName);
  ASSERT_TRUE(result.ok());
  auto graph = std::move(result.value());

  // Check n0, n1 and m
  ASSERT_EQ(graph->getNumNodes0(), 3);
  ASSERT_EQ(graph->getNumNodes1(), 4);
  ASSERT_EQ(graph->getNumEdges(), 5);
  std::remove(testFileName.c_str());
}

// Test case to check if the function handles errors correctly.
TEST(ReadGraphTest, ErrorHandling) {
  // Test with a non-existing file. Please, never write a non_existent.txt file
  const std::string nonExistentFileName = "non_existent.txt";
  auto result = readGraph<BipartiteGraph>(nonExistentFileName);
  ASSERT_FALSE(result.ok());
  ASSERT_EQ(result.status().code(), absl::StatusCode::kNotFound);
}