#include "oscm/grader.h"

#include "ds/bipartite_graph.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraderTest, SimpleTest) {
  Graph myGraph = Graph(0, 0, 0);

  Node node1 = Node(1, 2);
  Node node2 = Node(2, 1);
  Node node3 = Node(3, 1);
  Node node4 = Node(4, 1);

  myGraph.insertLeftNode(node1);
  myGraph.insertLeftNode(node2);
  myGraph.insertLeftNode(node3);
  myGraph.insertLeftNode(node4);

  myGraph.insertEdge(1, 1);
  myGraph.insertEdge(1, 2);
  myGraph.insertEdge(2, 2);
  myGraph.insertEdge(3, 1);
  myGraph.insertEdge(4, 3);

  EXPECT_EQ(crossGrader(myGraph), 2);
}