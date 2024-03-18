#include "oscm/grader.h"

#include "ds/bipartite_graph.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraderTest, SimpleTest) {
  BipartiteGraph myGraph = BipartiteGraph(4, 4, 4);

  myGraph.addEdge(0, 2);
  myGraph.addEdge(1, 1);
  myGraph.addEdge(2, 0);
  myGraph.addEdge(3, 2);
  int crossings = crossGrader(myGraph);
  EXPECT_EQ(crossings, 3);
}