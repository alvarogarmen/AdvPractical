#include "oscm/grader.h"

#include "ds/bipartite_graph.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraderTest, SimpleTest) {
  BipartiteGraph myGraph = BipartiteGraph(4, 4, 0);

  myGraph.addEdge(4, 3);
  myGraph.addEdge(5, 2);
  myGraph.addEdge(6, 1);
  myGraph.addEdge(7, 3);

  EXPECT_EQ(crossGrader(myGraph), 3);
}