#include <vector>

#include "fpt_graph.h"
#include "fpt_node.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraphTest, SimpleTest) {
  FptGraph myGraph = FptGraph<int>(3);
  std::vector<int> neighboursY0 = {0, 1};
  std::vector<int> neighboursY1 = {0, 2};  //{0}
  std::vector<int> neighboursY2 = {0, 1, 2};
  NonFixedNode y0 = NonFixedNode(0, neighboursY0, 0);
  NonFixedNode y1 = NonFixedNode(1, neighboursY1, 1);
  NonFixedNode y2 = NonFixedNode(2, neighboursY2, 2);
  std::vector<NonFixedNode<int>> neighboursX0 = {y0, y1, y2};
  std::vector<NonFixedNode<int>> neighboursX1 = {y0, y2};
  std::vector<NonFixedNode<int>> neighboursX2 = {y1, y2};  //{y2}
  FixedNode<int> x0 = FixedNode(0, neighboursX0);
  FixedNode<int> x1 = FixedNode(1, neighboursX1);
  FixedNode<int> x2 = FixedNode(1, neighboursX2);
  myGraph.insertFixedNode(x0);
  myGraph.insertFixedNode(x1);
  myGraph.insertFixedNode(x2);
  myGraph.insertNonFixedNode(y0);
  myGraph.insertNonFixedNode(y1);
  myGraph.insertNonFixedNode(y2);

  EXPECT_EQ(myGraph.getFixedNodesSize(), 3);
  myGraph.buildYx();
  myGraph.initCrossingMatrix();
  EXPECT_EQ(myGraph.getCrossing(0, 1), 0);
  myGraph.fillCrossingMatrix();
  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);

  // EXPECT_EQ((myGraph.getFixedNode(0)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getFixedNode(1)).yx.size(), 1);
  // EXPECT_EQ((myGraph.getFixedNode(2)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getFixedNode(1)).yx[0], 2);
  //  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);
}