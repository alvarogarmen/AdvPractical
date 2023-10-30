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
  YNode y0 = YNode(0, neighboursY0, 0);
  YNode y1 = YNode(1, neighboursY1, 1);
  YNode y2 = YNode(2, neighboursY2, 2);
  std::vector<YNode<int>> neighboursX0 = {y0, y1, y2};
  std::vector<YNode<int>> neighboursX1 = {y0, y2};
  std::vector<YNode<int>> neighboursX2 = {y1, y2};  //{y2}
  XNode<int> x0 = XNode(0, neighboursX0);
  XNode<int> x1 = XNode(1, neighboursX1);
  XNode<int> x2 = XNode(1, neighboursX2);
  myGraph.insertXNode(x0);
  myGraph.insertXNode(x1);
  myGraph.insertXNode(x2);
  myGraph.insertYNode(y0);
  myGraph.insertYNode(y1);
  myGraph.insertYNode(y2);

  EXPECT_EQ(myGraph.getXNodesSize(), 3);
  myGraph.buildYx();
  myGraph.initCrossingMatrix();
  EXPECT_EQ(myGraph.getCrossing(0, 1), 0);
  myGraph.fillCrossingMatrix();
  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);

  // EXPECT_EQ((myGraph.getXNode(0)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getXNode(1)).yx.size(), 1);
  // EXPECT_EQ((myGraph.getXNode(2)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getXNode(1)).yx[0], 2);
  //  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);
}