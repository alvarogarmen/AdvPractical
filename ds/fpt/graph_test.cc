#include <vector>

#include "fpt_graph.h"
#include "fpt_node.h"
#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(GraphTest, SimpleTest) {
  std::vector<std::vector<int>> freeNodes = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<std::vector<int>> fixedNodes = {{0, 1, 2}, {0, 2}, {2}};
  FptGraph myGraph = FptGraph<int>(freeNodes, fixedNodes);

  EXPECT_EQ(myGraph.freeNodes.size(), 3);
  EXPECT_EQ(myGraph.fixedNodes.size(), 3);
  EXPECT_EQ(myGraph.freeNodes[0].size(), 2);
  EXPECT_EQ(myGraph.freeNodes[0][1], 1);

  myGraph.buildYx();
  // EXPECT_EQ(myGraph.yx[0].size(), 0);
  // EXPECT_EQ(myGraph.yx[1].size(), 1);
  // EXPECT_EQ(myGraph.yx[2].size(), 0);
  myGraph.fillCrossingMatrix();
  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);
  EXPECT_EQ(myGraph.getCrossing(0, 2), 1);
  EXPECT_EQ(myGraph.getCrossing(1, 2), 0);
  EXPECT_EQ(myGraph.getCrossing(2, 0), 3);
  EXPECT_EQ(myGraph.getCrossing(2, 1), 2);
  EXPECT_EQ(myGraph.getCrossing(1, 0), 0);
  EXPECT_EQ(myGraph.getCrossing(2, 0), 3);

  // EXPECT_EQ((myGraph.getXNode(0)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getXNode(1)).yx.size(), 1);
  // EXPECT_EQ((myGraph.getXNode(2)).yx.size(), 0);
  // EXPECT_EQ((myGraph.getXNode(1)).yx[0], 2);
  //  EXPECT_EQ(myGraph.getCrossing(0, 1), 1);
}

TEST(GraphTest, FreeNodeTest) {
  std::vector<std::vector<int>> freeNodeData = {{0, 1}, {0}, {0, 1, 2}};
  std::vector<FreeNode<int>> freeNodes;
  for (auto i = 0; i < freeNodeData.size(); i++) {
    FreeNode free = FreeNode(i, freeNodeData[i], i);
    freeNodes.push_back(free);
    EXPECT_EQ(freeNodes[i].neighbours.size(), freeNodeData[i].size());
    EXPECT_EQ(freeNodes[i].neighbours, freeNodeData[i]);
  }
}