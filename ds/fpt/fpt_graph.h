#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "fpt_node.h"

template <typename NodeType>
class FptGraph {
 public:
  FptGraph(std::vector<YNode<NodeType>*> yNodes, std::vector<XNode<NodeType>*> xNodes)
      : yNodes(yNodes), xNodes(xNodes) {
    crossingMatrix =
        std::vector<std::vector<NodeType>>(yNodes.size(), std::vector<NodeType>(yNodes.size()));
  }
  FptGraph(NodeType yNodesSize) {
    crossingMatrix =
        std::vector<std::vector<NodeType>>(yNodesSize, std::vector<NodeType>(yNodesSize));
  }

  void insertYNode(YNode<NodeType> y);
  void insertXNode(XNode<NodeType> x);
  void adjustCrossingMatrix(NodeType uIndex, NodeType vIndex, NodeType sumOfCrossing);
  YNode<NodeType> getYNode(NodeType index);
  XNode<NodeType> getXNode(NodeType index);
  NodeType getCrossing(NodeType uIndex, NodeType vIndex);
  NodeType getYNodesSize();
  NodeType getXNodesSize();
  void buildYx();
  void initCrossingMatrix();
  void fillCrossingMatrix();

 private:
  std::vector<YNode<NodeType>> yNodes;
  std::vector<XNode<NodeType>> xNodes;
  std::vector<std::vector<NodeType>> crossingMatrix;
};

template <typename NodeType>
void FptGraph<NodeType>::insertYNode(YNode<NodeType> y) {
  yNodes.push_back(y);
}

template <typename NodeType>
void FptGraph<NodeType>::insertXNode(XNode<NodeType> x) {
  xNodes.push_back(x);
}

template <typename NodeType>
void FptGraph<NodeType>::adjustCrossingMatrix(NodeType uIndex, NodeType vIndex,
                                              NodeType sumOfCrossing) {
  crossingMatrix[uIndex][vIndex] = sumOfCrossing;
}

template <typename NodeType>
YNode<NodeType> FptGraph<NodeType>::getYNode(NodeType index) {
  return yNodes[index];
}

template <typename NodeType>
XNode<NodeType> FptGraph<NodeType>::getXNode(NodeType index) {
  return xNodes[index];
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getCrossing(NodeType uIndex, NodeType vIndex) {
  return crossingMatrix[uIndex][vIndex];
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getYNodesSize() {
  return yNodes.size();
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getXNodesSize() {
  return xNodes.size();
}

/*
template <typename NodeType>
void FptGraph<NodeType>::buildYx() {
  std::for_each(yNodes.begin(), yNodes.end(), [](YNode<NodeType>& y) {
    std::for_each(y.neighbours.begin() + 1, y.neighbours.end() - 1, [](NodeType xNeighbourI) {
      std::cout << "xNeighbourI: " << xNeighbourI << std::endl;
    });
  });
}
*/
template <typename NodeType>
void FptGraph<NodeType>::buildYx() {
  for (YNode<NodeType> y : yNodes) {
    for (NodeType i = 1; i < y.neighbours.size(); ++i) {  // y.neighbours.size() -1
      xNodes[y.neighbours[i]].yx.push_back(y.getNodeID());
    }
  }
}

template <typename NodeType>
void FptGraph<NodeType>::initCrossingMatrix() {
  for (XNode<NodeType> x : xNodes) {
    for (NodeType neighbourI = 0; neighbourI < x.neighbours.size(); ++neighbourI) {
      for (NodeType yxI = 0; yxI < x.yx.size(); ++yxI) {
        crossingMatrix[x.neighbours[neighbourI].getNodeID()][x.yx[yxI]] = 0;
      }
    }
  }
}

template <typename NodeType>
void FptGraph<NodeType>::fillCrossingMatrix() {
  for (XNode<NodeType> x : xNodes) {
    for (NodeType neighbourI = 0; neighbourI < x.neighbours.size(); ++neighbourI) {
      for (NodeType yxI = 0; yxI < x.yx.size(); ++yxI) {
        NodeType degreeSThenX;
        for (degreeSThenX = yNodes[x.yx[yxI]].getDxScannedIndex(); degreeSThenX < x.getNodeID();
             ++degreeSThenX) {
        }
        crossingMatrix[x.neighbours[neighbourI].getNodeID()][x.yx[yxI]] += degreeSThenX;
        yNodes[x.yx[yxI]].setDxScannedIndex(degreeSThenX);
        if (x.neighbours[neighbourI].neighbours[x.neighbours[neighbourI].neighbours.size() - 1] ==
            x.getNodeID()) {
          crossingMatrix[x.yx[yxI]][x.neighbours[neighbourI].getNodeID()] +=
              x.neighbours[neighbourI].neighbours.size() *
              (yNodes[x.yx[yxI]].neighbours.size() - degreeSThenX);
        }
      }
    }
  }
}