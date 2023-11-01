#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "fpt_node.h"

template <typename NodeType>
class FptGraph {
 public:
  FptGraph(std::vector<NonFixedNode<NodeType>*> NonFixedNodes,
           std::vector<FixedNode<NodeType>*> FixedNodes)
      : NonFixedNodes(NonFixedNodes), FixedNodes(FixedNodes) {
    crossingMatrix = std::vector<std::vector<NodeType>>(
        NonFixedNodes.size(), std::vector<NodeType>(NonFixedNodes.size()));
  }
  FptGraph(NodeType NonFixedNodesSize) {
    crossingMatrix = std::vector<std::vector<NodeType>>(NonFixedNodesSize,
                                                        std::vector<NodeType>(NonFixedNodesSize));
  }

  void insertNonFixedNode(NonFixedNode<NodeType> y);
  void insertFixedNode(FixedNode<NodeType> x);
  void adjustCrossingMatrix(NodeType uIndex, NodeType vIndex, NodeType sumOfCrossing);
  NonFixedNode<NodeType> getNonFixedNode(NodeType index);
  FixedNode<NodeType> getFixedNode(NodeType index);
  NodeType getCrossing(NodeType uIndex, NodeType vIndex);
  NodeType getNonFixedNodesSize();
  NodeType getFixedNodesSize();
  void buildYx();
  void initCrossingMatrix();
  void fillCrossingMatrix();

 private:
  std::vector<NonFixedNode<NodeType>> NonFixedNodes;
  std::vector<FixedNode<NodeType>> FixedNodes;
  std::vector<std::vector<NodeType>> crossingMatrix;
};

template <typename NodeType>
void FptGraph<NodeType>::insertNonFixedNode(NonFixedNode<NodeType> y) {
  NonFixedNodes.push_back(y);
}

template <typename NodeType>
void FptGraph<NodeType>::insertFixedNode(FixedNode<NodeType> x) {
  FixedNodes.push_back(x);
}

template <typename NodeType>
void FptGraph<NodeType>::adjustCrossingMatrix(NodeType uIndex, NodeType vIndex,
                                              NodeType sumOfCrossing) {
  crossingMatrix[uIndex][vIndex] = sumOfCrossing;
}

template <typename NodeType>
NonFixedNode<NodeType> FptGraph<NodeType>::getNonFixedNode(NodeType index) {
  return NonFixedNodes[index];
}

template <typename NodeType>
FixedNode<NodeType> FptGraph<NodeType>::getFixedNode(NodeType index) {
  return FixedNodes[index];
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getCrossing(NodeType uIndex, NodeType vIndex) {
  return crossingMatrix[uIndex][vIndex];
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getNonFixedNodesSize() {
  return NonFixedNodes.size();
}

template <typename NodeType>
NodeType FptGraph<NodeType>::getFixedNodesSize() {
  return FixedNodes.size();
}

/*
template <typename NodeType>
void FptGraph<NodeType>::buildYx() {
  std::for_each(NonFixedNodes.begin(), NonFixedNodes.end(), [](NonFixedNode<NodeType>& y) {
    std::for_each(y.neighbours.begin() + 1, y.neighbours.end() - 1, [](NodeType xNeighbourI) {
      std::cout << "xNeighbourI: " << xNeighbourI << std::endl;
    });
  });
}
*/
template <typename NodeType>
void FptGraph<NodeType>::buildYx() {
  for (NonFixedNode<NodeType> y : NonFixedNodes) {
    for (NodeType i = 1; i < y.neighbours.size(); ++i) {  // y.neighbours.size() -1
      FixedNodes[y.neighbours[i]].yx.push_back(y.getNodeID());
    }
  }
}

template <typename NodeType>
void FptGraph<NodeType>::initCrossingMatrix() {
  for (FixedNode<NodeType> x : FixedNodes) {
    for (NodeType neighbourI = 0; neighbourI < x.neighbours.size(); ++neighbourI) {
      for (NodeType yxI = 0; yxI < x.yx.size(); ++yxI) {
        crossingMatrix[x.neighbours[neighbourI].getNodeID()][x.yx[yxI]] = 0;
      }
    }
  }
}

template <typename NodeType>
void FptGraph<NodeType>::fillCrossingMatrix() {
  for (FixedNode<NodeType> x : FixedNodes) {
    for (NodeType neighbourI = 0; neighbourI < x.neighbours.size(); ++neighbourI) {
      for (NodeType yxI = 0; yxI < x.yx.size(); ++yxI) {
        NodeType degreeSThenX;
        for (degreeSThenX = NonFixedNodes[x.yx[yxI]].getDxScannedIndex();
             degreeSThenX < x.getNodeID(); ++degreeSThenX) {
        }
        crossingMatrix[x.neighbours[neighbourI].getNodeID()][x.yx[yxI]] += degreeSThenX;
        NonFixedNodes[x.yx[yxI]].setDxScannedIndex(degreeSThenX);
        if (x.neighbours[neighbourI].neighbours[x.neighbours[neighbourI].neighbours.size() - 1] ==
            x.getNodeID()) {
          crossingMatrix[x.yx[yxI]][x.neighbours[neighbourI].getNodeID()] +=
              x.neighbours[neighbourI].neighbours.size() *
              (NonFixedNodes[x.yx[yxI]].neighbours.size() - degreeSThenX);
        }
      }
    }
  }
}