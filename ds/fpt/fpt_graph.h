#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "fpt_node.h"

template <typename NodeType>
class FptGraph {
 public:
  FptGraph(std::vector<std::vector<NodeType>> freeNodes,
           std::vector<std::vector<NodeType>> fixedNodes)
      : freeNodes(freeNodes), fixedNodes(fixedNodes) {
    crossingMatrix = std::vector<std::vector<NodeType>>(freeNodes.size(),
                                                        std::vector<NodeType>(freeNodes.size()));
    yx = std::vector<std::vector<NodeType>>(freeNodes.size(), std::vector<NodeType>(0));
    dxScannedIndex = std::vector<NodeType>(freeNodes.size());
    ;
  }
  FptGraph(NodeType yNodesSize) {
    crossingMatrix =
        std::vector<std::vector<NodeType>>(yNodesSize, std::vector<NodeType>(yNodesSize));
  }

  std::vector<std::vector<NodeType>> yx;
  std::vector<std::vector<NodeType>> freeNodes;
  std::vector<std::vector<NodeType>> fixedNodes;

  void adjustCrossingMatrix(NodeType uIndex, NodeType vIndex, NodeType sumOfCrossing) {
    crossingMatrix[uIndex][vIndex] = sumOfCrossing;
  }
  NodeType getCrossing(NodeType uIndex, NodeType vIndex) { return crossingMatrix[uIndex][vIndex]; }

  void buildYx() {
    for (NodeType freeNodeI = 0; freeNodeI < freeNodes.size(); ++freeNodeI) {
      for (NodeType neighbourI = 1; neighbourI < freeNodes[freeNodeI].size(); ++neighbourI) {
        yx[freeNodes[freeNodeI][neighbourI]].push_back(freeNodeI);
      }
    }
  }

  void fillCrossingMatrix() {
    for (NodeType fixedNodeI = 0; fixedNodeI < fixedNodes.size(); ++fixedNodeI) {
      for (NodeType neighbourI = 0; neighbourI < fixedNodes[fixedNodeI].size(); ++neighbourI) {
        for (NodeType yxI = 0; yxI < yx[fixedNodeI].size(); ++yxI) {
          for (; freeNodes[yx[fixedNodeI][yxI]][dxScannedIndex[yx[fixedNodeI][yxI]]] < fixedNodeI;
               ++dxScannedIndex[yx[fixedNodeI][yxI]]) {
          }
          crossingMatrix[fixedNodes[fixedNodeI][neighbourI]][yx[fixedNodeI][yxI]] +=
              dxScannedIndex[yx[fixedNodeI][yxI]];
          if (freeNodes[neighbourI][freeNodes[neighbourI].size() - 1] == fixedNodeI) {
            NodeType u = yx[fixedNodeI][yxI];
            NodeType v = fixedNodes[fixedNodeI][neighbourI];
            crossingMatrix[u][v] +=
                freeNodes[v].size() * (freeNodes[u].size() - (dxScannedIndex[u] + 1));
          }
        }
        if (freeNodes[neighbourI][freeNodes[neighbourI].size() - 1] == fixedNodeI) {
          for (NodeType neighbour = 0; neighbour < fixedNodes[fixedNodeI].size(); ++neighbour) {
            // if it is the most left then do if its the most right (d(u)−d≤x(u)) would be 0
            if (freeNodes[fixedNodes[fixedNodeI][neighbour]][0] == fixedNodeI) {
              NodeType u = fixedNodes[fixedNodeI][neighbour];
              NodeType v = fixedNodes[fixedNodeI][neighbourI];
              crossingMatrix[u][v] += freeNodes[v].size() * (freeNodes[u].size() - 1);
            }
          }
        }
      }
    }
  }

 private:
  std::vector<std::vector<NodeType>> crossingMatrix;
  // dxScannedIndex is the index of the last scanned x
  // used to count amount of the node neighbours smaller than x
  std::vector<NodeType> dxScannedIndex;
};