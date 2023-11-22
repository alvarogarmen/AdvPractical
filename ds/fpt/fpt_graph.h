#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

template <typename NT>
class FptGraph {
 public:
  using NodeType = NT;
  FptGraph(const std::vector<std::vector<NodeType>>& freeNodes,
           const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes), fixedNodes(fixedNodes) {
    crossingMatrix = std::vector<std::vector<NodeType>>(freeNodes.size(),
                                                        std::vector<NodeType>(freeNodes.size(), 0));
    yx = std::vector<std::vector<NodeType>>(freeNodes.size(), std::vector<NodeType>(0));
    dxScannedIndex = std::vector<NodeType>(freeNodes.size());
    ;
  }

  FptGraph(NodeType numFreeNodes, NodeType numFixedNodes, NodeType numEdges) {
    freeNodes = std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0));
    fixedNodes = std::vector<std::vector<NodeType>>(numFixedNodes, std::vector<NodeType>(0));
    crossingMatrix =
        std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(numFreeNodes, 0));
    yx = std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0));
    dxScannedIndex = std::vector<NodeType>(numFreeNodes);
  }

  void adjustCrossingMatrix(NodeType uIndex, NodeType vIndex, NodeType sumOfCrossing) {
    crossingMatrix[uIndex][vIndex] = sumOfCrossing;
  }

  // return the number of crossings created by the edges from two free nodes (u, v)
  NodeType const getCrossing(NodeType const uIndex, NodeType const vIndex) {
    return crossingMatrix[uIndex][vIndex];
  }

  NodeType const getFixedNodesSize() { return fixedNodes.size(); }

  NodeType const getFixedNodeNeighboursSize(NodeType nodeID) { return fixedNodes[nodeID].size(); }

  NodeType const getFixedNodeNeighbour(NodeType fixedNodeID, NodeType neighbourI) {
    return fixedNodes[fixedNodeID][neighbourI];
  }

  NodeType const getFreeNodesSize() { return freeNodes.size(); }

  NodeType const getFreeNodeNeighboursSize(NodeType nodeID) { return freeNodes[nodeID].size(); }

  NodeType const getFreeNodeNeighbour(NodeType freeNodeID, NodeType neighbourI) {
    return fixedNodes[freeNodeID][neighbourI];
  }

  void addEdge(NodeType freeNode, NodeType fixedNode) {
    freeNodes[freeNode].push_back(fixedNode);
    fixedNodes[fixedNode].push_back(freeNode);
  }

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
        for (const auto freeNodeInYx : yx[fixedNodeI]) {
          for (; freeNodes[freeNodeInYx][dxScannedIndex[freeNodeInYx]] < fixedNodeI;
               ++dxScannedIndex[freeNodeInYx]) {
          }
          crossingMatrix[fixedNodes[fixedNodeI][neighbourI]][freeNodeInYx] +=
              dxScannedIndex[freeNodeInYx];
          if (freeNodes[neighbourI][freeNodes[neighbourI].size() - 1] == fixedNodeI) {
            NodeType u = freeNodeInYx;
            NodeType v = fixedNodes[fixedNodeI][neighbourI];
            crossingMatrix[u][v] +=
                freeNodes[v].size() * (freeNodes[u].size() - (dxScannedIndex[u] + 1));
          }
        }
        if (freeNodes[neighbourI][freeNodes[neighbourI].size() - 1] == fixedNodeI) {
          for (const auto u : fixedNodes[fixedNodeI]) {
            // if it is the most left then do if its the most right (d(u)−d≤x(u)) would be 0
            if (freeNodes[u][0] == fixedNodeI) {
              NodeType v = fixedNodes[fixedNodeI][neighbourI];
              crossingMatrix[u][v] += freeNodes[v].size() * (freeNodes[u].size() - 1);
            }
          }
        }
      }
    }
  }

 private:
  // save the crossing number that accur between two free nodes
  // (y1, y2) saves the crossing number of y1 and y2 assuming y1 is placed before y2
  std::vector<std::vector<NodeType>> crossingMatrix;
  // for each fixed node x, save all the free nodes y, such that x is between their most left and
  // right neighbours
  std::vector<std::vector<NodeType>> yx;
  // for each free node holds its neighbours
  std::vector<std::vector<NodeType>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NodeType>> fixedNodes;

  // dxScannedIndex is the index of the last scanned x
  // used to count amount of the node neighbours smaller than x
  std::vector<NodeType> dxScannedIndex;
};