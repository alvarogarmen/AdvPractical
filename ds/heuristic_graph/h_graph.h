#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

template <typename NodeType>
class HGraph {
 public:
  HGraph(const std::vector<std::vector<NodeType>>& freeNodes,
         const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes), fixedNodes(fixedNodes) {
    rightCrossingSum = std::vector<NodeType>(freeNodes.size(), 0);
    leftCrossingSum = std::vector<NodeType>(freeNodes.size(), 0);

    freeNodesPosition = std::vector<NodeType>(freeNodes.size());
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      freeNodesPosition[i] = i;
    }
  }

  HGraph(NodeType numFreeNodes, NodeType numFixedNodes) {
    freeNodes = std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0));
    fixedNodes = std::vector<std::vector<NodeType>>(numFixedNodes, std::vector<NodeType>(0));
    rightCrossingSum = std::vector<NodeType>(freeNodes.size(), 0);
    leftCrossingSum = std::vector<NodeType>(freeNodes.size(), 0);
    freeNodesPosition = std::vector<NodeType>(freeNodes.size());
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      freeNodesPosition[i] = i;
    }
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

  NodeType const getLeftCrossings(NodeType freeNodeID) { return leftCrossingSum[freeNodeID]; }

  NodeType const getRightCrossings(NodeType freeNodeID) { return rightCrossingSum[freeNodeID]; }

  void addEdge(NodeType freeNode, NodeType fixedNode) {
    freeNodes[freeNode].push_back(fixedNode);
    fixedNodes[fixedNode].push_back(freeNode);
  }

  // copmute the number of crossings created by the edges from two free nodes (u, v)
  NodeType const computeUVcrossing(NodeType u, NodeType v) {
    NodeType crossingSum = 0;
    for (const auto uNeighbour : freeNodes[u]) {
      for (const auto vNeighbour : freeNodes[v]) {
        if (vNeighbour < uNeighbour) {
          ++crossingSum;
        }
      }
    }
    return crossingSum;
  }

  // for two free nodes u, v switch there positions and update left and right crossings
  void switchNeighbours(NodeType u, NodeType v) {
    // the number of crossings created by the edges from (u, v)
    NodeType uvSum = computeUVcrossing(u, v);
    // the number of crossings created by the edges from (v, u)
    NodeType vuSum = computeUVcrossing(v, u);
    rightCrossingSum[u] -= uvSum;
    rightCrossingSum[v] += vuSum;
    leftCrossingSum[u] += vuSum;
    leftCrossingSum[v] -= uvSum;
    ++freeNodesPosition[u];
    --freeNodesPosition[v];
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right
  void computeCrossingSums() {
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      for (NodeType j = 0; j < freeNodes.size(); ++j) {
        if (freeNodesPosition[i] > freeNodesPosition[j]) {
          leftCrossingSum[i] += computeUVcrossing(j, i);
        } else if (freeNodesPosition[i] < freeNodesPosition[j]) {
          rightCrossingSum[i] += computeUVcrossing(i, j);
        }
      }
    }
  }

 private:
  // for each free node holds its neighbours
  std::vector<std::vector<NodeType>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NodeType>> fixedNodes;
  // free nodes current position in the permutation
  std::vector<NodeType> freeNodesPosition;
  // for each free node u, hold the crossing sum with free node positioned to the left of u
  std::vector<NodeType> leftCrossingSum;
  // for each free node u, hold the crossing sum with free node positioned to the right of u
  std::vector<NodeType> rightCrossingSum;
};