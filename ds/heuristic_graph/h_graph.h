#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

template <typename NodeType, typename CrossingCountType>
class HGraph {
  // for each free node holds its neighbours
  std::vector<std::vector<NodeType>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NodeType>> fixedNodes;
  // free nodes current position in the permutation
  std::vector<NodeType> freeNodesPosition;
  // for each free node u, the first place hold the crossing sum with free node positioned to the
  // left of u and the second hold the sum to its right
  std::vector<std::array<CrossingCountType, 2>> leftRightCrossingSum;

 public:
  HGraph(const std::vector<std::vector<NodeType>>& freeNodes,
         const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        freeNodesPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightCrossingSum(std::vector<std::array<CrossingCountType, 2>>(freeNodes.size())) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    computeCrossingSums();
  }

  HGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        fixedNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        freeNodesPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightCrossingSum(std::vector<std::array<CrossingCountType, 2>>(freeNodes.size())) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    computeCrossingSums();
  }

  NodeType const getFixedNodesSize() { return fixedNodes.size(); }

  NodeType const getFixedNodeNeighboursSize(NodeType nodeID) { return fixedNodes[nodeID].size(); }

  auto& getFixedNodeNeighbours(NodeType fixedNodeID) { return fixedNodes[fixedNodeID]; }

  NodeType const getFreeNodesSize() { return freeNodes.size(); }

  NodeType const getFreeNodeNeighboursSize(NodeType nodeID) { return freeNodes[nodeID].size(); }

  auto const getFreeNodeNeighbours(NodeType freeNodeID) { return fixedNodes[freeNodeID]; }

  CrossingCountType const getLeftCrossings(NodeType freeNodeID) {
    return leftRightCrossingSum[freeNodeID][0];
  }

  CrossingCountType const getRightCrossings(NodeType freeNodeID) {
    return leftRightCrossingSum[freeNodeID][1];
  }

  void addEdge(NodeType freeNode, NodeType fixedNode) {
    freeNodes[freeNode].push_back(fixedNode);
    fixedNodes[fixedNode].push_back(freeNode);
  }

  // copmute the number of crossings created by the edges from two free nodes (u, v)
  // when u is to the left of v
  CrossingCountType const computeUVcrossing(NodeType u, NodeType v) {
    NodeType crossingSum = 0;
    for (const auto uNeighbour : freeNodes[u]) {
      for (const auto vNeighbour : freeNodes[v]) {
        crossingSum += vNeighbour < uNeighbour;
      }
    }
    return crossingSum;
  }

  // for two free nodes u, v switch there positions and update left and right crossings
  // assume u is the left neighbour of v
  void switchNeighbours(NodeType u, NodeType v) {
    // the number of crossings created by the edges from (u, v)
    NodeType uvSum = computeUVcrossing(u, v);
    // the number of crossings created by the edges from (v, u)
    NodeType vuSum = computeUVcrossing(v, u);
    leftRightCrossingSum[u][1] -= uvSum;
    leftRightCrossingSum[v][1] += vuSum;
    leftRightCrossingSum[u][0] += vuSum;
    leftRightCrossingSum[v][0] -= uvSum;
    ++freeNodesPosition[u];
    --freeNodesPosition[v];
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right
  void computeCrossingSums() {
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      for (NodeType j = i + 1; j < freeNodes.size(); ++j) {
        if (freeNodesPosition[i] > freeNodesPosition[j]) {
          CrossingCountType crossing = computeUVcrossing(j, i);
          leftRightCrossingSum[i][0] += crossing;
          leftRightCrossingSum[j][1] += crossing;
        } else {
          CrossingCountType crossing = computeUVcrossing(i, j);
          leftRightCrossingSum[i][1] += crossing;
          leftRightCrossingSum[j][0] += crossing;
        }
      }
    }
  }
};