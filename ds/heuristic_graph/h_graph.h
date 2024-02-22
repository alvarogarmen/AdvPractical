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
  // hold free node id in its corrent position in the pemutation
  std::vector<NodeType> permutation;

 public:
  HGraph(const std::vector<std::vector<NodeType>>& freeNodes,
         const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        freeNodesPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightCrossingSum(std::vector<std::array<CrossingCountType, 2>>(freeNodes.size())),
        permutation(std::vector<NodeType>(freeNodes.size())) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }

  HGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        fixedNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        freeNodesPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightCrossingSum(std::vector<std::array<CrossingCountType, 2>>(freeNodes.size())),
        permutation(std::vector<NodeType>(freeNodes.size())) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }

  NodeType getFixedNodesSize() const { return fixedNodes.size(); }

  NodeType getFixedNodeNeighboursSize(NodeType nodeID) const { return fixedNodes[nodeID].size(); }

  const auto& getFixedNodeNeighbours(NodeType fixedNodeID) const { return fixedNodes[fixedNodeID]; }

  NodeType getFreeNodesSize() const { return freeNodes.size(); }

  NodeType getFreeNodeNeighboursSize(NodeType nodeID) const { return freeNodes[nodeID].size(); }

  const auto& getFreeNodeNeighbours(NodeType freeNodeID) const { return fixedNodes[freeNodeID]; }

  CrossingCountType getLeftCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][0];
  }

  CrossingCountType getRightCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][1];
  }

  const auto& getPermutation() const { return permutation; }

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
    std::swap(permutation[freeNodesPosition[u]], permutation[freeNodesPosition[v]]);
    std::swap(freeNodesPosition[u], freeNodesPosition[v]);
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