//
// Created by alvar on 25/10/2023.
//

#pragma once

#include <unordered_map>
#include <vector>

template <typename NodeType>
struct Edge {
  Edge(NodeType source, NodeType target);
  NodeType source;
  NodeType target;
};

template <typename NodeType>
Edge<NodeType>::Edge(NodeType source, NodeType target) {
  source = source;
  target = target;
}

template <typename NT>
struct BipartiteGraph {
  using NodeType = NT;
  BipartiteGraph() {
    freeNodes = std::vector<NodeType>(0);
    fixedNodes = std::vector<NodeType>(0);
  };
  BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize, NodeType edgeSize)
      : freeNodes(freeNodesSize), fixedNodes(fixedNodesSize), edges(freeNodesSize) {
    numEdges = edgeSize;

    // Write the (trivial) position of each free node in the edges matrix
    for (NodeType i = 0; i < freeNodesSize; ++i) {
      edges[i].push_back(i);
    }
    // Concrete edges will be added with push_back operations
  };
  NodeType numEdges;
  std::vector<NodeType> freeNodes;
  std::vector<NodeType> fixedNodes;
  std::vector<std::vector<NodeType>> edges;
  // sourceNodeIDs are the indices, targets are the vector's entries. The first entry is the
  // position tho.
  const std::vector<NodeType>& getFreeNodes() const { return freeNodes; }
  const std::vector<NodeType>& getFixedNodes() const { return fixedNodes; }
  const std::vector<std::vector<NodeType>>& getEdges() const { return edges; }

  void insertFreeNode(NodeType node) { freeNodes.push_back(node); };
  void insertFixedNode(NodeType node) { fixedNodes.push_back(node); };

  void addEdge(NodeType sourceID, NodeType targetID) {
    edges[sourceID - fixedNodes.size()].push_back(targetID);  // free Nodes come after the fixed
                                                              // ones
  };

  const NodeType getFreeNodesSize() { return freeNodes.size(); };
  const NodeType getFixedNodesSize() { return fixedNodes.size(); };
  const NodeType getEdgesSize() { return numEdges; };

  void switchNodes(NodeType firstNodeID, NodeType secondNodeID) {
    NodeType& firstPosition = edges[firstNodeID - fixedNodes.size()][0];
    NodeType& secondPosition = edges[secondNodeID - fixedNodes.size()][0];
    std::swap(firstPosition, secondPosition);
  };
};