//
// Created by alvar on 25/10/2023.
//

#pragma once

#include <unordered_map>
#include <vector>

template <typename NT>
class BipartiteGraph {
 public:
  using NodeType = NT;

  BipartiteGraph() : freeNodes({}), fixedNodes({}){};

  BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize, NodeType edgeSize)
      : freeNodes(freeNodesSize), fixedNodes(fixedNodesSize), edges(freeNodesSize) {
    numEdges = edgeSize;
    // Write trivial positions in both vectors
    for (NodeType i = 0; i < freeNodesSize; i++) {
      freeNodes[i] = i;
    }
    for (NodeType i = 0; i < fixedNodesSize; i++) {
      fixedNodes[i] = i;
    }
    // Concrete edges will be added with push_back operations
  };

  const std::vector<NodeType>& getFreeNodes() const { return freeNodes; }
  const std::vector<NodeType>& getFixedNodes() const { return fixedNodes; }
  const std::vector<std::vector<NodeType>>& getEdges() const { return edges; }

  void insertFreeNode(NodeType node) { freeNodes.push_back(node); };
  void insertFixedNode(NodeType node) { fixedNodes.push_back(node); };

  void addEdge(NodeType sourceID, NodeType targetID) {
    edges[sourceID].push_back(targetID);  // free Nodes come after the fixed
                                          // ones
  };
  NodeType getEdge(NT nodeID, NT number) { return auto const & edges[nodeIT] }
  const NodeType getFreeNodesSize() { return freeNodes.size(); };
  const NodeType getFixedNodesSize() { return fixedNodes.size(); };
  const NodeType getEdgesSize() { return numEdges; };

  void switchNodes(NodeType firstNodeID, NodeType secondNodeID) {
    NodeType& firstPosition = freeNodes[firstNodeID];
    NodeType& secondPosition = freeNodes[secondNodeID];
    std::swap(firstPosition, secondPosition);
  };

 private:
  NodeType numEdges;
  std::vector<NodeType> freeNodes;  // Stores the position of each node, e.g. freeNodes[0]=position
  std::vector<NodeType> fixedNodes;
  std::vector<std::vector<NodeType>> edges;
  // sourceNodeIDs are the indices (FreeNodes), targets are the vector's entries (FixedNodes). The
  // first entry is the position tho.
};