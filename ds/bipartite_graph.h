//
// Created by alvar on 25/10/2023.
//

#pragma once

#include <unordered_map>
#include <vector>

template <typename NodeType>
struct Node {
  Node(NodeType nodeID, NodeType outDegree);
  int nodeID{};
  int outDegree{};
};
template <typename NodeType>
Node<NodeType>::Node(NodeType nodeID, NodeType outDegree) {
  this->nodeID = nodeID;
  this->outDegree = outDegree;
}

template <typename NodeType>
struct Edge {
  Edge(NodeType source, NodeType target);
  NodeType source;
  NodeType target;
};

template <typename NodeType>
Edge<NodeType>::Edge(NodeType source, NodeType target) {
  this->source = source;
  this->target = target;
}

template <typename NT>
struct BipartiteGraph {
  using NodeType = NT;
  BipartiteGraph() {
    this->freeNodes = std::vector<NodeType>(0);
    this->fixedNodes = std::vector<NodeType>(0);
  };
  BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize, NodeType edgeSize)
      : freeNodes(freeNodesSize),
        fixedNodes(fixedNodesSize){
            // Edges will be added with push_back operations
        };

  std::vector<NodeType> freeNodes;
  std::vector<NodeType> fixedNodes;
  std::vector<Edge<NodeType>> edges;  // Not sure if actually needed
  std::vector<NodeType, std::vector<NodeType>>
      edgeHash;  // nodeIDs are the indices. The first entry in the array is the
  // position of the nodeID in the freeNodes graph
  const std::vector<NodeType>& getFreeNodes() { return this->freeNodes&; };
  const std::vector<NodeType>& getFixedNodes() { return this->fixedNodes&; };
  const std::vector<Edge<NodeType>>& getEdges() { return this->edges&; };

  void insertFreeNode(Node<NodeType> node) {
    this->freeNodes.push_back((this->freeNodes.empty()) ? node.outDegree
                                                        : this->freeNodes.back() + node.outDegree);
    if (this->edgeHash[node.nodeID].empty()) {
      this->edgeHash[node.nodeID].push_back(this->freeNodes.size() - 1);
    } else {
      this->edgeHash[node.nodeID][0] = this->freeNodes.size() - 1;
    }
  };
  void insertFixedNode(Node<NodeType> node) { this->fixedNodes.push_back(node.outDegree); };
  void insertEdgeAtIndex(NodeType sourceID, NodeType targetID, NodeType index) {
    this->edges[index].source = sourceID;
    this->edges[index].target = targetID;
  };
  void addEdge(NodeType sourceID, NodeType targetID) {
    this->edges.push_back(Edge(sourceID, targetID));
    this->edgeHash[sourceID].push_back(targetID);
  };  // nodeID represents the target of the edge. The index in
      // combination with the left
  // vector gives the source of the edge.
  const NodeType getFreeNodesSize() { return this->freeNodes.size(); };
  const NodeType getFixedNodesSize() { return this->fixedNodes.size(); };
  const NodeType getEdgesSize() { return this->edges.size(); };
  void switchNodes(NodeType firstNodeID, NodeType secondNodeID) {
    NodeType& firstPosition = this->edgeHash[firstNodeID][0];
    NodeType& secondPosition = this->edgeHash[secondNodeID][0];
    std::swap(firstPosition, secondPosition);
  };
};