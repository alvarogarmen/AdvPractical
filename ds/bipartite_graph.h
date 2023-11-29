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
  BipartiteGraph();
  BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize, NodeType edgeSize);
  std::vector<NodeType> freeNodes;
  std::vector<NodeType> fixedNodes;
  std::vector<Edge<NodeType>> edges;  // Not sure if actually needed
  std::unordered_map<NodeType, std::vector<NodeType>>
      edgeHash;  // nodeIDs are used as handles. The first entry in the array is the
  // position of the nodeID in the freeNodes graph
  std::vector<NodeType> getFreeNodes();
  std::vector<NodeType> getFixedNodes();
  std::vector<Edge<NodeType>> getEdges();

  void insertFreeNode(Node<NodeType> node);
  void insertFixedNode(Node<NodeType> node);
  void insertEdgeAtIndex(NodeType sourceID, NodeType targetID, NodeType index);
  void addEdge(NodeType sourceID,
               NodeType targetID);  // nodeID represents the target of the edge. The index in
                                    // combination with the left
  // vector gives the source of the edge.
  NodeType getFreeNodesSize();
  NodeType getFixedNodesSize();
  NodeType getEdgesSize();
  void switchNodes(NodeType firstNodeID, NodeType secondNodeID);
};

template <typename NodeType>
BipartiteGraph<NodeType>::BipartiteGraph() {
  this->freeNodes = std::vector<NodeType>(0);
  this->fixedNodes = std::vector<NodeType>(0);
}

template <typename NodeType>
BipartiteGraph<NodeType>::BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize,
                                         NodeType edgeSize) {
  this->freeNodes = std::vector<NodeType>(freeNodesSize);
  this->fixedNodes = std::vector<NodeType>(fixedNodesSize);
  // Edges will be added with push_back operations
}

template <typename NodeType>
std::vector<NodeType> BipartiteGraph<NodeType>::getFreeNodes() {
  return this->freeNodes;
}

template <typename NodeType>
std::vector<NodeType> BipartiteGraph<NodeType>::getFixedNodes() {
  return this->fixedNodes;
}

template <typename NodeType>
std::vector<Edge<NodeType>> BipartiteGraph<NodeType>::getEdges() {
  return this->edges;
}

template <typename NodeType>
void BipartiteGraph<NodeType>::insertFreeNode(Node<NodeType> node) {
  this->freeNodes.push_back((this->freeNodes.empty()) ? node.outDegree
                                                      : this->freeNodes.back() + node.outDegree);
  if (this->edgeHash[node.nodeID].empty()) {
    this->edgeHash[node.nodeID].push_back(this->freeNodes.size() - 1);
  } else {
    this->edgeHash[node.nodeID][0] = this->freeNodes.size() - 1;
  }
}

template <typename NodeType>
void BipartiteGraph<NodeType>::insertFixedNode(Node<NodeType> node) {
  this->fixedNodes.push_back(node.outDegree);
}

template <typename NodeType>
void BipartiteGraph<NodeType>::insertEdgeAtIndex(NodeType sourceID, NodeType targetID,
                                                 NodeType index) {
  this->edges[index].source = sourceID;
  this->edges[index].target = targetID;
}

template <typename NodeType>
void BipartiteGraph<NodeType>::addEdge(NodeType sourceID, NodeType targetID) {
  this->edges.push_back(Edge(sourceID, targetID));
  this->edgeHash[sourceID].push_back(targetID);
}

template <typename NodeType>
NodeType BipartiteGraph<NodeType>::getFixedNodesSize() {
  return this->fixedNodes.size();
}
template <typename NodeType>
NodeType BipartiteGraph<NodeType>::getFreeNodesSize() {
  return this->freeNodes.size();
}

template <typename NodeType>
NodeType BipartiteGraph<NodeType>::getEdgesSize() {
  return this->edges.size();
}

template <typename NodeType>
void BipartiteGraph<NodeType>::switchNodes(NodeType firstNodeID, NodeType secondNodeID) {
  NodeType& firstPosition = this->edgeHash[firstNodeID][0];
  NodeType& secondPosition = this->edgeHash[secondNodeID][0];
  std::swap(firstPosition, secondPosition);
}