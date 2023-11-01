#pragma once

#include <vector>

template <typename NodeType>
class FptNode {
 public:
  FptNode(NodeType nodeID);
  NodeType getNodeID();

 private:
  NodeType nodeID;
};
template <typename NodeType>
FptNode<NodeType>::FptNode(NodeType nodeID) {
  this->nodeID = nodeID;
}

template <typename NodeType>
NodeType FptNode<NodeType>::getNodeID() {
  return nodeID;
}

template <typename NodeType>
class NonFixedNode : public FptNode<NodeType> {
 public:
  NonFixedNode(NodeType nodeID, std::vector<NodeType> neighbours, NodeType position)
      : FptNode<NodeType>(nodeID), neighbours(neighbours), position(position) {
    dxScannedIndex = 0;
  }
  std::vector<NodeType> neighbours;
  NodeType getDxScannedIndex();
  void setDxScannedIndex(NodeType newVal);

 private:
  NodeType position;
  // dxScannedIndex is the index of the last scanned x
  // used to count amount of the node neighbours smaller than x
  NodeType dxScannedIndex;
};

template <typename NodeType>
NodeType NonFixedNode<NodeType>::getDxScannedIndex() {
  return dxScannedIndex;
}

template <typename NodeType>
void NonFixedNode<NodeType>::setDxScannedIndex(NodeType newVal) {
  dxScannedIndex = newVal;
}

template <typename NodeType>
class FixedNode : public FptNode<NodeType> {
 public:
  FixedNode(NodeType nodeID, std::vector<NonFixedNode<NodeType>> neighbours)
      : FptNode<NodeType>(nodeID), neighbours(neighbours) {}
  std::vector<NonFixedNode<NodeType>> neighbours;
  // saves the node ids of y nodes thet x is positioned between there smallest and largest
  // neighbours
  std::vector<NodeType> yx;
};