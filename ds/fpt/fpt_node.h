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
class YNode : public FptNode<NodeType> {
 public:
  YNode(NodeType nodeID, std::vector<NodeType> neighbours, NodeType position)
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
NodeType YNode<NodeType>::getDxScannedIndex() {
  return dxScannedIndex;
}

template <typename NodeType>
void YNode<NodeType>::setDxScannedIndex(NodeType newVal) {
  dxScannedIndex = newVal;
}

template <typename NodeType>
class XNode : public FptNode<NodeType> {
 public:
  XNode(NodeType nodeID, std::vector<YNode<NodeType>> neighbours)
      : FptNode<NodeType>(nodeID), neighbours(neighbours) {}
  std::vector<YNode<NodeType>> neighbours;
  // saves the node ids of y nodes thet x is positioned between there smallest and largest
  // neighbours
  std::vector<NodeType> yx;
};