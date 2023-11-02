#pragma once

#include <vector>

template <typename NodeType>
class FptNode {
 public:
  FptNode(NodeType nodeID);
  NodeType getNodeID();

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
class FreeNode : public FptNode<NodeType> {
 public:
  FreeNode(NodeType nodeID, std::vector<NodeType>& neighbours, NodeType position)
      : FptNode<NodeType>(nodeID), neighbours(neighbours), position(position) {
    this->dxScannedIndex = 0;
  }
  std::vector<NodeType> neighbours;
  NodeType getDxScannedIndex();
  void setDxScannedIndex(NodeType newVal);

  NodeType position;
  // dxScannedIndex is the index of the last scanned x
  // used to count amount of the node neighbours smaller than x
  NodeType dxScannedIndex;
};

template <typename NodeType>
NodeType FreeNode<NodeType>::getDxScannedIndex() {
  return dxScannedIndex;
}

template <typename NodeType>
void FreeNode<NodeType>::setDxScannedIndex(NodeType newVal) {
  dxScannedIndex = newVal;
}

template <typename NodeType>
class FixedNode : public FptNode<NodeType> {
 public:
  FixedNode(NodeType nodeID, std::vector<FreeNode<NodeType>>& neighbours) {
    this->nodeID = nodeID;
    this->neighbours = neighbours;
  }
  std::vector<FreeNode<NodeType>> neighbours;  // why put the whole FreeNode here? A FixedNode would
                                               // have its neighbours and its neighbours neighbours
  // saves the node ids of y nodes thet x is positioned between there smallest and largest
  // neighbours
  std::vector<NodeType> yx;
};