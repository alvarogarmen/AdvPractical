//
// Created by alvar on 25/10/2023.
//
#ifndef NODE_H
#define NODE_H

template <typename NodeType>
struct Node{
    Node(NodeType nodeID, NodeType outDegree);
    int nodeID{};
    int outDegree{};


};
template <typename NodeType>
Node<NodeType>::Node(NodeType nodeID, NodeType outDegree) {
    this->nodeID=nodeID;
    this->outDegree=outDegree;
}

#endif