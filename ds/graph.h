//
// Created by alvar on 25/10/2023.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>

#include "edge.h"
#include "node.h"

template<typename SizeType>
struct Graph{
    Graph();
    Graph(SizeType leftSize, SizeType rightSize, SizeType edgeSize);
    std::vector<SizeType> left;
    std::vector<SizeType> right;
    std::vector<Edge<SizeType>> edges;         //Not sure if actually needed
    std::unordered_map<SizeType, std::vector<SizeType>> edgeHash ;  //nodeIDs are used as handles. The first entry in the array is the
    //position of the nodeID in the left graph
    std::vector<SizeType> getLeftNodes();
    std::vector<SizeType> getRightNodes();
    std::vector<Edge<SizeType>> getEdges();

    void insertLeftNode(Node<SizeType> node);
    void insertRightNode(Node<SizeType> node);
    void insertEdgeAtIndex(SizeType sourceID, SizeType targetID, SizeType index);
    void insertEdge(SizeType sourceID, SizeType targetID);        //nodeID represents the target of the edge. The index in combination with the left
    //vector gives the source of the edge.

    void switchNodes(SizeType firstNodeID, SizeType secondNodeID);
};

template<typename SizeType>
Graph<SizeType>::Graph() {
    this->left = std::vector<SizeType>(0);
    this->right = std::vector<SizeType>(0);
}

template<typename SizeType>
Graph<SizeType>::Graph(SizeType leftSize, SizeType rightSize, SizeType edgeSize) {
    this->left = std::vector<SizeType>(leftSize);
    this->right = std::vector<SizeType>(rightSize);
    this->edges = std::vector<Edge<SizeType>>(edgeSize, Edge(0, 0));
}

template<typename SizeType>
std::vector<SizeType> Graph<SizeType>::getLeftNodes() {
    return this->left;
}

template<typename SizeType>
std::vector<SizeType> Graph<SizeType>::getRightNodes() {
    return this->right;
}

template<typename SizeType>
std::vector<Edge<SizeType>> Graph<SizeType>::getEdges() {
    return this->edges;
}

template<typename SizeType>
void Graph<SizeType>::insertLeftNode(Node<SizeType> node) {
    this->left.push_back((this->left.empty()) ? node.outDegree : this->left.back() + node.outDegree);
    if (this->edgeHash[node.nodeID].empty()) {
        this->edgeHash[node.nodeID].push_back(this->left.size() - 1);
    }
    else {
        this->edgeHash[node.nodeID][0] = this->left.size() - 1;
    }
}

template<typename SizeType>
void Graph<SizeType>::insertRightNode(Node<SizeType> node) {
    this->right.push_back(node.outDegree);
}

template<typename SizeType>
void Graph<SizeType>::insertEdgeAtIndex(SizeType sourceID, SizeType targetID, SizeType index) {
    this->edges[index].source = sourceID;
    this->edges[index].target = targetID;
}

template<typename SizeType>
void Graph<SizeType>::insertEdge(SizeType sourceID, SizeType targetID) {
    this->edges.push_back(Edge(sourceID, targetID));
    this->edgeHash[sourceID].push_back(targetID);
}

template<typename SizeType>
void Graph<SizeType>::switchNodes(SizeType firstNodeID, SizeType secondNodeID) {
    SizeType &firstPosition = this->edgeHash[firstNodeID][0];
    SizeType &secondPosition = this->edgeHash[secondNodeID][0];
    std::swap(firstPosition, secondPosition);
}


#endif