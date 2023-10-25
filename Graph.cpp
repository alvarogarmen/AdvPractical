//
// Created by alvar on 23/10/2023.
//


#include <vector>
#include <unordered_map>
#include <iostream>
#include "Node.h"
#include "Edge.h"
#include "Graph.h"
//TODO: end implementing the edgeHash for the switching



Graph::Graph() {
    this->left=std::vector<int>(0);
    this->right=std::vector<int>(0);
    this->edges=std::vector<Edge>(0,Edge(0,0));
}
Graph::Graph(int leftSize, int rightSize, int edgeSize) {
    this->left=std::vector<int>(leftSize);
    this->right=std::vector<int>(rightSize);
    this->edges=std::vector<Edge>(edgeSize, Edge(0,0));
}
std::vector<int> Graph::getLeftNodes() {
    return this->left;
}
std::vector<int> Graph::getRightNodes() {
    return this->right;
}
std::vector<Edge> Graph::getEdges() {
    return this->edges;
}

void Graph::insertLeftNode(Node node) {
    this->left.push_back((left.size() == 0) ? node.outDegree : left[left.size() - 1] + node.outDegree);
    if (this->edgeHash[node.nodeID].empty()){
        this->edgeHash[node.nodeID].push_back(left.size()-1);
    }
    else{
        this->edgeHash[node.nodeID][0]=left.size()-1;
    }
}

void Graph::insertRightNode(Node node) {
    this->right.push_back(node.outDegree);
}

void Graph::insertEdgeAtIndex(int &sourceID, int &targetID, int index) {
    this->edges[index].source=sourceID;
    this->edges[index].target=targetID;
}

void Graph::insertEdge(int sourceID, int targetID) {
    this->edges.push_back(Edge(sourceID, targetID));
    this->edgeHash[sourceID].push_back(targetID);
}

void Graph::switchNodes(int firstNodeID, int secondNodeID) {
    int& firstPosition = this->edgeHash[firstNodeID][0];
    int& secondPosition = this->edgeHash[secondNodeID][0];

    std::swap(firstPosition, secondPosition);
}
