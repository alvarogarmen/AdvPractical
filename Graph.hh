//
// Created by alvar on 23/10/2023.
//

#ifndef ADVPRACTICAL_GRAPH_H
#define ADVPRACTICAL_GRAPH_H
#include <vector>
#include "Node.hh"
#include "Edge.hh"

struct Graph{
    Graph();
    Graph(int leftSize, int rightSize, int edgeSize);
    std::vector<int> left;
    std::vector<int> right;
    std::vector<Edge> edges;         //Not sure if actually needed
    std::vector<int> getLeftNodes();
    std::vector<int> getRightNodes();
    std::vector<Edge> getEdges();

    void insertLeftNode(Node node);
    void insertRightNode(Node node);
    void insertEdgeAtIndex(int& sourceID, int& targetID, int index);
    void insertEdge(int sourceID, int targetID);        //nodeID represents the target of the edge. The index in combination with the left
};                                                        //vector gives the source of the edge.

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
}
#endif //ADVPRACTICAL_GRAPH_H
