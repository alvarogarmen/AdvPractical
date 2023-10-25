//
// Created by alvar on 25/10/2023.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
#include "Node.h"
#include "Edge.h"

struct Graph{
    Graph();
    Graph(int leftSize, int rightSize, int edgeSize);
    std::vector<int> left;
    std::vector<int> right;
    std::vector<Edge> edges;         //Not sure if actually needed
    std::unordered_map<int, std::vector<int>> edgeHash ;  //nodeIDs are used as handles. The first entry in the array is the
    //position of the nodeID in the left graph
    std::vector<int> getLeftNodes();
    std::vector<int> getRightNodes();
    std::vector<Edge> getEdges();

    void insertLeftNode(Node node);
    void insertRightNode(Node node);
    void insertEdgeAtIndex(int& sourceID, int& targetID, int index);
    void insertEdge(int sourceID, int targetID);        //nodeID represents the target of the edge. The index in combination with the left
    //vector gives the source of the edge.

    void switchNodes(int firstNodeID, int secondNodeID);
};

#endif