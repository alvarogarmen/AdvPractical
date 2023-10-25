//
// Created by alvar on 25/10/2023.
//
#ifndef NODE_H
#define NODE_H
struct Node{
    Node(int nodeID, int outDegree);
    int nodeID;
    int outDegree;
};

#endif