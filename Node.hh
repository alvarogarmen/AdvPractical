//
// Created by alvar on 24/10/2023.
//

#ifndef ADVPRACTICAL_NODE_H
#define ADVPRACTICAL_NODE_H

struct Node{
    Node(int nodeID, int outDegree);
    int nodeID;
    int outDegree;
};

Node::Node(int nodeID, int outDegree) {
    this->nodeID=nodeID;
    this->outDegree=outDegree;
}
#endif //ADVPRACTICAL_NODE_H
