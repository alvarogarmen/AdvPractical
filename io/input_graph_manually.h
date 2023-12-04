//
// Created by alvar on 25/10/2023.
//

#include <iostream>

#include "ds/bipartite_graph.h"

template <typename BipartiteGraph>
BipartiteGraph inputGraphManually() {
  BipartiteGraph myGraph = BipartiteGraph();
  int source;
  int target;
  int degree;
  int rightSize;
  int nodeID = 1;
  bool out = true;
  int leftRight = 0;
  while (out) {
    switch (leftRight) {
      case 0:
        std::cout << "Input degree of node: " << nodeID << std::endl;
        std::cin >> degree;
        myGraph.insertFreeNode(Node(nodeID, degree));
        while (myGraph.edges.size() < myGraph.freeNodes[nodeID - 1]) {
          std::cout << "Input target of edge" << std::endl;
          std::cin >> target;
          myGraph.addEdge(nodeID, target);
        }
        nodeID++;
        std::cout << "Input 0 if new node on left side, 1 if go to right side" << std::endl;
        std::cin >> leftRight;
        break;
      case 1:
        nodeID = 1;
        std::cout << "Input number of nodes on right side" << std::endl;
        std::cin >> rightSize;
        while (myGraph.fixedNodes.size() < rightSize) {
          myGraph.fixedNodes.push_back(0);
        }
        out = false;
        break;
    }
  }
  return myGraph;
}