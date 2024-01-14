//
// Created by alvar on 25/10/2023.
//

#include <iostream>

#include "ds/bipartite_graph.h"

template <typename SizeType>
BipartiteGraph<SizeType> inputGraphManually() {
  BipartiteGraph myGraph = BipartiteGraph<SizeType>();
  int source;
  int target;
  int rightSize;
  int nodeID = 0;
  bool out = true;
  int leftRight = 0;
  while (out) {
    switch (leftRight) {
      case 0:
        std::cout << "Input source of edge: " << nodeID << std::endl;
        std::cin >> source;
        myGraph.insertFreeNode(nodeID);
        while (static_cast<SizeType>(myGraph.edges.size()) < myGraph.freeNodes[nodeID - 1]) {
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
        while (static_cast<SizeType>(myGraph.fixedNodes.size()) < rightSize) {
          myGraph.fixedNodes.push_back(0);
        }
        out = false;
        break;
    }
  }
  return myGraph;
}