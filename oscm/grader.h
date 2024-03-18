#pragma once
#include <cstddef>
#include <iostream>
template <typename BipartiteGraphType>
int crossGrader(BipartiteGraphType& myGraph) {  // Compare from left to right if edges cross
  int crossings = 0;
  using NodeType = typename BipartiteGraphType::NodeType;

  for (NodeType i = myGraph.getFreeNodes()[0]; i < myGraph.getFreeNodesSize(); i++) {

    for (size_t j = 0; j < (myGraph.getOutEdges(i).size()); j++) {
      for (NodeType k = i; k < myGraph.getFreeNodesSize(); k++) {
        for (size_t l = 0; l < (myGraph.getOutEdges(k).size()); l++) {
          crossings += (myGraph.getFreeNodes()[i] < myGraph.getFreeNodes()[k] &&
                        myGraph.getEdge(i, j) > myGraph.getEdge(k, l));
        }
      }
    }
  }
  return crossings;
}