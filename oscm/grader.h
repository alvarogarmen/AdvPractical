#pragma once
#include <cstddef>
#include <iostream>
template <typename BipartiteGraphType>
int crossGrader(BipartiteGraphType& myGraph) {  // Looks very bad but it is needed to iterate over
                                                // all edge combinations
  int crossings = 0;
  using NodeType = typename BipartiteGraphType::NodeType;

  for (NodeType i = myGraph.getFreeNodes()[0]; i < myGraph.getFreeNodesSize(); i++) {
    for (size_t j = 1; j < (myGraph.edges[i].size());
         j++) {  // First entry in "edges" is the position of the node in the FreeNodes array
      for (NodeType k = i + 1; k < myGraph.getFreeNodesSize(); k++) {
        for (size_t l = 1; l < (myGraph.edges[k].size()); l++) {
          crossings += (myGraph.getEdge(i, 0) < myGraph.getEdge(k, 0) &&
                        myGraph.getEdge(i, j) > myGraph.getEdge(k, l));
        }
      }
    }
  }

  return crossings;
}