//
// Created by alvar on 25/10/2023.
//
#pragma once

#include "ds/bipartite_graph.h"

template <typename SizeType>
bool edgeCross(SizeType& sourceEdge1, SizeType& targetEdge1, SizeType& sourceEdge2,
               SizeType& targetEdge2) {
  if (((sourceEdge1 < sourceEdge2 && targetEdge1 > targetEdge2)) ||
      (sourceEdge1 > sourceEdge2 && targetEdge1 < sourceEdge2)) {
    return true;
  }
  return false;
}
template <typename BipartiteGraphType, typename SizeType>
int crossGrader(BipartiteGraphType& myGraph) {  // Looks very bad but it is needed to iterate over
                                                // all edge combinations
  int crossings = 0;
  for (SizeType i = 0; i < myGraph.getFreeNodesSize(); i++) {                        // O(freeNodes)
    for (SizeType j = 0; j < static_cast<SizeType>(myGraph.edges[i].size()); j++) {  // O(maxDegree)
      for (SizeType k = 1; k < myGraph.getFreeNodesSize(); k++) {                    // O(freeNodes)
        for (SizeType l = 0;                                                         // O(maxDegree)
             l < static_cast<SizeType>(myGraph.edges[k].size()); l++) {
          if (edgeCross(i, myGraph.edges[i][j], k, myGraph.edges[k][l])) {
            crossings++;
          }
        }
      }
    }
  }

  return crossings;
}
// Add the bools instead of conditional addition ->std::accumulate or transform reduce
// Henrik does not like for loops, use algorithm lib for c++