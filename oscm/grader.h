//
// Created by alvar on 25/10/2023.
//
#ifndef GRADER_H
#define GRADER_H

#include "../ds/edge.h"
#include "../ds/graph.h"

template <typename SizeType>
bool edgeCross(Edge<SizeType>& edge1, Edge<SizeType>& edge2) {
  if ((edge1.source < edge2.source && edge1.target > edge2.target) ||
      edge1.source > edge2.source && edge1.target < edge2.source) {
    return true;
  }
  return false;
}
template <typename SizeType>
int crossGrader(Graph<SizeType>& myGraph) {
  int crossings = 0;
  for (int i = 0; i < myGraph.edges.size(); i++) {
    for (int j = i; j < myGraph.edges.size(); j++) {
      if (edgeCross(myGraph.edges[i], myGraph.edges[j])) {
        crossings++;
      }
    }
  }
  return crossings;
}
//Add the bools instead of conditional addition ->std::accumulate or transform reduce
//Henrik does not like for loops, use algorithm lib for c++

#endif