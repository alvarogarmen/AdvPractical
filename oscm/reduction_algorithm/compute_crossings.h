#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/reduction_graph.h"
#include "ds/reduction_graph/undo_algorithm_step.h"

namespace reductionalgorithms {
// compute the number of crossings created by the edges from two free nodes (u, v)
// when u is to the left of v
template <class Graph>
typename Graph::CrossingCountType const computeUVcrossing(Graph& graph, typename Graph::NodeType u,
                                                          typename Graph::NodeType v) {
  using CrossingCountType = typename Graph::CrossingCountType;
  CrossingCountType crossingSum = 0;
  for (const auto uNeighbour : graph.getFreeNodeNeighbours(u)) {
    for (const auto vNeighbour : graph.getFreeNodeNeighbours(v)) {
      crossingSum += vNeighbour < uNeighbour;
    }
  }
  return crossingSum;
}

// copmute for each free node u the number of crossings created by the edges from nodes to the
// left of u  and to its right
template <class Graph, class Undo>
void computeCrossingSums(Graph& graph) {
  using NodeType = typename Graph::NodeType;
  using CrossingCountType = typename Graph::CrossingCountType;
  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    for (NodeType v = u + 1; v < graph.getFreeNodesSize(); ++v) {
      CrossingCountType crossingUV = computeUVcrossing(graph, u, v);
      graph.addCrossing(u, v, crossingUV);
      // reduction RR1: For each pair of vertices {u, v} ⊆ free nodes that forms a 0/j pattern
      // with j > 0, commit u < v
      CrossingCountType crossingVU = computeUVcrossing(graph, v, u);
      graph.addCrossing(v, u, crossingVU);
    }
  }
}
}  // namespace reductionalgorithms