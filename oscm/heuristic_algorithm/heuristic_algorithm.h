#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

#include "ds/heuristic_graph/heuristic_graph.h"
/*
 * For each  pair of neighbours  u , v with u < v, if the left crossings of u equal 0
 * and the left crossings of v are bigger than 0,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r1(Graph& graph) {
  using NodeType = typename Graph::NodeType;
  std::vector<NodeType> permutation = graph.getPermutation();
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < permutation.size() - 1; ++i) {
      NodeType nodeId = permutation[i];
      if (graph.getLeftCrossings(nodeId) == 0) {
        NodeType neighbourId = permutation[i + 1];
        if (graph.getLeftCrossings(neighbourId) > 0) {
          if (graph.switchNeighbours(nodeId, neighbourId, true)) {
            didChange = true;
          }
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = permutation.size() - 1; i > 0; --i) {
      NodeType nodeId = permutation[i];
      if (graph.getRightCrossings(nodeId) == 0) {
        NodeType neighbourId = permutation[i - 1];
        if (graph.getRightCrossings(neighbourId) > 0) {
          if (graph.switchNeighbours(nodeId, neighbourId, true)) {
            didChange = true;
          }
        }
      }
    }
    if (didChange) {
      madeSwitch = true;
    }
  }
  return madeSwitch;
}