#pragma once

#include <stdio.h>

#include <algorithm>
#include <atomic>
#include <csignal>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#include <vector>

#include "ds/heuristic_graph/heuristic_graph.h"
#include "oscm/median_algorithm/median_algorithm.h"
namespace heuristic_algorithm {
/*
 * For each  pair of neighbours  u , v with u < v, if the left crossings of u equal 0
 * and the left crossings of v are bigger than 0,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r1(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getLeftCrossings(nodeId) == 0 && graph.getLeftCrossings(neighbourId) > 0;
  } else {
    return graph.getRightCrossings(nodeId) == 0 && graph.getRightCrossings(neighbourId) > 0;
  }
}

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of u
 * and the left crossings of v are bigger than right crossings of v,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r2(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getRightCrossings(nodeId) > graph.getLeftCrossings(nodeId) &&
           graph.getLeftCrossings(neighbourId) > graph.getRightCrossings(neighbourId);
  } else {
    return graph.getLeftCrossings(nodeId) > graph.getRightCrossings(nodeId) &&
           graph.getRightCrossings(neighbourId) > graph.getLeftCrossings(neighbourId);
  }
}

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of v
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r3(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getRightCrossings(nodeId) > graph.getLeftCrossings(neighbourId);
  } else {
    return graph.getLeftCrossings(nodeId) > graph.getRightCrossings(neighbourId);
  }
}

template <class Graph>
bool heuristicAlgorithm(Graph& graph, bool runR1, bool runR2, bool runR3,
                        std::atomic<bool>& terminationRequested) {
  using NodeType = typename Graph::NodeType;

  median_algorithm::medianAlgorithm(graph);

  for (int i = 0; i < graph.getFreeNodesSize(); ++i) {
    auto permutation = graph.getPermutation();
    graph.setFreeNodePosition(permutation[i], i);
  }
  graph.computeCrossingSums();
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if ((runR1 && r1(graph, nodeId, neighbourId)) || (runR2 && r2(graph, nodeId, neighbourId)) ||
          (runR3 && r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if ((runR1 && r1(graph, nodeId, neighbourId)) || (runR2 && r2(graph, nodeId, neighbourId)) ||
          (runR3 && r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }

    if (didChange) {
      madeSwitch = true;
    }
  }

  int bestSolution = graph.getCrossings();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, graph.getFreeNodesSize() - 1);

  while (!terminationRequested.load()) {
    NodeType node1 = dis(gen);
    NodeType node2 = dis(gen);

    if (node1 != node2) {
      graph.switchNodes(node1, node2, bestSolution);
    }
  }

  return madeSwitch;
}
}  // namespace heuristic_algorithm