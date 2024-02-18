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
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      if (graph.getLeftCrossings(nodeId) == 0) {
        NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
        if (graph.getLeftCrossings(neighbourId) > 0) {
          if (graph.switchNeighbours(nodeId, neighbourId, true)) {
            didChange = true;
          }
        }
      }
    }

    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      if (graph.getRightCrossings(nodeId) == 0) {
        NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
        if (graph.getRightCrossings(neighbourId) > 0) {
          if (graph.switchNeighbours(neighbourId, nodeId, true)) {
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

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of u
 * and the left crossings of v are bigger than right crossings of v,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r2(Graph& graph) {
  using NodeType = typename Graph::NodeType;
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if (graph.getRightCrossings(nodeId) > graph.getLeftCrossings(nodeId) &&
          graph.getLeftCrossings(neighbourId) > graph.getRightCrossings(neighbourId)) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if (graph.getLeftCrossings(nodeId) > graph.getRightCrossings(nodeId) &&
          graph.getRightCrossings(neighbourId) > graph.getLeftCrossings(neighbourId)) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }
    if (didChange) {
      madeSwitch = true;
    }
  }
  return madeSwitch;
}

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of v
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r3(Graph& graph) {
  using NodeType = typename Graph::NodeType;
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if (graph.getRightCrossings(nodeId) > graph.getLeftCrossings(neighbourId)) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if (graph.getLeftCrossings(nodeId) > graph.getRightCrossings(neighbourId)) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }
    if (didChange) {
      madeSwitch = true;
    }
  }
  return madeSwitch;
}