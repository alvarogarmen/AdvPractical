#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/UndoAlgorithmStep.h"
#include "ds/reduction_graph/reduction_graph.h"
#include "oscm/reduction_algorithm/compute_crossings.h"

/*
 * For each decision u < v that we make adjust the leftRightSet of u and v
 * Add all the transitiv decisions that follow the u < v decision
 */
template <class Graph, class Undo>
void parameterAccounting(Graph& graph, typename Graph::NodeType u, typename Graph::NodeType v,
                         typename Graph::CrossingCountType* currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  if (u != v) {
    if (graph.getRightNodes(u).find(v) == graph.getRightNodes(u).end()) {
      graph.insertRightNode(u, v);
      graph.insertLeftNode(v, u);
      *currentSolution += graph.getCrossing(u, v);
      if (undo) {
        undo->addParameterAccountingUndo(u, v, graph.getCrossing(u, v), graph.getCrossing(v, u));
      }
      graph.deleteCrossings(u, v);
      for (NodeType smallerThanU : graph.getLeftNodes(u)) {
        parameterAccounting<Graph, Undo>(graph, smallerThanU, v, currentSolution, undo);
        for (NodeType biggerThanV : graph.getRightNodes(v)) {
          parameterAccounting<Graph, Undo>(graph, smallerThanU, biggerThanV, currentSolution, undo);
          parameterAccounting<Graph, Undo>(graph, u, biggerThanV, currentSolution, undo);
        }
      }
    }
  }
}

//  if a free node v is comparable  with all other free nodes, then put v in its right fixed
//  position
template <class Graph, class Undo>
void rrlo1(Graph& graph, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  for (NodeType v = 0; v < graph.getFreeNodesSize(); ++v) {
    if (graph.getNodeCrossing(v).size() == 0) {
      NodeType endIndex = graph.getLeftNodes(v).size();
      graph.setFixedPosition(v, endIndex);
      if (undo) {
        undo->addSetPositionUndo(endIndex);
      }
    }
  }
}

// If {u, v} is an incomparable pair in which, u and v are comparable with all other nodes, with
// c(u, v) <= c(v, u) , then commit u < v, and do the parameter accounting.
template <class Graph, class Undo>
bool rrlo2(Graph& graph, typename Graph::CrossingCountType* currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;

  NodeType n = graph.getFreeNodesSize();
  bool didChange = false;
  for (NodeType u = 0; u < n; ++u) {
    // only one node is missing for u
    if (graph.getLeftNodes(u).size() + graph.getRightNodes(u).size() == n - 2) {
      for (NodeType v = u + 1; v < n; ++v) {
        // only one node missing for v
        // same number of left and right nodes for v and u ->
        // v and u are neighbours and v is missing for u and the other way around
        if (graph.getLeftNodes(v).size() + graph.getRightNodes(v).size() == n - 2 &&
            graph.getLeftNodes(v).size() == graph.getLeftNodes(u).size()) {
          if (graph.getCrossing(u, v) <= graph.getCrossing(v, u)) {
            parameterAccounting<Graph, Undo>(graph, u, v, currentSolution, undo);
            didChange = true;
          } else {
            parameterAccounting<Graph, Undo>(graph, v, u, currentSolution, undo);
            didChange = true;
          }
        }
      }
    }
  }
  return didChange;
}

// if we have c(u, v) = 1 and c(v, u) = 2 with d(u) == 2  d(v) == 2 then commit u < v and do
// parameter accounting
template <class Graph, class Undo>
void rr3(Graph& graph, typename Graph::CrossingCountType* currentSolution, Undo* undo = nullptr) {
  // Create a list of pairs to be modified
  using NodeType = typename Graph::NodeType;

  std::vector<std::pair<NodeType, NodeType>> pairsToModify;

  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    for (auto [v, crossingValue] : graph.getNodeCrossing(u)) {
      if (crossingValue == 1 && graph.getNodeCrossing(v).at(u) == 2 &&
          graph.getFreeNodeNeighboursSize(u) == 2 && graph.getFreeNodeNeighboursSize(v) == 2) {
        pairsToModify.emplace_back(u, v);
      }
    }
  }

  // Modify the pairs outside the loop
  for (auto [u, v] : pairsToModify) {
    parameterAccounting<Graph, Undo>(graph, u, v, currentSolution, undo);
  }
}

// For each pair of free nodes u, v with N(u) = N(v),(arbitrarily) commit a < b, and do parameter
// accounting.
template <class Graph, class Undo>
void rr2(Graph& graph, typename Graph::CrossingCountType* currentSolution) {
  using NodeType = typename Graph::NodeType;
  NodeType n = graph.getFreeNodesSize();
  for (NodeType u = 0; u < n; ++u) {
    for (NodeType v = u + 1; v < n; ++v) {
      if (std::equal(graph.getFreeNodeNeighbours(u).begin(), graph.getFreeNodeNeighbours(u).end(),
                     graph.getFreeNodeNeighbours(v).begin(),
                     graph.getFreeNodeNeighbours(v).end())) {
        parameterAccounting<Graph, Undo>(graph, u, v, currentSolution);
      }
    }
  }
}

// If c(u, v) > k, then commit v < u and do the parameter accounting.
template <class Graph, class Undo>
bool rrLarge(Graph& graph, typename Graph::CrossingCountType crossingsLeft,
             typename Graph::CrossingCountType* currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  std::vector<std::pair<NodeType, NodeType>> pairsToModify;
  bool didChange = false;
  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    for (auto [v, crossingValue] : graph.getNodeCrossing(u)) {
      if (crossingValue > crossingsLeft) {
        pairsToModify.emplace_back(v, u);
        didChange = true;
      }
    }
  }
  for (auto [v, u] : pairsToModify) {
    parameterAccounting<Graph, Undo>(graph, v, u, currentSolution, undo);
  }
  return didChange;
}

// Check if there is  an incomparable i/j pattern {u, v} with i + j >= 4
template <class Graph>
bool IJBiggerThenFour(const Graph& graph, typename Graph::NodeType* u,
                      typename Graph::NodeType* v) {
  using NodeType = typename Graph::NodeType;
  using CrossingCountType = typename Graph::CrossingCountType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      CrossingCountType secondFirstcrossingValue = graph.getCrossing(secondNode, firstNode);
      if (FirstSecondcrossingValue + secondFirstcrossingValue >= 4) {
        *u = firstNode;
        *v = secondNode;
        return true;
      }
    }
  }
  return false;
}

// Check if there is an there is a dependent 2/1 pattern {u, v} with c(u, v) + c(v, u) = 3
template <class Graph>
bool IJEqualToThree(const Graph& graph, typename Graph::NodeType* u, typename Graph::NodeType* v) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      if (FirstSecondcrossingValue == 2) {
        *u = firstNode;
        *v = secondNode;
        return true;
      }
    }
  }
  return false;
}

// Check if there is a 1/1 pattern {u, v}
template <class Graph>
bool IJEqualToTwo(const Graph& graph, typename Graph::NodeType* u, typename Graph::NodeType* v) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      *u = firstNode;
      *v = secondNode;
      return true;
    }
  }
  return false;
}
template <class Graph, class Undo>
void algorithmStep(Graph& graph, typename Graph::CrossingCountType currentSolution, bool isInitStep,
                   typename Graph::NodeType leftNode, typename Graph::NodeType rightNode) {
  using NodeType = typename Graph::NodeType;
  Undo undo;
  if (isInitStep) {
    parameterAccounting<Graph, Undo>(graph, leftNode, rightNode, &currentSolution, &undo);
  }
  bool didChangeRrlo2 = true;
  bool didChangeRrLarge = true;
  while (didChangeRrlo2 || didChangeRrLarge) {
    didChangeRrlo2 = rrlo2(graph, &currentSolution, &undo);
    didChangeRrLarge =
        rrLarge(graph, graph.getBestSolution() - currentSolution, &currentSolution, &undo);
    rrlo1(graph, &undo);
  }
  if (graph.getBestSolution() <= currentSolution) {
    graph.doUndo(undo);
    return;
  }
  NodeType u;
  NodeType v;
  if (IJBiggerThenFour(graph, &u, &v)) {
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u);
  } else if (IJEqualToThree(graph, &u, &v)) {
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u);
  } else if (IJEqualToTwo(graph, &u, &v)) {
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v);
  }
  graph.setBestSolution(currentSolution, graph.getFixedPosition());
  graph.doUndo(undo);
  return;
}

template <class Graph, class Undo>
void algorithm(Graph& graph) {
  using CrossingCountType = typename Graph::CrossingCountType;
  computeCrossingSums<Graph, Undo>(graph);
  CrossingCountType currentSolution = 0;
  rr2<Graph, Undo>(graph, &currentSolution);
  bool didChangeRrlo2 = true;
  bool didChangeRrLarge = true;

  while (didChangeRrlo2 || didChangeRrLarge) {
    didChangeRrlo2 = rrlo2<Graph, Undo>(graph, &currentSolution);
    didChangeRrLarge =
        rrLarge<Graph, Undo>(graph, graph.getBestSolution() - currentSolution, &currentSolution);
    rrlo1<Graph, Undo>(graph);
  }

  rr3<Graph, Undo>(graph, &currentSolution);
  algorithmStep<Graph, Undo>(graph, currentSolution, false, 0, 0);
}