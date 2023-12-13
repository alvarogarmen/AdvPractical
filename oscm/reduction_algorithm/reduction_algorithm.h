#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <tuple>
#include <vector>

#include "ds/reduction_graph/reduction_graph.h"
#include "ds/reduction_graph/undo_algorithm_step.h"
#include "oscm/reduction_algorithm/compute_crossings.h"

/*
 * For each decision u < v that we make adjust the leftRightSet of u and v
 * Add all the transitiv decisions that follow the u < v decision
 */
template <class Graph, class Undo>
void parameterAccounting(Graph& graph, typename Graph::NodeType u, typename Graph::NodeType v,
                         typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  if (u != v) {
    if (graph.getRightNodes(u).find(v) == graph.getRightNodes(u).end()) {
      graph.insertRightNode(u, v);
      graph.insertLeftNode(v, u);
      currentSolution += graph.getCrossing(u, v);
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
      if (graph.getFixedPosition()[endIndex] != v) {
        graph.setFixedPosition(v, endIndex);
        if (undo) {
          undo->addSetPositionUndo(endIndex);
        }
      }
    }
  }
}

// If {u, v} is an incomparable pair in which, u and v are comparable with all other nodes, with
// c(u, v) <= c(v, u) , then commit u < v, and do the parameter accounting.
template <class Graph, class Undo>
bool rrlo2(Graph& graph, typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
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
void rr3(Graph& graph, typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
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
void rr2(Graph& graph, typename Graph::CrossingCountType& currentSolution) {
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
             typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
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
std::tuple<bool, typename Graph::NodeType, typename Graph::NodeType> IJBiggerThenFour(
    const Graph& graph) {
  using NodeType = typename Graph::NodeType;
  using CrossingCountType = typename Graph::CrossingCountType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      CrossingCountType secondFirstcrossingValue = graph.getCrossing(secondNode, firstNode);
      if (FirstSecondcrossingValue + secondFirstcrossingValue >= 4) {
        return std::make_tuple(true, firstNode, secondNode);
      }
    }
  }
  return std::make_tuple(false, 0, 0);
}

// Check if there is an there is a dependent 2/1 pattern {u, v} with c(u, v) + c(v, u) = 3
template <class Graph>
std::tuple<bool, typename Graph::NodeType, typename Graph::NodeType> IJEqualToThree(
    const Graph& graph) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      if (FirstSecondcrossingValue == 2) {
        return std::make_tuple(true, firstNode, secondNode);
      }
    }
  }
  return std::make_tuple(false, 0, 0);
}

// Check if there is a 1/1 pattern {u, v}
template <class Graph>
std::tuple<bool, typename Graph::NodeType, typename Graph::NodeType> IJEqualToTwo(
    const Graph& graph) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
      // Because we do not have any more IJBiggerThanFour and IJEqualToThree, and we delete the
      // crossing entries for each order that we set. All the crossings that are left are in the
      // form IJEqualToTwo.
      return std::make_tuple(true, firstNode, secondNode);
    }
  }
  return std::make_tuple(false, 0, 0);
}
template <class Graph, class Undo>
void algorithmStep(Graph& graph, typename Graph::CrossingCountType currentSolution, bool isInitStep,
                   typename Graph::NodeType leftNode, typename Graph::NodeType rightNode,
                   typename Graph::CrossingCountType& bestSolution,
                   std::vector<typename Graph::NodeType>& bestOrder) {
  using NodeType = typename Graph::NodeType;
  Undo undo;
  std::cout << "went in algorithmStep " << std::endl;

  if (isInitStep) {
    parameterAccounting<Graph, Undo>(graph, leftNode, rightNode, currentSolution, &undo);
  }
  bool didChangeRrlo2 = true;
  bool didChangeRrLarge = true;
  while (didChangeRrlo2 || didChangeRrLarge) {
    didChangeRrlo2 = rrlo2(graph, currentSolution, &undo);
    didChangeRrLarge = rrLarge(graph, bestSolution - currentSolution, currentSolution, &undo);
    rrlo1(graph, &undo);
  }
  if (bestSolution <= currentSolution) {
    graph.doUndo(undo);
    return;
  }
  std::tuple<bool, NodeType, NodeType> tupelBiggerThenFour = IJBiggerThenFour(graph);
  std::tuple<bool, NodeType, NodeType> EqualToThree = IJEqualToThree(graph);
  std::tuple<bool, NodeType, NodeType> EqualToTwo = IJEqualToTwo(graph);
  if (std::get<0>(tupelBiggerThenFour)) {
    NodeType u = std::get<1>(tupelBiggerThenFour);
    NodeType v = std::get<2>(tupelBiggerThenFour);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u, bestSolution, bestOrder);
  } else if (std::get<0>(EqualToThree)) {
    NodeType u = std::get<1>(EqualToThree);
    NodeType v = std::get<2>(EqualToThree);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u, bestSolution, bestOrder);
  } else if (std::get<0>(EqualToTwo)) {
    NodeType u = std::get<1>(EqualToTwo);
    NodeType v = std::get<2>(EqualToTwo);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    return;
  }
  bestSolution = currentSolution;
  bestOrder = graph.getFixedPosition();
  graph.doUndo(undo);
  return;
}

template <class Graph, class Undo>
std::tuple<typename Graph::CrossingCountType, std::vector<typename Graph::NodeType>> algorithm(
    Graph& graph) {
  using CrossingCountType = typename Graph::CrossingCountType;
  using NodeType = typename Graph::NodeType;

  // holds the crossing number of the best solution so far
  CrossingCountType bestSolution = INT_MAX;
  // holds the order of the solution best so far
  std::vector<NodeType> bestOrder;

  computeCrossingSums<Graph, Undo>(graph);
  CrossingCountType currentSolution = 0;
  rr2<Graph, Undo>(graph, currentSolution);
  bool didChangeRrlo2 = true;

  while (didChangeRrlo2) {
    didChangeRrlo2 = rrlo2<Graph, Undo>(graph, currentSolution);
    rrlo1<Graph, Undo>(graph);
  }

  rr3<Graph, Undo>(graph, currentSolution);
  algorithmStep<Graph, Undo>(graph, currentSolution, false, 0, 0, bestSolution, bestOrder);
  return std::make_tuple(bestSolution, bestOrder);
}