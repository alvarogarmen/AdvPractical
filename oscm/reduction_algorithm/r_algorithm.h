#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/r_graph.h"
#include "ds/reduction_graph/undo.h"

template <typename NodeType, typename CrossingCountType>
class ReductionAlgorithm {
 public:
  //  if a free node v is comparable  with all other free nodes, then put v in its right fixed
  //  position
  static void rrlo1(RGraph<NodeType, CrossingCountType>& graph) {
    for (NodeType v = 0; v < graph.getFreeNodesSize(); ++v) {
      if (graph.getNodeCrossing(v).size() == 0) {
        NodeType endIndex = graph.getLeftNodes(v).size();
        graph.setFixedPosition(v, endIndex);
      }
    }
  }

  // If {u, v} is an incomparable pair in which, u and v are comparable with all other nodes, with
  // c(u, v) <= c(v, u) , then commit u < v, and do the parameter accounting.
  static bool rrlo2(RGraph<NodeType, CrossingCountType>& graph,
                    CrossingCountType& currentSolution) {
    NodeType n = graph.getFreeNodesSize();
    bool didChange = false;
    for (NodeType u = 0; u < n; ++u) {
      if (graph.getLeftNodes(u).size() + graph.getRightNodes(u).size() == n - 2) {
        for (NodeType v = u + 1; v < n; ++v) {
          if (graph.getLeftNodes(v).size() + graph.getRightNodes(v).size() == n - 2) {
            if (graph.getLeftNodes(v).size() == graph.getLeftNodes(u).size()) {
              if (graph.getCrossing(u, v) <= graph.getCrossing(v, u)) {
                graph.parameterAccounting(u, v, currentSolution);
                didChange = true;
              } else {
                graph.parameterAccounting(v, u, currentSolution);
                didChange = true;
              }
            }
          }
        }
      }
    }
    return didChange;
  }

  // if we have c(u, v) = 1 and c(v, u) = 2 with d(u) == 2  d(v) == 2 then commit u < v and do
  // parameter accounting
  static void rr3(RGraph<NodeType, CrossingCountType>& graph, CrossingCountType& currentSolution) {
    // Create a list of pairs to be modified
    std::vector<std::pair<NodeType, NodeType>> pairsToModify;

    for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
      for (const auto& pair : graph.getNodeCrossing(u)) {
        NodeType v = pair.first;
        CrossingCountType crossingValue = pair.second;

        if (crossingValue == 1) {
          if (graph.getNodeCrossing(v).at(u) == 2) {
            if (graph.getFreeNodeNeighboursSize(u) == 2 &&
                graph.getFreeNodeNeighboursSize(v) == 2) {
              pairsToModify.emplace_back(u, v);
            }
          }
        }
      }
    }

    // Modify the pairs outside the loop
    for (const auto& pair : pairsToModify) {
      graph.parameterAccounting(pair.first, pair.second, currentSolution);
    }
  }

  static bool areVectorsEqual(const std::vector<NodeType>& v1, const std::vector<NodeType>& v2) {
    if (v1.size() != v2.size()) {
      return false;
    }

    for (NodeType i = 0; i < v1.size(); ++i) {
      if (v1[i] != v2[i]) {
        return false;
      }
    }

    return true;
  }
  // For each pair of free nodes u, v with N(u) = N(v),(arbitrarily) commit a < b, and do parameter
  // accounting.
  static void rr2(RGraph<NodeType, CrossingCountType>& graph, CrossingCountType& currentSolution) {
    NodeType n = graph.getFreeNodesSize();
    for (NodeType u = 0; u < n; ++u) {
      for (NodeType v = u + 1; v < n; ++v) {
        if (areVectorsEqual(graph.getFreeNodeNeighbours(u), graph.getFreeNodeNeighbours(v))) {
          graph.parameterAccounting(u, v, currentSolution);
        }
      }
    }
  }

  // If c(u, v) > k, then commit v < u and do the parameter accounting.
  static bool rrLarge(RGraph<NodeType, CrossingCountType>& graph, CrossingCountType crossingsLeft,
                      CrossingCountType& currentSolution) {
    std::vector<std::pair<NodeType, NodeType>> pairsToModify;
    bool didChange = false;
    for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
      for (const auto& pair : graph.getNodeCrossing(u)) {
        NodeType v = pair.first;
        CrossingCountType crossingValue = pair.second;

        if (crossingValue > crossingsLeft) {
          pairsToModify.emplace_back(v, u);
          didChange = true;
        }
      }
    }
    for (const auto& pair : pairsToModify) {
      graph.parameterAccounting(pair.first, pair.second, currentSolution);
    }
    return didChange;
  }

  // Check if there is  an incomparable i/j pattern {u, v} with i + j >= 4,
  static bool IJBiggerThenFour(const RGraph<NodeType, CrossingCountType>& graph, NodeType* u,
                               NodeType* v) {
    for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
      for (const auto& pair : graph.getNodeCrossing(firstNode)) {
        NodeType secondNode = pair.first;
        CrossingCountType FirstSecondcrossingValue = pair.second;
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
  static bool IJEqualToThree(const RGraph<NodeType, CrossingCountType>& graph, NodeType* u,
                             NodeType* v) {
    for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
      for (const auto& pair : graph.getNodeCrossing(firstNode)) {
        NodeType secondNode = pair.first;
        CrossingCountType FirstSecondcrossingValue = pair.second;
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
  static bool IJEqualToTwo(const RGraph<NodeType, CrossingCountType>& graph, NodeType* u,
                           NodeType* v) {
    for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
      for (const auto& pair : graph.getNodeCrossing(firstNode)) {
        NodeType secondNode = pair.first;
        *u = firstNode;
        *v = secondNode;
        return true;
      }
    }
    return false;
  }

  static void algorithmStep(RGraph<NodeType, CrossingCountType>& graph,
                            CrossingCountType currentSolution, bool isInitStep, NodeType leftNode,
                            NodeType rightNode) {
    Undo undo = Undo<NodeType, CrossingCountType>();
    if (isInitStep) {
      graph.parameterAccounting(leftNode, rightNode, currentSolution, undo);
    }
    bool didChangeRrlo2 = true;
    bool didChangeRrLarge = true;
    while (didChangeRrlo2 || didChangeRrLarge) {
      didChangeRrlo2 = rrlo2(graph, currentSolution);
      didChangeRrLarge = rrLarge(graph, graph.getBestSolution() - currentSolution, currentSolution);
      rrlo1(graph);
    }
    if (graph.getBestSolution() <= currentSolution) {
      graph.doUndo(undo);
      return;
    }
    NodeType u;
    NodeType v;
    if (IJBiggerThenFour(graph, &u, &v)) {
      algorithmStep(graph, currentSolution, true, u, v);
      algorithmStep(graph, currentSolution, true, v, u);
    } else if (IJEqualToThree(graph, &u, &v)) {
      algorithmStep(graph, currentSolution, true, u, v);
      algorithmStep(graph, currentSolution, true, v, u);
    } else if (IJEqualToTwo(graph, &u, &v)) {
      algorithmStep(graph, currentSolution, true, u, v);
    }
    graph.setBestSolution(currentSolution);
    graph.setBestOrder(graph.fixedPosition);
    graph.doUndo(undo);
    return;
  }
};