#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/reduction_graph.h"
#include "ds/reduction_graph/undo.h"

template <typename NodeType, typename CrossingCountType>
class ReductionAlgorithm {
 public:
  //  if a free node v is comparable  with all other free nodes, then put v in its right fixed
  //  position
  static void rrlo1(ReductionGraph<NodeType, CrossingCountType>& graph,
                    Undo<NodeType, CrossingCountType> undo) {
    for (NodeType v = 0; v < graph.getFreeNodesSize(); ++v) {
      if (graph.getNodeCrossing(v).size() == 0) {
        NodeType endIndex = graph.getLeftNodes(v).size();
        graph.setFixedPosition(v, endIndex);
        undo.addSetPositionUndo(endIndex);
      }
    }
  }

  // If {u, v} is an incomparable pair in which, u and v are comparable with all other nodes, with
  // c(u, v) <= c(v, u) , then commit u < v, and do the parameter accounting.
  static bool rrlo2(ReductionGraph<NodeType, CrossingCountType>& graph,
                    CrossingCountType* currentSolution) {
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
    return didChange;
  }

  // if we have c(u, v) = 1 and c(v, u) = 2 with d(u) == 2  d(v) == 2 then commit u < v and do
  // parameter accounting
  static void rr3(ReductionGraph<NodeType, CrossingCountType>& graph,
                  CrossingCountType* currentSolution) {
    // Create a list of pairs to be modified
    std::vector<std::pair<NodeType, NodeType>> pairsToModify;

    for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
      for (const auto& pair : graph.getNodeCrossing(u)) {
        NodeType v = pair.first;
        CrossingCountType crossingValue = pair.second;

        if (crossingValue == 1 && graph.getNodeCrossing(v).at(u) == 2 &&
            graph.getFreeNodeNeighboursSize(u) == 2 && graph.getFreeNodeNeighboursSize(v) == 2) {
          pairsToModify.emplace_back(u, v);
        }
      }
    }

    // Modify the pairs outside the loop
    for (const auto& pair : pairsToModify) {
      graph.parameterAccounting(pair.first, pair.second, currentSolution);
    }
  }

  // For each pair of free nodes u, v with N(u) = N(v),(arbitrarily) commit a < b, and do parameter
  // accounting.
  static void rr2(ReductionGraph<NodeType, CrossingCountType>& graph,
                  CrossingCountType* currentSolution) {
    NodeType n = graph.getFreeNodesSize();
    for (NodeType u = 0; u < n; ++u) {
      for (NodeType v = u + 1; v < n; ++v) {
        if (std::equal(graph.getFreeNodeNeighbours(u).begin(), graph.getFreeNodeNeighbours(u).end(),
                       graph.getFreeNodeNeighbours(v).begin(),
                       graph.getFreeNodeNeighbours(v).end())) {
          graph.parameterAccounting(u, v, currentSolution);
        }
      }
    }
  }

  // If c(u, v) > k, then commit v < u and do the parameter accounting.
  static bool rrLarge(ReductionGraph<NodeType, CrossingCountType>& graph,
                      CrossingCountType crossingsLeft, CrossingCountType* currentSolution) {
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
    for (auto [v, crossingValue] : pairsToModify) {
      graph.parameterAccounting(v, crossingValue, currentSolution);
    }
    return didChange;
  }

  // Check if there is  an incomparable i/j pattern {u, v} with i + j >= 4,
  static bool IJBiggerThenFour(const ReductionGraph<NodeType, CrossingCountType>& graph,
                               NodeType* u, NodeType* v) {
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
  static bool IJEqualToThree(const ReductionGraph<NodeType, CrossingCountType>& graph, NodeType* u,
                             NodeType* v) {
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
  static bool IJEqualToTwo(const ReductionGraph<NodeType, CrossingCountType>& graph, NodeType* u,
                           NodeType* v) {
    for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
      for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
        *u = firstNode;
        *v = secondNode;
        return true;
      }
    }
    return false;
  }

  static void algorithmStep(ReductionGraph<NodeType, CrossingCountType>& graph,
                            CrossingCountType currentSolution, bool isInitStep, NodeType leftNode,
                            NodeType rightNode) {
    Undo undo = Undo<NodeType, CrossingCountType>();
    if (isInitStep) {
      graph.parameterAccounting(leftNode, rightNode, &currentSolution, &undo);
    }
    bool didChangeRrlo2 = true;
    bool didChangeRrLarge = true;
    while (didChangeRrlo2 || didChangeRrLarge) {
      didChangeRrlo2 = rrlo2(graph, &currentSolution);
      didChangeRrLarge =
          rrLarge(graph, graph.getBestSolution() - currentSolution, &currentSolution);
      rrlo1(graph, undo);
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
    graph.setBestOrder(graph.getFixedPosition());
    graph.doUndo(undo);
    return;
  }

  static void algorithm(ReductionGraph<NodeType, CrossingCountType>& graph) {
    CrossingCountType currentSolution = 0;
    rr2(graph, &currentSolution);
    rr3(graph, &currentSolution);
    algorithmStep(graph, currentSolution, false, 0, 0);
  }
};