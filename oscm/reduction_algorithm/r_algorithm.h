#pragma once

#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/r_graph.h"

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
  static void rrlo2(RGraph<NodeType, CrossingCountType>& graph) {
    NodeType n = graph.getFreeNodesSize();
    for (NodeType u = 0; u < n; ++u) {
      if (graph.getLeftNodes(u).size() + graph.getRightNodes(u).size() == n - 2) {
        for (NodeType v = u + 1; v < n; ++v) {
          if (graph.getLeftNodes(v).size() + graph.getRightNodes(v).size() == n - 2) {
            if (graph.getLeftNodes(v).size() == graph.getLeftNodes(u).size()) {
              if (graph.getCrossing(u, v) <= graph.getCrossing(v, u)) {
                graph.parameterAccounting(u, v);
              } else {
                graph.parameterAccounting(v, u);
              }
            }
          }
        }
      }
    }
  }

  // if we have c(u, v) = 1 and c(v, u) = 2 with d(u) == 2  d(v) == 2 then commit u < v and do
  // parameter accounting
  static void rr3(RGraph<NodeType, CrossingCountType>& graph) {
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
      graph.parameterAccounting(pair.first, pair.second);
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

  static void rr2(RGraph<NodeType, CrossingCountType>& graph) {
    NodeType n = graph.getFreeNodesSize();
    for (NodeType u = 0; u < n; ++u) {
      for (NodeType v = u + 1; v < n; ++v) {
        if (areVectorsEqual(graph.getFreeNodeNeighbours(u), graph.getFreeNodeNeighbours(v))) {
          graph.parameterAccounting(u, v);
        }
      }
    }
  }
};