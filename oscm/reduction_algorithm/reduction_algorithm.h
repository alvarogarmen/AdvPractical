#pragma once

#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "ds/reduction_graph/reduction_graph.h"
#include "ds/reduction_graph/undo_algorithm_step.h"
#include "oscm/reduction_algorithm/compute_crossings.h"

namespace reductionalgorithms {
std::map<std::string, std::chrono::duration<double, std::milli>> functionDurations;

void printTotalDurations() {
  // Create a vector of pairs from the map
  std::vector<std::pair<std::string, std::chrono::duration<double, std::milli>>> durationVector(
      functionDurations.begin(), functionDurations.end());

  // Sort the vector by the duration
  std::sort(durationVector.begin(), durationVector.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

  // Print the sorted durations
  for (const auto& entry : durationVector) {
    std::cout << entry.first << " " << entry.second.count() << std::endl;
  }
}

void updateDuration(const std::string& functionName,
                    std::chrono::high_resolution_clock::time_point start) {
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration = end - start;
  functionDurations[functionName] += duration;
}

/*
 * For each decision u < v that we make adjust the leftRightSet of u and v
 * Add all the transitiv decisions that follow the u < v decision
 */
template <class Graph, class Undo>
void parameterAccounting(Graph& graph, typename Graph::NodeType u, typename Graph::NodeType v,
                         typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  if (u != v) {
    if (!graph.getRightNodesBit(u, v) && !graph.getLeftNodesBit(u, v)) {
      // if (graph.getNodeCrossing(u).find(v) != graph.getNodeCrossing(u).end()) {
      graph.insertRightNode(u, v);
      graph.insertLeftNode(v, u);
      // auto start = std::chrono::high_resolution_clock::now();
      currentSolution += graph.getCrossing(u, v);
      // updateDuration("Hash", start);

      if (undo) {
        undo->addParameterAccountingUndo(u, v, graph.getCrossing(u, v), graph.getCrossing(v, u));
      }
      // start = std::chrono::high_resolution_clock::now();
      //  graph.deleteCrossings(u, v);
      // updateDuration("Hash", start);

      for (auto smallerThanU : graph.getLeftNodes(u)) {
        parameterAccounting<Graph, Undo>(graph, smallerThanU, v, currentSolution, undo);
        for (auto biggerThanV : graph.getRightNodes(v)) {
          parameterAccounting<Graph, Undo>(graph, smallerThanU, biggerThanV, currentSolution, undo);
        }
      }
      for (auto biggerThanV : graph.getRightNodes(v)) {
        parameterAccounting<Graph, Undo>(graph, u, biggerThanV, currentSolution, undo);
      }
    }
  }
}

//  For each pair of free nodes a, b  that forms a 0/j pattern with j > 0, commit a < b.
template <class Graph, class Undo>
void rr1(Graph& graph, typename Graph::CrossingCountType& currentSolution) {
  using NodeType = typename Graph::NodeType;
  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    for (NodeType v = u + 1; v < graph.getFreeNodesSize(); ++v) {
      if ((!graph.getRightNodesBit(u, v)) && (!graph.getLeftNodesBit(u, v))) {
        if (graph.getCrossing(u, v) == 0) {
          parameterAccounting<Graph, Undo>(graph, u, v, currentSolution);
        } else if (graph.getCrossing(v, u) == 0) {
          parameterAccounting<Graph, Undo>(graph, v, u, currentSolution);
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
    auto uCrossings = graph.getNodeCrossing(u);
    for (size_t i = 0; i < uCrossings.size(); ++i) {
      if ((graph.getCrossing(u, uCrossings[i]) == 1) &&
          (graph.getCrossing(uCrossings[i], u) == 2) && (graph.getFreeNodeNeighboursSize(u) == 2) &&
          (graph.getFreeNodeNeighboursSize(uCrossings[i]) == 2)) {
        pairsToModify.emplace_back(u, uCrossings[i]);
      }
    }
    /*
    for (auto [v, crossingValue] : graph.getNodeCrossing(u)) {
      if (crossingValue == 1 && graph.getNodeCrossing(v).at(u) == 2 &&
          graph.getFreeNodeNeighboursSize(u) == 2 && graph.getFreeNodeNeighboursSize(v) == 2) {
        pairsToModify.emplace_back(u, v);
      }
    }*/
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
      if (graph.getNeighbourhoodHash(u) == graph.getNeighbourhoodHash(v)) {
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
  using CrossingCountType = typename Graph::CrossingCountType;

  auto pairsToModify = graph.getHeapTop(crossingsLeft);

  bool didChange = false;

  if (pairsToModify.size() > 0) {
    didChange = true;
    // std::cout << pairsToModify.size() << std::endl;
  }

  /*
  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    auto unsetNodesOfU = graph.getUnsetNodes(u);
    for (size_t i = 0; i < unsetNodesOfU.size(); ++i) {
      NodeType v = unsetNodesOfU[i];
      if (v != u) {
        if (graph.getCrossing(u, v) > crossingsLeft) {
          if (u < v) {
            pairsToModify.emplace_back(v, u);
            didChange = true;
          } else if (graph.getCrossing(v, u) < crossingsLeft) {
            pairsToModify.emplace_back(v, u);
            didChange = true;
          }
        }
      }
    }
  }
  * / /*
for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
for (auto [v, crossingValue] : graph.getNodeCrossing(u)) {
 if (crossingValue > crossingsLeft) {
   if (u < v) {
     pairsToModify.emplace_back(v, u);
     didChange = true;
   } else if (graph.getNodeCrossing(v).at(u) < crossingsLeft) {
     pairsToModify.emplace_back(v, u);
     didChange = true;
   }
 }
}
}*/
  for (auto [u, v, crossing] : pairsToModify) {
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
    auto unsetNodesOfU = graph.getUnsetNodes(firstNode);

    for (size_t i = 0; i < unsetNodesOfU.size(); ++i) {
      NodeType secondNode = unsetNodesOfU[i];
      auto FirstSecondcrossingValue = graph.getCrossing(firstNode, secondNode);
      auto secondFirstcrossingValue = graph.getCrossing(secondNode, firstNode);
      if (FirstSecondcrossingValue + secondFirstcrossingValue >= 4) {
        return std::make_tuple(true, firstNode, secondNode);
      }
    } /*
     for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
       CrossingCountType secondFirstcrossingValue = graph.getCrossing(secondNode, firstNode);
       if (FirstSecondcrossingValue + secondFirstcrossingValue >= 4) {
         return std::make_tuple(true, firstNode, secondNode);
       }*/
  }
  return std::make_tuple(false, 0, 0);
}

// Check if there is an there is a dependent 2/1 pattern {u, v} with c(u, v) + c(v, u) = 3
template <class Graph>
std::tuple<bool, typename Graph::NodeType, typename Graph::NodeType> IJEqualToThree(
    const Graph& graph) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    auto unsetNodesOfU = graph.getUnsetNodes(firstNode);

    for (size_t i = 0; i < unsetNodesOfU.size(); ++i) {
      NodeType secondNode = unsetNodesOfU[i];
      auto FirstSecondcrossingValue = graph.getCrossing(firstNode, secondNode);
      if (FirstSecondcrossingValue == 2) {
        return std::make_tuple(true, firstNode, secondNode);
      }
    }
    /*
        for (auto [secondNode, FirstSecondcrossingValue] : graph.getNodeCrossing(firstNode)) {
          if (FirstSecondcrossingValue == 2) {
            return std::make_tuple(true, firstNode, secondNode);
          }
        }
      }*/
    return std::make_tuple(false, 0, 0);
  }
}

// Check if there is a 1/1 pattern {u, v}
template <class Graph>
std::tuple<bool, typename Graph::NodeType, typename Graph::NodeType> IJEqualToTwo(
    const Graph& graph) {
  using NodeType = typename Graph::NodeType;
  for (NodeType firstNode = 0; firstNode < graph.getFreeNodesSize(); ++firstNode) {
    auto unsetNodesOfU = graph.getUnsetNodes(firstNode);

    for (size_t i = 0; i < unsetNodesOfU.size(); ++i) {
      // Because we do not have any more IJBiggerThanFour and IJEqualToThree, and we delete
      // the crossing entries for each order that we set. All the crossings that are left are
      // in the form IJEqualToTwo.
      return std::make_tuple(true, firstNode, unsetNodesOfU[i]);
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
  if (isInitStep) {
    // auto start = std::chrono::high_resolution_clock::now();
    parameterAccounting<Graph, Undo>(graph, leftNode, rightNode, currentSolution, &undo);
    // updateDuration("parameterAccounting", start);
  }

  bool didChangeRrlo2 = true;
  bool didChangeRrLarge = true;

  while (didChangeRrlo2 || didChangeRrLarge) {
    // auto start = std::chrono::high_resolution_clock::now();
    didChangeRrlo2 = rrlo2(graph, currentSolution, &undo);
    // updateDuration("rrlo2", start);
    // start = std::chrono::high_resolution_clock::now();
    didChangeRrLarge = rrLarge(graph, bestSolution - currentSolution, currentSolution, &undo);
    // updateDuration("rrLarge", start);

    // start = std::chrono::high_resolution_clock::now();
    rrlo1(graph, &undo);
    // updateDuration("rrlo1", start);
  }

  if (bestSolution <= currentSolution) {
    graph.doUndo(undo);
    return;
  }

  // auto start = std::chrono::high_resolution_clock::now();
  std::tuple<bool, NodeType, NodeType> tupelBiggerThenFour = IJBiggerThenFour(graph);
  // updateDuration("IJBiggerThenFour", start);

  // start = std::chrono::high_resolution_clock::now();
  std::tuple<bool, NodeType, NodeType> EqualToThree = IJEqualToThree(graph);
  // updateDuration("IJEqualToThree", start);

  // start = std::chrono::high_resolution_clock::now();
  std::tuple<bool, NodeType, NodeType> EqualToTwo = IJEqualToTwo(graph);
  // updateDuration("IJEqualToTwo", start);

  if (std::get<0>(tupelBiggerThenFour)) {
    NodeType u = std::get<1>(tupelBiggerThenFour);
    NodeType v = std::get<2>(tupelBiggerThenFour);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u, bestSolution, bestOrder);
    graph.doUndo(undo);
    return;
  } else if (std::get<0>(EqualToThree)) {
    NodeType u = std::get<1>(EqualToThree);
    NodeType v = std::get<2>(EqualToThree);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, v, u, bestSolution, bestOrder);
    graph.doUndo(undo);
    return;
  } else if (std::get<0>(EqualToTwo)) {
    NodeType u = std::get<1>(EqualToTwo);
    NodeType v = std::get<2>(EqualToTwo);
    algorithmStep<Graph, Undo>(graph, currentSolution, true, u, v, bestSolution, bestOrder);
    graph.doUndo(undo);
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
  CrossingCountType bestSolution = 1482;
  // holds the order of the solution best so far
  std::vector<NodeType> bestOrder;

  // auto start = std::chrono::high_resolution_clock::now();
  computeCrossingSums<Graph, Undo>(graph);
  // updateDuration("computeCrossingSums", start);
  CrossingCountType currentSolution = 0;
  // start = std::chrono::high_resolution_clock::now();
  rr1<Graph, Undo>(graph, currentSolution);
  // updateDuration("rr1", start);
  // start = std::chrono::high_resolution_clock::now();
  rr2<Graph, Undo>(graph, currentSolution);
  // updateDuration("rr2", start);
  bool didChangeRrlo2 = true;

  while (didChangeRrlo2) {
    // start = std::chrono::high_resolution_clock::now();
    didChangeRrlo2 = rrlo2<Graph, Undo>(graph, currentSolution);
    // updateDuration("rrlo2", start);
    // start = std::chrono::high_resolution_clock::now();
    rrlo1<Graph, Undo>(graph);
    // updateDuration("rrlo1", start);
  }

  // start = std::chrono::high_resolution_clock::now();
  rr3<Graph, Undo>(graph, currentSolution);
  // updateDuration("rr3", start);
  graph.createHeap();
  // start = std::chrono::high_resolution_clock::now();
  algorithmStep<Graph, Undo>(graph, currentSolution, false, 0, 0, bestSolution, bestOrder);
  // updateDuration("algorithmStep", start);
  graph.setFixedPositions(bestOrder);
  // printTotalDurations();
  return std::make_tuple(bestSolution, bestOrder);
}
}  // namespace reductionalgorithms