#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "undo.h"

template <typename NodeType, typename CrossingCountType>
class RGraph {
  // for each free node holds its neighbours
  std::vector<std::vector<NodeType>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NodeType>> fixedNodes;
  // holds the end position of free nodes
  std::vector<NodeType> fixedPosition;
  // for each free node u, the first place hold free node that are knowned to be positioned to the
  // left of u
  // and the second hold the free nodes that are knowned to be positioned to its right
  std::vector<std::array<std::set<NodeType>, 2>> leftRightSet;
  // save the crossing number that accur between two free nodes (u, v) ,that do not have < order
  // yet.
  // saves the crossing number of u and v assuming u is placed before v
  std::vector<std::map<NodeType, CrossingCountType>> crossings;
  // for each new branch save all the changes in order to do back tracking
  std::vector<NodeType> undo;
  // holds the crossing number of the solution so far
  CrossingCountType bestSolution;
  // holds the crossing number of the current selution
  CrossingCountType currentSolution;

 public:
  RGraph(const std::vector<std::vector<NodeType>>& freeNodes,
         const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        fixedPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightSet(std::vector<std::array<std::set<NodeType>, 2>>(freeNodes.size())),
        crossings(std::vector<std::map<NodeType, CrossingCountType>>(freeNodes.size())) {
    currentSolution = 0;
    computeCrossingSums();
  }

  RGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        fixedNodes(std::vector<std::vector<NodeType>>(numFreeNodes, std::vector<NodeType>(0))),
        fixedPosition(std::vector<NodeType>(freeNodes.size())),
        leftRightSet(std::vector<std::array<std::set<NodeType>, 2>>(freeNodes.size())),
        crossings(std::vector<std::map<NodeType, CrossingCountType>>(freeNodes.size())) {
    currentSolution = 0;
    computeCrossingSums();
  }

  NodeType getFixedNodesSize() const { return fixedNodes.size(); }

  NodeType getFixedNodeNeighboursSize(NodeType nodeID) const { return fixedNodes[nodeID].size(); }

  const auto& getFixedNodeNeighbours(NodeType fixedNodeID) const { return fixedNodes[fixedNodeID]; }

  NodeType getFreeNodesSize() const { return freeNodes.size(); }

  NodeType getFreeNodeNeighboursSize(NodeType nodeID) const { return freeNodes[nodeID].size(); }

  const auto& getFreeNodeNeighbours(NodeType freeNodeID) const { return fixedNodes[freeNodeID]; }

  const auto& getNodeCrossing(NodeType u) const { return crossings[u]; }

  const auto& getCrossing(NodeType u, NodeType v) const { return crossings[u].at(v); }

  const auto& getLeftNodes(NodeType u) const { return leftRightSet[u][0]; }

  const auto& getRightNodes(NodeType u) const { return leftRightSet[u][1]; }

  void setLeftNodes(NodeType u, std::set<NodeType> leftNodes) { leftRightSet[u][0] = leftNodes; }

  void setRightNodes(NodeType u, std::set<NodeType> rightNodes) { leftRightSet[u][1] = rightNodes; }

  const auto& getFixedPosition() const { return fixedPosition; }

  void setFixedPosition(NodeType u, NodeType index) { fixedPosition[index] = u; }

  const auto& getCurrentSolution() const { return currentSolution; }

  void setCurrentSolution(CrossingCountType newCurrent) { currentSolution = newCurrent; }

  const auto& getBestSolution() const { return bestSolution; }

  void setBestSolution(CrossingCountType newBest) { bestSolution = newBest; }

  void setCrossings(std::vector<std::map<NodeType, CrossingCountType>> m) { crossings = m; }

  void setFreeNodes(std::vector<std::vector<NodeType>> newfreeNodes) { freeNodes = newfreeNodes; }

  void deleteLeftNode(NodeType node, NodeType leftNode) { leftRightSet[node][0].erase(leftNode); }

  void deleteRightNode(NodeType node, NodeType rightNode) {
    leftRightSet[node][1].erase(rightNode);
  }

  void addCrossing(NodeType leftNode, NodeType rightNode, CrossingCountType crossingSum) {
    crossings[leftNode][rightNode] = crossingSum;
  }

  void doUndo(Undo<NodeType, CrossingCountType>& undo) {
    for (const auto& operation : undo.getParameterAccountingUndo()) {
      deleteLeftNode(operation.rightNode, operation.leftNode);
      deleteRightNode(operation.leftNode, operation.rightNode);
      addCrossing(operation.leftNode, operation.rightNode, operation.leftRightCrossing);
      addCrossing(operation.rightNode, operation.leftNode, operation.rightLeftCrossing);
    }
    for (const auto& position : undo.getSetPositionUndo()) {
      setFixedPosition(0, position);
    }
  }

  void clearLeftRightSet() {
    for (NodeType u = 0; u < freeNodes.size(); ++u) {
      leftRightSet[u][0].clear();
      leftRightSet[u][1].clear();
    }
  }

  void parameterAccounting(NodeType u, NodeType v,
                           Undo<NodeType, CrossingCountType>* undo = nullptr) {
    if (u != v) {
      if (leftRightSet[u][1].find(v) == leftRightSet[u][1].end()) {
        leftRightSet[u][1].insert(v);
        leftRightSet[v][0].insert(u);
        currentSolution += crossings[u][v];
        if (undo) {
          (*undo).addParameterAccountingUndo(u, v, crossings[u][v], crossings[v][u]);
        }
        crossings[u].erase(v);
        crossings[v].erase(u);
        for (const auto& smallerThanU : leftRightSet[u][0]) {
          parameterAccounting(smallerThanU, v, undo);
          for (const auto& biggerThanV : leftRightSet[v][1]) {
            parameterAccounting(smallerThanU, biggerThanV, undo);
            parameterAccounting(u, biggerThanV, undo);
          }
        }
      }
    }
  }

  void addEdge(NodeType freeNode, NodeType fixedNode) {
    freeNodes[freeNode].push_back(fixedNode);
    fixedNodes[fixedNode].push_back(freeNode);
  }

  // copmute the number of crossings created by the edges from two free nodes (u, v)
  // when u is to the left of v
  CrossingCountType const computeUVcrossing(NodeType u, NodeType v) {
    NodeType crossingSum = 0;
    for (const auto uNeighbour : freeNodes[u]) {
      for (const auto vNeighbour : freeNodes[v]) {
        crossingSum += vNeighbour < uNeighbour;
      }
    }
    return crossingSum;
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right

  void computeCrossingSums() {
    for (NodeType u = 0; u < freeNodes.size(); ++u) {
      for (NodeType v = u + 1; v < freeNodes.size(); ++v) {
        CrossingCountType crossingUV = computeUVcrossing(u, v);
        // reduction RR1: For each pair of vertices {u, v} ⊆ free nodes that forms a 0/j pattern
        // with j > 0, commit u < v
        if (crossingUV == 0) {
          parameterAccounting(u, v);
        } else {
          CrossingCountType crossingVU = computeUVcrossing(v, u);
          if (crossingVU == 0) {
            parameterAccounting(v, u);
          } else {
            crossings[u][v] = crossingUV;
            crossings[v][u] = crossingVU;
          }
        }
      }
    }
  }
};