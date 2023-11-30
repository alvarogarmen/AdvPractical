#pragma once
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "ds/reduction_graph/UndoAlgorithmStep.h"

template <typename NT, typename CCT>
class ReductionGraph {
  // for each free node holds its neighbours
  std::vector<std::vector<NT>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NT>> fixedNodes;
  // holds the end position of free nodes
  std::vector<NT> fixedPosition;
  // for each free node u, the first place hold free node that are knowned to be positioned to the
  // left of u
  // and the second hold the free nodes that are knowned to be positioned to its right
  std::vector<std::array<std::set<NT>, 2>> leftRightSet;
  // save the crossing number that accur between two free nodes (u, v) ,that do not have < order
  // yet.
  // saves the crossing number of u and v assuming u is placed before v
  std::vector<std::map<NT, CCT>> crossings;
  // holds the crossing number of the best solution so far
  CCT bestSolution;
  // holds the order of the solution best so far
  std::vector<NT> bestOrder;

 public:
  using NodeType = NT;
  using CrossingCountType = CCT;
  ReductionGraph(const std::vector<std::vector<NodeType>>& freeNodes,
                 const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size()) {
    computeCrossingSums();
  }

  ReductionGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size()) {
    computeCrossingSums();
  }

  NodeType getFixedNodesSize() const { return fixedNodes.size(); }

  NodeType getFixedNodeNeighboursSize(NodeType nodeID) const { return fixedNodes[nodeID].size(); }

  const auto& getFixedNodeNeighbours(NodeType fixedNodeID) const { return fixedNodes[fixedNodeID]; }

  NodeType getFreeNodesSize() const { return freeNodes.size(); }

  NodeType getFreeNodeNeighboursSize(NodeType nodeID) const { return freeNodes[nodeID].size(); }

  const auto& getFreeNodeNeighbours(NodeType freeNodeID) const { return freeNodes[freeNodeID]; }

  const auto& getNodeCrossing(NodeType u) const { return crossings[u]; }

  const auto& getCrossing(NodeType u, NodeType v) const { return crossings[u].at(v); }

  const auto& getLeftNodes(NodeType u) const { return leftRightSet[u][0]; }

  const auto& getRightNodes(NodeType u) const { return leftRightSet[u][1]; }

  void setLeftNodes(NodeType u, std::set<NodeType> leftNodes) { leftRightSet[u][0] = leftNodes; }

  void setRightNodes(NodeType u, std::set<NodeType> rightNodes) { leftRightSet[u][1] = rightNodes; }

  const auto& getFixedPosition() const { return fixedPosition; }

  void setFixedPosition(NodeType u, NodeType index) { fixedPosition[index] = u; }

  const auto& getBestSolution() const { return bestSolution; }

  const auto& getBestOrder() const { return bestOrder; }

  void setBestSolution(CrossingCountType newBest) { bestSolution = newBest; }

  void setCrossings(std::vector<std::map<NodeType, CrossingCountType>> m) { crossings = m; }

  void setFreeNodes(std::vector<std::vector<NodeType>> newfreeNodes) { freeNodes = newfreeNodes; }

  void deleteLeftNode(NodeType node, NodeType leftNode) { leftRightSet[node][0].erase(leftNode); }

  void deleteRightNode(NodeType node, NodeType rightNode) {
    leftRightSet[node][1].erase(rightNode);
  }

  void setBestOrder(const std::vector<NodeType>& bestOrderSoFar) { bestOrder = bestOrderSoFar; }

  void addCrossing(NodeType leftNode, NodeType rightNode, CrossingCountType crossingSum) {
    crossings[leftNode][rightNode] = crossingSum;
  }

  void doUndo(UndoAlgorithmStep<NodeType, CrossingCountType>& undo) {
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

  /*
   * For each decision u < v that we make adjust the leftRightSet of u and v
   * Add all the transitiv decisions that follow the u < v decision
   */
  void parameterAccounting(NodeType u, NodeType v, CrossingCountType* currentSolution,
                           UndoAlgorithmStep<NodeType, CrossingCountType>* undo = nullptr) {
    if (u != v) {
      if (leftRightSet[u][1].find(v) == leftRightSet[u][1].end()) {
        leftRightSet[u][1].insert(v);
        leftRightSet[v][0].insert(u);
        *currentSolution += crossings[u][v];
        if (undo) {
          undo->addParameterAccountingUndo(u, v, crossings[u][v], crossings[v][u]);
        }
        crossings[u].erase(v);
        crossings[v].erase(u);
        for (NodeType smallerThanU : leftRightSet[u][0]) {
          parameterAccounting(smallerThanU, v, currentSolution, undo);
          for (NodeType biggerThanV : leftRightSet[v][1]) {
            parameterAccounting(smallerThanU, biggerThanV, currentSolution, undo);
            parameterAccounting(u, biggerThanV, currentSolution, undo);
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
    CrossingCountType currentSolution = 0;
    for (NodeType u = 0; u < freeNodes.size(); ++u) {
      for (NodeType v = u + 1; v < freeNodes.size(); ++v) {
        CrossingCountType crossingUV = computeUVcrossing(u, v);
        // reduction RR1: For each pair of vertices {u, v} ⊆ free nodes that forms a 0/j pattern
        // with j > 0, commit u < v
        if (crossingUV == 0) {
          parameterAccounting(u, v, &currentSolution);
        } else {
          CrossingCountType crossingVU = computeUVcrossing(v, u);
          if (crossingVU == 0) {
            parameterAccounting(v, u, &currentSolution);
          } else {
            crossings[u][v] = crossingUV;
            crossings[v][u] = crossingVU;
          }
        }
      }
    }
  }
};