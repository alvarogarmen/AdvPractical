#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "ds/reduction_graph/undo_algorithm_step.h"

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

 public:
  NT currentNumNodes() { return getFreeNodesSize() + getFixedNodesSize(); }
  NT currentNumEdges() { return 0; }
  using EdgeType = NT;
  using NodeType = NT;
  using WeightType = NT;
  using CrossingCountType = CCT;
  ReductionGraph(const std::vector<std::vector<NodeType>>& freeNodes,
                 const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size()) {}

  ReductionGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size()) {}

  ReductionGraph(NodeType numFixedNodes, NodeType numFreeNodes, CrossingCountType numEdges)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFixedNodes, std::vector<NodeType>(0)),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size()) {}

  void addEdge(NodeType source,
               NodeType target) {  // where source is the freeNode and target is the fixedNode
    freeNodes[source].push_back(target);
    fixedNodes[target].push_back(source);
    return;
  }

  NodeType getEdge(NodeType source, NodeType index) { return freeNodes[source][index]; }
  /**
  Should be moved to io
  This function write the result order into a file.
  @param os The output stream.
  @param solution The solution order.
*/
  absl::Status writeResultsToFile(std::ostream& os, const std::vector<NodeType>& solution) {
    // Write each element of the vector to the output stream
    for (NodeType element : solution) {
      NodeType n = fixedNodes.size();
      element = element + n + 1;
      os << element << std::endl;

      // Check for errors after writing each element
      if (!os) {
        return absl::UnknownError(absl::StrCat("Error writing element to output stream"));
      }
    }
    return absl::OkStatus();
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

  void insertRightNode(NodeType u, NodeType v) { leftRightSet[u][1].insert(v); }

  void insertLeftNode(NodeType u, NodeType v) { leftRightSet[u][0].insert(v); }

  const auto& getRightNodes(NodeType u) const { return leftRightSet[u][1]; }

  void setLeftNodes(NodeType u, std::set<NodeType> leftNodes) { leftRightSet[u][0] = leftNodes; }

  void setRightNodes(NodeType u, std::set<NodeType> rightNodes) { leftRightSet[u][1] = rightNodes; }

  const auto& getFixedPosition() const { return fixedPosition; }

  void setFixedPosition(NodeType u, NodeType index) { fixedPosition[index] = u; }

  void setFixedPositions(const std::vector<NodeType>& bestOrder) { fixedPosition = bestOrder; }


  void setCrossings(const std::vector<std::map<NodeType, CrossingCountType>>& m) { crossings = m; }

  void setFreeNodes(const std::vector<std::vector<NodeType>>& newfreeNodes) {
    freeNodes = newfreeNodes;
  }

  void deleteLeftNode(NodeType node, NodeType leftNode) { leftRightSet[node][0].erase(leftNode); }

  void deleteRightNode(NodeType node, NodeType rightNode) {
    leftRightSet[node][1].erase(rightNode);
  }

  void addCrossing(NodeType leftNode, NodeType rightNode, CrossingCountType crossingSum) {
    crossings[leftNode][rightNode] = crossingSum;
  }

  void deleteCrossings(NodeType u, NodeType v) {
    crossings[u].erase(v);
    crossings[v].erase(u);
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

  CrossingCountType const computeUVcrossing(NodeType u, NodeType v) {
    NodeType crossingSum = 0;
    for (const auto uNeighbour : freeNodes[u]) {
      for (const auto vNeighbour : freeNodes[v]) {
        crossingSum += vNeighbour < uNeighbour;
      }
    }
    return crossingSum;
  }

CrossingCountType getCrossings() {
  CrossingCountType crossingCount = 0;
  for (NodeType u = 0; u < freeNodes.size(); ++u) {
    for (NodeType v = u + 1; v < freeNodes.size(); ++v) {
      crossingCount += computeUVcrossing(fixedPosition[u], fixedPosition[v]);
    }
  }
  return crossingCount;
}
};