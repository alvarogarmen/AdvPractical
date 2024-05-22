#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

// #include "absl/container/node_hash_map.h"
//  #include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "ds/reduction_graph/bit_vector.h"
#include "ds/reduction_graph/max_heap.h"
#include "ds/reduction_graph/undo_algorithm_step.h"
// #include "external_hash/unordered_dense/include/ankerl/unordered_dense.h"

template <typename NT, typename CCT>
class ReductionGraph {
  // using CrossingMap = ankerl::unordered_dense::map<NT, CCT>;
  //  for each free node holds its neighbours
  std::vector<std::vector<NT>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NT>> fixedNodes;
  // holds the end position of free nodes
  std::vector<NT> fixedPosition;
  // for each free node u, the first place hold free node that are knowned to be positioned to the
  // left of u
  // and the second hold the free nodes that are knowned to be positioned to its right
  std::vector<std::array<BitVector<NT>, 2>> leftRightSet;
  // save the crossing number that accur between two free nodes (u, v) ,that do not have < order
  // yet.
  // saves the crossing number of u and v assuming u is placed before v
  // std::vector<absl::node_hash_map<NT, CCT>> crossings;
  // std::vector<std::map<NT, CCT>> crossings;
  // std::vector<CrossingMap> crossings;
  std::vector<std::vector<unsigned short>> crossings;
  //  stores the hash values of neighbourhood
  std::vector<NT> neighbourhoodHash;
  // Max heap for the crossing of each pair of nodes
  MaxHeap<NT, CCT> maxHeap;

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
        crossings(freeNodes.size()),
        neighbourhoodHash(freeNodes.size()) {
    for (auto& bitVectorArray :
         leftRightSet) {  // Initialize each BitVector in the array to be the right size
      for (auto& bitVector : bitVectorArray) {
        bitVector.resize(freeNodes.size());
      }
    }
  }

  ReductionGraph(NodeType numFreeNodes, NodeType numFixedNodes)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size(), std::vector<unsigned short>(freeNodes.size())),
        neighbourhoodHash(freeNodes.size()) {
    for (auto& bitVectorArray : leftRightSet) {
      // Initialize each BitVector in the array to be the right size
      for (auto& bitVector : bitVectorArray) {
        bitVector.resize(freeNodes.size());
      }
    }
  }

  ReductionGraph(NodeType numFixedNodes, NodeType numFreeNodes, CrossingCountType numEdges)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFixedNodes, std::vector<NodeType>(0)),
        fixedPosition(freeNodes.size()),
        leftRightSet(freeNodes.size()),
        crossings(freeNodes.size(), std::vector<unsigned short>(freeNodes.size())),
        neighbourhoodHash(freeNodes.size()) {
    for (auto& bitVectorArray : leftRightSet) {
      // Initialize each BitVector in the array to be the right size
      for (auto& bitVector : bitVectorArray) {
        bitVector.resize(freeNodes.size());
      }
    }
  }

  void addEdge(NodeType source,
               NodeType target) {  // where source is the freeNode and target is the fixedNode
    freeNodes[source].push_back(target);
    fixedNodes[target].push_back(source);
    neighbourhoodHash[source] += target * target;
    return;
  }
  NodeType getNeighbourhoodHash(NodeType index) { return neighbourhoodHash[index]; }
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

  const auto getNodeCrossing(NodeType u) const {
    return leftRightSet[u][0].findCommonUnsetBits(leftRightSet[u][1], freeNodes.size());
  }

  const auto& getCrossing(NodeType u, NodeType v) const { return crossings[u][v]; }

  const auto getLeftNodes(NodeType u) const { return leftRightSet[u][0].findSetBits(); }
  const auto getUnsetNodes(NodeType u) const {
    return leftRightSet[u][0].findCommonUnsetBits(leftRightSet[u][1], freeNodes.size());
  }

  void insertRightNode(NodeType u, NodeType v) { leftRightSet[u][1].insert(v); }

  void insertLeftNode(NodeType u, NodeType v) { leftRightSet[u][0].insert(v); }

  const auto getRightNodesBit(NodeType u, NodeType v) const { return leftRightSet[u][1].find(v); }
  const auto getLeftNodesBit(NodeType u, NodeType v) const { return leftRightSet[u][0].find(v); }

  const auto getRightNodes(NodeType u) const { return leftRightSet[u][1].findSetBits(); }

  const auto& getFixedPosition() const { return fixedPosition; }

  void setFixedPosition(NodeType u, NodeType index) { fixedPosition[index] = u; }

  void setFixedPositions(const std::vector<NodeType>& bestOrder) { fixedPosition = bestOrder; }

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
  void addCrossingToHeap(NodeType leftNode, NodeType rightNode, CrossingCountType crossingSum) {
    maxHeap.push(leftNode, rightNode, crossingSum);
  }
  const auto getHeapSize() { return maxHeap.getSize(); }

  const auto getHeapTop(CrossingCountType crossingsLeft) {
    std::vector<std::tuple<NodeType, NodeType, CrossingCountType>> pairsToModify;
    bool isTopINHeap = false;
    auto topElement = maxHeap.top();
    NodeType firstNode = std::get<0>(topElement);
    NodeType secondNode = std::get<1>(topElement);
    while (!isTopINHeap) {
      if (std::get<2>(topElement) < crossingsLeft) {
        isTopINHeap = true;
      } else {
        if (!getRightNodesBit(firstNode, secondNode) && !getLeftNodesBit(firstNode, secondNode)) {
          pairsToModify.push_back(topElement);
        }
        maxHeap.pop();
        topElement = maxHeap.top();
      }
      return pairsToModify;
    }
  }

  void createHeap() {
    for (size_t firstNode = 0; firstNode < freeNodes.size(); ++firstNode) {
      auto unsetNodesOfU = getUnsetNodes(firstNode);
      for (size_t i = 0; i < unsetNodesOfU.size(); ++i) {
        NodeType secondNode = unsetNodesOfU[i];
        auto FirstSecondcrossingValue = getCrossing(firstNode, secondNode);
        addCrossingToHeap(firstNode, secondNode, FirstSecondcrossingValue);
      }
    }
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
      addCrossingToHeap(operation.leftNode, operation.rightNode, operation.leftRightCrossing);
      addCrossingToHeap(operation.rightNode, operation.leftNode, operation.rightLeftCrossing);
    }
    for (const auto& position : undo.getSetPositionUndo()) {
      setFixedPosition(0, position);
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
    for (size_t u = 0; u < freeNodes.size(); ++u) {
      for (size_t v = u + 1; v < freeNodes.size(); ++v) {
        crossingCount += computeUVcrossing(fixedPosition[u], fixedPosition[v]);
      }
    }
    return crossingCount;
  }
};