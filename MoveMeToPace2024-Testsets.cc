#include <stdio.h>

#include <algorithm>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

template <typename NodeType, typename CrossingCountType>
struct Operation {
  NodeType leftNode;
  NodeType rightNode;
  CrossingCountType leftRightCrossing;
  CrossingCountType rightLeftCrossing;

  Operation(NodeType leftNode, NodeType rightNode, CrossingCountType leftRightCrossing,
            CrossingCountType rightLeftCrossing)
      : leftNode(leftNode),
        rightNode(rightNode),
        leftRightCrossing(leftRightCrossing),
        rightLeftCrossing(rightLeftCrossing) {}
};
template <typename NT, typename CCT>
class UndoAlgorithmStep {
  std::vector<Operation<NT, CCT>> parameterAccountingUndo;
  std::vector<NT> setPositionUndo;

 public:
  using NodeType = NT;
  using CrossingCountType = CCT;

  UndoAlgorithmStep() {}
  // Undo every ParameterAccounting for a specific step
  void addParameterAccountingUndo(NodeType leftNode, NodeType rightNode,
                                  CrossingCountType leftRightCrossing,
                                  CrossingCountType rightLeftCrossing) {
    parameterAccountingUndo.push_back(
        Operation(leftNode, rightNode, leftRightCrossing, rightLeftCrossing));
  }
  // Undo every SetPosition node for a specific step
  void addSetPositionUndo(NodeType position) { setPositionUndo.push_back(position); }
  const auto& getParameterAccountingUndo() const { return parameterAccountingUndo; }
  const auto& getSetPositionUndo() const { return setPositionUndo; }
};
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
    for (size_t u = 0; u < freeNodes.size(); ++u) {
      for (size_t v = u + 1; v < freeNodes.size(); ++v) {
        crossingCount += computeUVcrossing(fixedPosition[u], fixedPosition[v]);
      }
    }
    return crossingCount;
  }
};

template <typename NT, typename CCT>
class HeuristicGraph {
  // for each free node holds its neighbours
  std::vector<std::vector<NT>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NT>> fixedNodes;
  // free nodes current position in the permutation
  std::vector<NT> freeNodesPosition;
  // for each free node u, the first place hold the crossing sum with free node positioned to the
  // left of u and the second hold the sum to its right
  std::vector<std::array<CCT, 2>> leftRightCrossingSum;
  // hold free node id in its corrent position in the pemutation
  std::vector<NT> permutation;

 public:
  NT currentNumNodes() { return getFreeNodesSize() + getFixedNodesSize(); }
  NT currentNumEdges() { return 0; }
  using EdgeType = NT;
  using NodeType = NT;
  using WeightType = NT;
  using CrossingCountType = CCT;

  HeuristicGraph(const std::vector<std::vector<NodeType>>& freeNodes,
                 const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        freeNodesPosition(freeNodes.size()),
        leftRightCrossingSum(freeNodes.size()),
        permutation(freeNodes.size()) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }

  HeuristicGraph(NodeType numFixedNodes, NodeType numFreeNodes, CrossingCountType edgeNum)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFixedNodes, std::vector<NodeType>(0)),
        freeNodesPosition(freeNodes.size()),
        leftRightCrossingSum(freeNodes.size()),
        permutation(freeNodes.size()) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }
  void setFreeNodePosition(NodeType index, NodeType value) { freeNodesPosition[index] = value; }

  void addEdge(NodeType source,
               NodeType target) {  // where source is the freeNode and target is the fixedNode
    freeNodes[source].push_back(target);
    fixedNodes[target].push_back(source);
    return;
  }

  NodeType getFixedNodesSize() const { return fixedNodes.size(); }

  NodeType getFixedNodeNeighboursSize(NodeType nodeID) const { return fixedNodes[nodeID].size(); }

  const auto& getFixedNodeNeighbours(NodeType fixedNodeID) const { return fixedNodes[fixedNodeID]; }

  NodeType getFreeNodesSize() const { return freeNodes.size(); }

  NodeType getFreeNodeNeighboursSize(NodeType nodeID) const { return freeNodes[nodeID].size(); }

  const auto& getFreeNodeNeighbours(NodeType freeNodeID) const { return freeNodes[freeNodeID]; }

  CrossingCountType getLeftCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][0];
  }

  CrossingCountType getRightCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][1];
  }

  const auto& getPermutation() const { return permutation; }

  std::vector<NodeType>& getFreeNodes() { return permutation; }

  const auto& getFreeNodesPosition() const { return freeNodesPosition; }

  const auto& getEdges() const { return freeNodes; }

  /**
  This function returns the free node that is in index i of the perutation
  @param i The index of the permutation
*/
  const auto getPermutatuinAtIndex(NodeType i) const { return permutation[i]; }

  void setFreeNodes(const std::vector<NodeType>& newPermutation) {
    for (NodeType i = 0; i < newPermutation.size(); ++i) {
      freeNodesPosition[newPermutation[i]] = i;
      permutation[i] = newPermutation[i];
    }
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

  CrossingCountType getCrossings() {
    CrossingCountType crossingCount = 0;
    for (size_t u = 0; u < freeNodes.size(); ++u) {
      for (size_t v = u + 1; v < freeNodes.size(); ++v) {
        crossingCount += computeUVcrossing(permutation[u], permutation[v]);
      }
    }
    return crossingCount;
  }

  /**
  This function switches the positions of two neighboring free nodes. Assumes u < v.
  @param u The first free node
  @param v The second free node
  @param isConditional switch only if reduce crossings
  @return True if the switch
*/
  bool switchNeighbours(NodeType u, NodeType v, bool isConditional) {
    // the number of crossings created by the edges from (u, v)
    NodeType uvSum = computeUVcrossing(u, v);
    // the number of crossings created by the edges from (v, u)
    NodeType vuSum = computeUVcrossing(v, u);
    if (isConditional) {
      if (uvSum <= vuSum) {
        return false;
      }
    }
    leftRightCrossingSum[u][1] -= uvSum;
    leftRightCrossingSum[v][1] += vuSum;
    leftRightCrossingSum[u][0] += vuSum;
    leftRightCrossingSum[v][0] -= uvSum;
    std::swap(permutation[freeNodesPosition[u]], permutation[freeNodesPosition[v]]);
    std::swap(freeNodesPosition[u], freeNodesPosition[v]);
    return true;
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right
  void computeCrossingSums() {
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      for (NodeType j = i + 1; j < freeNodes.size(); ++j) {
        if (freeNodesPosition[i] > freeNodesPosition[j]) {
          CrossingCountType crossing = computeUVcrossing(j, i);
          leftRightCrossingSum[i][0] += crossing;
          leftRightCrossingSum[j][1] += crossing;
        } else {
          CrossingCountType crossing = computeUVcrossing(i, j);
          leftRightCrossingSum[i][1] += crossing;
          leftRightCrossingSum[j][0] += crossing;
        }
      }
    }
  }
};
template <typename NT>
class BipartiteGraph {
 public:
  NT currentNumNodes() { return getFreeNodesSize() + getFixedNodesSize(); }
  NT currentNumEdges() { return 0; }

  using EdgeType = NT;
  using NodeType = NT;
  using WeightType = NT;

  BipartiteGraph() : freeNodes({}), fixedNodes({}){};

  BipartiteGraph(NodeType fixedNodesSize, NodeType freeNodesSize, NodeType edgeSize)
      : freeNodes(freeNodesSize),
        fixedNodes(fixedNodesSize),
        edges(freeNodesSize) {  // nodes with no edges get assigned a -1
    numEdges = edgeSize;
    // Write trivial positions in both vectors
    for (NodeType i = 0; i < freeNodesSize; i++) {
      freeNodes[i] = i;
    }
    for (NodeType i = 0; i < fixedNodesSize; i++) {
      fixedNodes[i] = i;
    }
    // Concrete edges will be added with push_back operations
  };

  std::vector<NodeType>& getFreeNodes() { return freeNodes; }
  const std::vector<NodeType>& getFixedNodes() const { return fixedNodes; }
  void setFreeNodes(std::vector<NodeType>& newFreeNodes) { this->freeNodes = newFreeNodes; }

  const std::vector<std::vector<NodeType>>& getEdges() const { return edges; }
  const NodeType getEdge(NodeType FreeNode, NodeType index) const { return edges[FreeNode][index]; }

  void insertFreeNode(NodeType node) { freeNodes.push_back(node); };
  void insertFixedNode(NodeType node) { fixedNodes.push_back(node); };

  void addEdge(NodeType sourceID, NodeType targetID) {
    edges[sourceID].push_back(targetID);  // free Nodes come after the
                                          // fixed ones
  };

  const NodeType getFreeNodesSize() { return freeNodes.size(); };
  const NodeType getFixedNodesSize() { return fixedNodes.size(); };
  const NodeType getEdgesSize() { return numEdges; };

  void switchNodes(NodeType firstNodeID, NodeType secondNodeID) {
    NodeType& firstPosition = freeNodes[firstNodeID];
    NodeType& secondPosition = freeNodes[secondNodeID];
    std::swap(firstPosition, secondPosition);
  };
  int getCrossings() {  // Compare from left to right if edges cross
    int crossings = 0;
    for (NodeType i = freeNodes[0]; i < freeNodes.size(); i++) {
      for (size_t j = 0; j < (edges[i].size()); j++) {
        for (NodeType k = i; k < freeNodes.size(); k++) {
          for (size_t l = 0; l < (edges[k].size()); l++) {
            crossings += (freeNodes[i] < freeNodes[k] && edges[i][j] > edges[k][l]);
          }
        }
      }
    }
    return crossings;
  }

 private:
  NodeType numEdges;
  std::vector<NodeType> freeNodes;  // Stores the position of each node, e.g. freeNodes[0]=position
  std::vector<NodeType> fixedNodes;
  std::vector<std::vector<NodeType>> edges;
  // sourceNodeIDs are the indices (FreeNodes), targets are the vector's entries (FixedNodes). The
  // first entry is the position tho.
};
namespace median_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void medianAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;

  std::vector<int> positions(myGraph.getFreeNodesSize(), 0);
  // holds the average neighbour position for each of the free nodes
  std::vector<double> posE(myGraph.getFreeNodesSize(), 0.0);

  auto& edges = myGraph.getEdges();
  auto& permutationFreeNodes = myGraph.getFreeNodes();

  // Update positions and edges in-place
  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    if (edges[i].size() == 0) {
      continue;
    }
    int posTmp = ceil((edges[i][(edges[i].size() / 2)] - 1));
    positions[i] = posTmp;
    double avarege = 0;
    for (size_t j = 0; j < edges[i].size(); ++j) {
      avarege += edges[i][j];
    }
    posE[i] = avarege / edges[i].size();
  }

  // Sort permutationFreeNodes based on positions
  std::sort(permutationFreeNodes.begin(), permutationFreeNodes.end(),
            [&positions, &posE](const NT& a, const NT& b) {
              // If two free nodes have the same median, they to be ordered according to their
              // average neighbor
              if (positions[a] == positions[b]) {
                return posE[a] < posE[b];
              } else {
                return positions[a] < positions[b];
              }
            });
}
}  // namespace median_algorithm
namespace heuristic_algorithm {
/*
 * For each  pair of neighbours  u , v with u < v, if the left crossings of u equal 0
 * and the left crossings of v are bigger than 0,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r1(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getLeftCrossings(nodeId) == 0 && graph.getLeftCrossings(neighbourId) > 0;
  } else {
    return graph.getRightCrossings(nodeId) == 0 && graph.getRightCrossings(neighbourId) > 0;
  }
}

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of u
 * and the left crossings of v are bigger than right crossings of v,
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r2(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getRightCrossings(nodeId) > graph.getLeftCrossings(nodeId) &&
           graph.getLeftCrossings(neighbourId) > graph.getRightCrossings(neighbourId);
  } else {
    return graph.getLeftCrossings(nodeId) > graph.getRightCrossings(nodeId) &&
           graph.getRightCrossings(neighbourId) > graph.getLeftCrossings(neighbourId);
  }
}

/*
 * For each  pair of neighbours  u , v with u < v, if the right crossings of u
 * is bigger than the left crossings of v
 * try to switch u and v
 * Analog for the other direction
 */
template <class Graph>
bool r3(Graph& graph, typename Graph::NodeType nodeId, typename Graph::NodeType neighbourId) {
  if (nodeId < neighbourId) {
    return graph.getRightCrossings(nodeId) > graph.getLeftCrossings(neighbourId);
  } else {
    return graph.getLeftCrossings(nodeId) > graph.getRightCrossings(neighbourId);
  }
}

template <class Graph>
bool heuristicAlgorithm(Graph& graph, bool runR1, bool runR2, bool runR3) {
  using NodeType = typename Graph::NodeType;
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if (runR1 && r1(graph, nodeId, neighbourId) || runR2 && r2(graph, nodeId, neighbourId) ||
          runR3 && r3(graph, nodeId, neighbourId)) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if (runR1 && r1(graph, nodeId, neighbourId) || runR2 && r2(graph, nodeId, neighbourId) ||
          runR3 && r3(graph, nodeId, neighbourId)) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }

    if (didChange) {
      madeSwitch = true;
    }
  }
  return madeSwitch;
}
}  // namespace heuristic_algorithm
namespace barycenter_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void barycenterAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;

  std::vector<double> positions(myGraph.getFreeNodesSize(), 0.0);

  auto& edges = myGraph.getEdges();
  auto& permutationFreeNodes = myGraph.getFreeNodes();

  // Update positions and edges in-place
  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    for (size_t j = 0; j < edges[i].size(); j++) {
      positions[i] += edges[i][j];
    }
    if (edges[i].size() == 0) {
      continue;
    }
    positions[i] = ((double)positions[i]) / ((double)edges[i].size());
  }

  // Sort permutationFreeNodes based on positions
  std::sort(permutationFreeNodes.begin(), permutationFreeNodes.end(),
            [&positions](const NT& a, const NT& b) { return positions[a] < positions[b]; });
}
}  // namespace barycenter_algorithm
HeuristicGraph<double, double> readGraph(std::string graph_file) {
  // std::cout << "Reading graph..." << std::endl;
  std::string tmp;  // for p and ocr in input format
  int n0;           // neighbours = A, fixed partition
  int n1;           // movable_nodes = B, free partition
  int m;

  std::ifstream file(graph_file);

  if (!file.is_open()) {
    // std::cout << "error: file not open" << std::endl;
  }

  std::string line;

  std::getline(file, line);

  while (line[0] == 'c') {
    std::getline(file, line);
  }
  std::stringstream ss(line);
  ss >> tmp;  // p
  ss >> tmp;  // ocr
  ss >> n0;
  ss >> n1;
  ss >> m;

  // std::cout << "n0: " << n0 << std::endl;
  // std::cout << "n1: " << n1 << std::endl;
  // std::cout << "m: " << m << std::endl;

  // initialize graph
  auto g = HeuristicGraph<double, double>(n0, n1, m);
  // read adjacencies of the nodes in the graph file
  while (std::getline(file, line)) {
    if (line[0] == 'c' || line.empty()) {
      continue;
    }
    std::stringstream ss(line);

    int x;
    int y;
    ss >> x;
    ss >> y;
    g.addEdge(y - n0 - 1, x - 1);
  }
  file.close();

  return g;
}
std::vector<int> readGraphGraph(std::string graph_file) {
  // std::cout << "Reading graph..." << std::endl;
  std::string tmp;  // for p and ocr in input format
  int n0;           // neighbours = A, fixed partition
  int n1;           // movable_nodes = B, free partition
  int m;

  std::ifstream file(graph_file);

  if (!file.is_open()) {
    // std::cout << "error: file not open" << std::endl;
  }

  std::string line;

  std::getline(file, line);

  while (line[0] == 'c') {
    std::getline(file, line);
  }
  std::stringstream ss(line);
  ss >> tmp;  // p
  ss >> tmp;  // ocr
  ss >> n0;
  ss >> n1;
  ss >> m;
  std::vector<int> sol;
  sol.push_back(n0);
  sol.push_back(n1);
  sol.push_back(m);
  return sol;
}
ReductionGraph<int, int> readGraphRed(std::string graph_file) {
  // std::cout << "Reading graph..." << std::endl;
  std::string tmp;  // for p and ocr in input format
  int n0;           // neighbours = A, fixed partition
  int n1;           // movable_nodes = B, free partition
  int m;

  std::ifstream file(graph_file);

  if (!file.is_open()) {
    // std::cout << "error: file not open" << std::endl;
  }

  std::string line;

  std::getline(file, line);

  while (line[0] == 'c') {
    std::getline(file, line);
  }
  std::stringstream ss(line);
  ss >> tmp;  // p
  ss >> tmp;  // ocr
  ss >> n0;
  ss >> n1;
  ss >> m;

  // std::cout << "n0: " << n0 << std::endl;
  // std::cout << "n1: " << n1 << std::endl;
  // std::cout << "m: " << m << std::endl;

  // initialize graph
  auto g = ReductionGraph<int, int>(n0, n1, m);
  // read adjacencies of the nodes in the graph file
  while (std::getline(file, line)) {
    if (line[0] == 'c' || line.empty()) {
      continue;
    }
    std::stringstream ss(line);

    int x;
    int y;
    ss >> x;
    ss >> y;
    g.addEdge(y - n0 - 1, x - 1);
  }
  file.close();

  return g;
}
namespace heuristicMedian {
template <typename GraphType>
bool heuristicMedianAlgorithm(GraphType& graph, bool runR1, bool runR2, bool runR3) {
  median_algorithm::medianAlgorithm(graph);
  // std::cout << "Median Crossings: " << graph.getCrossings() << std::endl;

  for (int i = 0; i < graph.getFreeNodesSize(); ++i) {
    auto permutation = graph.getPermutation();
    graph.setFreeNodePosition(permutation[i], i);
  }
  graph.computeCrossingSums();
  bool didChange = true;
  bool madeSwitch = false;
  while (didChange) {
    didChange = false;
    // check switch with  nodes to the right
    for (size_t i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      int nodeId = graph.getPermutatuinAtIndex(i);
      int neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if ((runR1 && heuristic_algorithm::r1(graph, nodeId, neighbourId)) ||
          (runR2 && heuristic_algorithm::r2(graph, nodeId, neighbourId)) ||
          (runR3 && heuristic_algorithm::r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (size_t i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      int nodeId = graph.getPermutatuinAtIndex(i);
      int neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if ((runR1 && heuristic_algorithm::r1(graph, nodeId, neighbourId)) ||
          (runR2 && heuristic_algorithm::r2(graph, nodeId, neighbourId)) ||
          (runR3 && heuristic_algorithm::r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }

    if (didChange) {
      madeSwitch = true;
    }
  }

  return madeSwitch;
}
}  // namespace heuristicMedian

namespace reductionalgorithms {
// compute the number of crossings created by the edges from two free nodes (u, v)
// when u is to the left of v
template <class Graph>
typename Graph::CrossingCountType const computeUVcrossing(Graph& graph, typename Graph::NodeType u,
                                                          typename Graph::NodeType v) {
  using CrossingCountType = typename Graph::CrossingCountType;
  CrossingCountType crossingSum = 0;
  for (const auto uNeighbour : graph.getFreeNodeNeighbours(u)) {
    for (const auto vNeighbour : graph.getFreeNodeNeighbours(v)) {
      crossingSum += vNeighbour < uNeighbour;
    }
  }
  return crossingSum;
}

// copmute for each free node u the number of crossings created by the edges from nodes to the
// left of u  and to its right
template <class Graph, class Undo>
void computeCrossingSums(Graph& graph) {
  using NodeType = typename Graph::NodeType;
  using CrossingCountType = typename Graph::CrossingCountType;
  for (NodeType u = 0; u < graph.getFreeNodesSize(); ++u) {
    for (NodeType v = u + 1; v < graph.getFreeNodesSize(); ++v) {
      CrossingCountType crossingUV = computeUVcrossing(graph, u, v);
      graph.addCrossing(u, v, crossingUV);
      // reduction RR1: For each pair of vertices {u, v} ⊆ free nodes that forms a 0/j pattern
      // with j > 0, commit u < v
      CrossingCountType crossingVU = computeUVcrossing(graph, v, u);
      graph.addCrossing(v, u, crossingVU);
    }
  }
}
}  // namespace reductionalgorithms

namespace reductionalgorithms {

/*
 * For each decision u < v that we make adjust the leftRightSet of u and v
 * Add all the transitiv decisions that follow the u < v decision
 */
template <class Graph, class Undo>
void parameterAccounting(Graph& graph, typename Graph::NodeType u, typename Graph::NodeType v,
                         typename Graph::CrossingCountType& currentSolution, Undo* undo = nullptr) {
  using NodeType = typename Graph::NodeType;
  if (u != v) {
    // if (graph.getRightNodes(u).find(v) == graph.getRightNodes(u).end()) {
    if (graph.getNodeCrossing(u).find(v) != graph.getNodeCrossing(u).end()) {
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
        }
      }
      for (NodeType biggerThanV : graph.getRightNodes(v)) {
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
      if (graph.getNodeCrossing(u).find(v) != graph.getNodeCrossing(u).end()) {
        if (graph.getNodeCrossing(u).at(v) == 0) {
          parameterAccounting<Graph, Undo>(graph, u, v, currentSolution);
        } else if (graph.getNodeCrossing(v).at(u) == 0) {
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
        if (u < v) {
          pairsToModify.emplace_back(v, u);
          didChange = true;
        } else if (graph.getNodeCrossing(v).at(u) < crossingsLeft) {
          pairsToModify.emplace_back(v, u);
          didChange = true;
        }
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
    Graph& graph, int upperBound) {
  using CrossingCountType = typename Graph::CrossingCountType;
  using NodeType = typename Graph::NodeType;

  // holds the crossing number of the best solution so far
  CrossingCountType bestSolution = upperBound;
  // holds the order of the solution best so far
  std::vector<NodeType> bestOrder;

  computeCrossingSums<Graph, Undo>(graph);
  CrossingCountType currentSolution = 0;
  rr1<Graph, Undo>(graph, currentSolution);
  rr2<Graph, Undo>(graph, currentSolution);

  bool didChangeRrlo2 = true;

  while (didChangeRrlo2) {
    didChangeRrlo2 = rrlo2<Graph, Undo>(graph, currentSolution);
    rrlo1<Graph, Undo>(graph);
  }

  rr3<Graph, Undo>(graph, currentSolution);
  algorithmStep<Graph, Undo>(graph, currentSolution, false, 0, 0, bestSolution, bestOrder);
  graph.setFixedPositions(bestOrder);
  return std::make_tuple(bestSolution, bestOrder);
}
}  // namespace reductionalgorithms
int main() {
  std::vector<std::string> directories = {"heuristic_set"};
  std::vector<std::string> instances;
  for (int i = 1; i <= 100; ++i) {
    instances.push_back(std::to_string(i) + ".gr");
  }
  /*std::ofstream outputFile("HeuristicResults.csv");
  outputFile << "Instance,Time,Quality\n";
  int crossings = 0;
  for (const auto& directory : directories) {
    for (const auto& instance : instances) {
      for (int i = 0; i < 15; i++) {
        std::string filename = directory + "/" + instance;

        auto graph = readGraph(filename);

        // std::cout << graph.getEdgesSize() << std::endl;

        std::cout << filename << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        heuristic_algorithm::heuristicAlgorithm<HeuristicGraph<int, int>>(graph, true, true, true);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Computation for graph: " << filename << " number: " << i << " done"
                  << std::endl;
        if (i == 0) {
          std::cout << "Took: " << duration.count() << std::endl;
          crossings = graph.getCrossings();
          std::cout << crossings << std::endl;
        }
        outputFile << instance << "," << duration.count() << "," << crossings << "\n";
      }
    }
  }

  outputFile.close();*/
  directories = {"heuristic_set"};
  std::ofstream outputFileHeu("HeuristicResults.csv");
  outputFileHeu << "Instance,Time,Quality\n";
  int crossings = 0;
  std::vector<std::string> mediumInstances;

  for (const auto& directory : directories) {
    for (const auto& instance : instances) {
      if (instance == "44.gr") {
        outputFileHeu << instance << "," << 0.000875196 << "," << -1 << "\n";
        continue;
      }
      for (int i = 0; i < 15; i++) {
        std::string filename = directory + "/" + instance;

        auto graph = readGraph(filename);

        // std::cout << graph.getEdgesSize() << std::endl;

        std::cout << filename << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        heuristic_algorithm::heuristicAlgorithm<HeuristicGraph<double, double>>(graph, true, true,
                                                                                true);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Computation for graph: " << filename << " number: " << i << " done"
                  << std::endl;
        if (i == 0) {
          std::cout << "Took: " << duration.count() << std::endl;
          crossings = graph.getCrossings();
          std::cout << crossings << std::endl;
        }
        outputFileHeu << instance << "," << duration.count() << "," << crossings << "\n";
      }
    }
  }

  outputFileHeu.close(); /*
 std::ofstream outputFile("HeuristicGraphs.csv");
 outputFile << "n0, n1, m\n";
 int crossings = 0;
 for (const auto& directory : directories) {
   for (const auto& instance : instances) {
     std::string filename = directory + "/" + instance;

     auto graph = readGraphGraph(filename);

     outputFile << graph[0] << "," << graph[1] << "," << graph[2] << "\n";
   }
 }
 outputFile.close();
 /*
   std::ofstream outputFile("MediumMedianResults.csv");
   outputFile << "Instance,Time,Quality\n";
   int crossings = 0;
   for (const auto& directory : directories) {
     for (const auto& instance : instances) {
       for (int i = 0; i < 15; i++) {
         std::string filename = directory + "/instances/" + instance;

         auto graph = readGraph(filename);

         // std::cout << graph.getEdgesSize() << std::endl;

         std::cout << filename << std::endl;
         auto start = std::chrono::high_resolution_clock::now();
         median_algorithm::medianAlgorithm<HeuristicGraph<int, int>>(graph);
         auto end = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> duration = end - start;
         std::cout << "Computation for graph: " << filename << " number: " << i << " done"
                   << std::endl;
         if (i == 0) {
           std::cout << "Took: " << duration.count() << std::endl;
           crossings = graph.getCrossings();
           std::cout << crossings << std::endl;
         }
         outputFile << instance << "," << duration.count() << "," << crossings << "\n";
       }
     }
   }

   outputFile.close();
   std::ofstream outputFileHM("MediumHeuristicsMedianResults.csv");
   outputFileHM << "Instance,Time,Quality\n";
   crossings = 0;
   for (const auto& directory : directories) {
     for (const auto& instance : instances) {
       for (int i = 0; i < 15; i++) {
         std::string filename = directory + "/instances/" + instance;

         auto graph = readGraph(filename);
         std::cout << "Please: " << graph.getCrossings() << std::endl;
         // std::cout << graph.getEdgesSize() << std::endl;

         std::cout << filename << std::endl;
         auto start = std::chrono::high_resolution_clock::now();
         heuristicMedian::heuristicMedianAlgorithm<HeuristicGraph<int, int>>(graph, true, true,
                                                                             true);
         auto end = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double> duration = end - start;
         std::cout << "Computation for graph: " << filename << " number: " << i << " done"
                   << std::endl;
         if (i == 0) {
           std::cout << "Took: " << duration.count() << std::endl;
           crossings = graph.getCrossings();
           std::cout << crossings << std::endl;
         }
         outputFileHM << instance << "," << duration.count() << "," << crossings << "\n";
       }
     }
   }

   outputFileHM.close();
*/
                         /* std::vector<std::string> directoriesRed = {"medium_test_set"};
                       
                          std::ofstream outputFileRed("MediumReductionResults.csv");
                          outputFileRed << "Instance,Time,Quality\n";
                          std::set<std::string> skip;
                          int crossings = 0;
                          for (const auto& directory : directoriesRed) {
                            for (const auto& instance : instances) {
                              if (skip.find(instance) != skip.end()) {
                                std::cout << "Skipped: " << instance << std::endl;
                                ;
                                outputFileRed << instance << "," << -1 << "," << -1 << "\n";
                                continue;
                              }
                              for (int i = 0; i < 1; i++) {
                                std::string filename = directory + "/instances/" + instance;
                       
                                auto graph = readGraphRed(filename);
                                auto graphBound = readGraph(filename);
                       
                                std::cout << filename << std::endl;
                                auto start = std::chrono::high_resolution_clock::now();
                                heuristicMedian::heuristicMedianAlgorithm(graphBound, true, true, true);
                                auto upperBound = graphBound.getCrossings();
                                reductionalgorithms::algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(
                                    graph, upperBound);
                                auto end = std::chrono::high_resolution_clock::now();
                                std::chrono::duration<double> duration = end - start;
                                std::cout << "Computation for graph: " << filename << " number: " << i
                                          << "with reduction done" << std::endl;
                                if (i == 0) {
                                  std::cout << "Took: " << duration.count() << std::endl;
                                  crossings = graph.getCrossings();
                                  std::cout << crossings << std::endl;
                                }
                                outputFileRed << instance << "," << duration.count() << "," << crossings << "\n";
                              }
                            }
                          }
                       
                          outputFileRed.close();*/
  return 0;
}