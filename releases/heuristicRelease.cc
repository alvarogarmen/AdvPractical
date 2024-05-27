#include <stdio.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
#include <vector>

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

  const auto& getFreeNodesPosition() const { return freeNodesPosition; }

  const auto& getEdges() const { return freeNodes; }

  std::vector<NT>& getFreeNodes() { return permutation; }

  void setFreeNodePosition(NodeType index, NodeType value) { freeNodesPosition[index] = value; }

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

  void switchNodes(NodeType u, NodeType v, int& bestSolution) {
    std::swap(permutation[freeNodesPosition[u]], permutation[freeNodesPosition[v]]);
    std::swap(freeNodesPosition[u], freeNodesPosition[v]);
    int currentSolution = this->getCrossings();
    if (currentSolution < bestSolution) {
      bestSolution = currentSolution;
      return;
    }
    bestSolution = currentSolution;
    // If not improved, swap back
    std::swap(permutation[freeNodesPosition[u]], permutation[freeNodesPosition[v]]);
    std::swap(freeNodesPosition[u], freeNodesPosition[v]);
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right
  void computeCrossingSums() {
    for (size_t i = 0; i < freeNodes.size(); ++i) {
      for (size_t j = i + 1; j < freeNodes.size(); ++j) {
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
template <typename BipartiteGraph>
std::unique_ptr<BipartiteGraph> readGraph(std::string pathToGraph) {
  std::string trash;
  typename BipartiteGraph::NodeType n0;
  typename BipartiteGraph::NodeType n1;
  typename BipartiteGraph::NodeType m;
  std::ifstream file(pathToGraph);

  std::string line;

  std::getline(file, line);

  while (line[0] == 'c') {
    std::getline(file, line);
  }

  std::stringstream ss(line);

  ss >> trash;

  ss >> trash;

  ss >> n0;

  ss >> n1;

  ss >> m;

  std::unique_ptr bipartiteGraph = std::make_unique<BipartiteGraph>(n0, n1, m);

  typename BipartiteGraph::NodeType source, target;

  while (std::getline(file, line)) {
    if (line[0] == 'c' || line.empty()) {
      continue;
    }

    std::stringstream ss(line);

    ss >> target;

    ss >> source;

    bipartiteGraph->addEdge(source - n0 - 1, target - 1);
  }
  file.close();
  return bipartiteGraph;
}

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
bool heuristicAlgorithm(Graph& graph, bool runR1, bool runR2, bool runR3,
                        std::atomic<bool>& terminationRequested) {
  using NodeType = typename Graph::NodeType;

  median_algorithm::medianAlgorithm(graph);

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
    for (NodeType i = 0; i < graph.getFreeNodesSize() - 1; ++i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i + 1);
      if ((runR1 && r1(graph, nodeId, neighbourId)) || (runR2 && r2(graph, nodeId, neighbourId)) ||
          (runR3 && r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(nodeId, neighbourId, true)) {
          didChange = true;
        }
      }
    }
    // check switch with  nodes to the left
    for (NodeType i = graph.getFreeNodesSize() - 1; i > 0; --i) {
      NodeType nodeId = graph.getPermutatuinAtIndex(i);
      NodeType neighbourId = graph.getPermutatuinAtIndex(i - 1);
      if ((runR1 && r1(graph, nodeId, neighbourId)) || (runR2 && r2(graph, nodeId, neighbourId)) ||
          (runR3 && r3(graph, nodeId, neighbourId))) {
        if (graph.switchNeighbours(neighbourId, nodeId, true)) {
          didChange = true;
        }
      }
    }

    if (didChange) {
      madeSwitch = true;
    }
  }

  int bestSolution = graph.getCrossings();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, graph.getFreeNodesSize() - 1);
  while (!terminationRequested.load()) {
    NodeType node1 = dis(gen);
    NodeType node2 = dis(gen);

    if (node1 != node2) {
      graph.switchNodes(node1, node2, bestSolution);
    }
  }

  return madeSwitch;
}
}  // namespace heuristic_algorithm
std::atomic<bool> terminationRequested(false);

void signalHandler(int signum) {
  if (signum == SIGTERM) {
    terminationRequested.store(true);
  }
}

int main(int argc, char* argv[]) {
  std::signal(SIGTERM, signalHandler);
  std::ifstream inputFile(argv[1]);
  std::unique_ptr graph = readGraph<HeuristicGraph<int, int>>(argv[1]);
  // auto graph = std::move(result.move());
  auto vector = graph->getPermutation();

  HeuristicGraph<int, int>& myGraph = *graph;
  // Create an output file stream
  heuristic_algorithm::heuristicAlgorithm(myGraph, true, true, true, terminationRequested);
  auto solution = myGraph.getPermutation();
  // Check if the file is open
  std::ofstream outputFile(argv[2]);
  int n0 = myGraph.getFixedNodesSize();
  for (size_t i = 0; i < solution.size(); i++) {
    outputFile << solution[i] + 1 + n0 << std::endl;
  }
  // Close the file stream
  outputFile.close();

  return 0;
}