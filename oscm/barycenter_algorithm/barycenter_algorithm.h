#pragma once
#include <algorithm>
#include <iostream>
#include <vector>
namespace barycenter_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void barycenterAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;
  std::vector<std::pair<NT, int>> meanPositions(myGraph.getFreeNodesSize());
  for (size_t i = 0; i < meanPositions.size(); ++i) {
    meanPositions[i] = std::make_pair(0, i);
  }

  auto& edges = myGraph.getEdges();

  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    if (edges[i].size() == 0) {
      meanPositions[i].first = -1;
      continue;
    }
    for (size_t j = 0; j < edges[i].size(); j++) {
      meanPositions[i].first += edges[i][j];
    }
    meanPositions[i].first = meanPositions[i].first / edges[i].size();
  }

  std::sort(meanPositions.begin(), meanPositions.end(),
            [](const std::pair<NT, int>& a, const std::pair<NT, int>& b) {
              // Compare the first elements (mean positions) for sorting
              return a.first < b.first;
            });
  std::vector<NT> newFreeNodes(meanPositions.size());

  for (size_t i = 0; i < edges.size(); i++) {
    newFreeNodes[i] = meanPositions[i].second;  // Move the nodes to their new positions
  }

  myGraph.setFreeNodes(newFreeNodes);
}

// Here we also copy the edges over, so that it can work as a starting point for another algorithm
template <typename GraphType>
void barycenterAlgorithmEdges(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;
  std::vector<std::pair<NT, int>> meanPositions(myGraph.getFreeNodesSize());
  for (size_t i = 0; i < meanPositions.size(); ++i) {
    meanPositions[i] = std::make_pair(0, i);
  }

  auto& edges = myGraph.getEdges();

  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    for (size_t j = 0; j < edges[i].size(); j++) {
      if (edges[i].size() == 0) {
        meanPositions[i].first = -1;
        continue;
      }
      meanPositions[i].first += edges[i][j];
    }
    meanPositions[i].first = meanPositions[i].first / edges[i].size();
  }

  std::sort(meanPositions.begin(), meanPositions.end(),
            [](const std::pair<NT, int>& a, const std::pair<NT, int>& b) {
              // Compare the first elements (mean positions) for sorting
              return a.first < b.first;
            });
  std::vector<NT> newFreeNodes(meanPositions.size());

  for (size_t i = 0; i < edges.size(); i++) {
    newFreeNodes[i] = meanPositions[i].second;  // Move the nodes to their new positions
  }

  myGraph.setFreeNodes(newFreeNodes);
  myGraph.setEdges(newFreeNodes);
}
}  // namespace barycenter_algorithm