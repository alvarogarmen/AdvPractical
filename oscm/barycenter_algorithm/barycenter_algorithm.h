#pragma once
#include <algorithm>
#include <iostream>
#include <vector>
namespace barycenter_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void barycenterAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;
  std::vector<std::pair<NT, doubles>> meanPositions(myGraph.getFreeNodesSize());
  for (size_t i = 0; i < meanPositions.size(); ++i) {
    meanPositions[i] = std::make_pair(0, i);
  }

  auto& edges = myGraph.getEdges();

  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    for (size_t j = 0; j < edges[i].size(); j++) {
      meanPositions[i].first += edges[i][j];
    }
    meanPositions[i].first = (double)(meanPositions[i].first) / (double)(edges[i].size());
  }

  std::sort(meanPositions.begin(), meanPositions.end(),
            [](const std::pair<NT, int>& a, const std::pair<NT, int>& b) {
              // Compare the first elements (mean positions) for sorting
              return a.first < b.first;
            });

  // TODO: DO it in place
  std::vector<NT> newFreeNodes(meanPositions.size());

  for (size_t i = 0; i < edges.size(); i++) {
    newFreeNodes[i] = meanPositions[i].second;  // Move the nodes to their new positions
  }

  myGraph.setFreeNodes(newFreeNodes);
}
}  // namespace barycenter_algorithm
