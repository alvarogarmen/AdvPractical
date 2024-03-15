#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

namespace barycenter_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void barycenterAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;

  std::vector<double> positions(myGraph.getFreeNodesSize(), 0.0);

  auto& edges = myGraph.getEdges();
  auto& freeNodes = myGraph.getFreeNodes();

  // Update positions and edges in-place
  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    for (size_t j = 0; j < edges[i].size(); j++) {
      positions[i] += edges[i][j];
    }
    if (edges[i].size() == 0) {
      continue;
    }
    positions[i] = (double)(positions[i]) / (double)edges[i].size();
  }

  // Sort freeNodes based on positions
  std::sort(freeNodes.begin(), freeNodes.end(),
            [&positions](const NT& a, const NT& b) { return positions[a] < positions[b]; });
}
}  // namespace barycenter_algorithm