#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

namespace median_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void medianAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;

  std::vector<int> positions(myGraph.getFreeNodesSize(), 0);
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
    for (NT j = 0; j < edges[i].size(); ++j) {
      avarege += edges[i][j];
    }
    posE[i] = avarege / edges[i].size();

    if (edges[i].size() % 2 == 0) {
      // positions[i] += ;
    }
  }

  // Sort permutationFreeNodes based on positions
  std::sort(permutationFreeNodes.begin(), permutationFreeNodes.end(),
            [&positions, &posE](const NT& a, const NT& b) {
              if (positions[a] == positions[b]) {
                return posE[a] < posE[b];
              } else {
                return positions[a] < positions[b];
              }
            });
}
}  // namespace median_algorithm
   // namespace median_algorithm