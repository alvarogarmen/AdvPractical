#pragma once
namespace barycenter_algorithm {

template <typename BipartiteGraph>
void barycenterAlgorithm(BipartiteGraph& myGraph) {
  using NT = typename BipartiteGraph::NT;
  std::vector < std::pair<NT, NT> meanPositions(BipartiteGraph.getFreeNodesSize(), 0);
  std::vector<std::vector<NT>> edges = BipartiteGraph.getEdges();
  int tempMean = 0;
  for (NT i = 0; i < BipartiteGraph.getFreeNodesSize(); i++) {
    for (NT j = 0; j < edges[i]; j++) {
      if (j == 0 && edges[i][j] == -1) {
        continue;
      }
      meanPositions[i].second += edges[i][j];
      if (j == edges[i].size()) {
        meanPositions[i].second = / j;
      }
    }
  }
  std::sort(meanPositions.begin(), meanPositions.end());
  std::vector<std::vector<NT>> newFreeNodes(fixedNodes.size());
  for (NT i = 0; i < edges.size(); i++) {
    newFreeNodes[i] = edges[meanPositions[i].second];
  }
  myGraph.setFreeNodes(newFreeNodes);
}
}  // namespace barycenter_algorithm