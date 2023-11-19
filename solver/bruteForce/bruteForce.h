#pragma once

#include <algorithm>
#include <limits>
#include <vector>

template <typename Edge, typename SizeType>
std::vector<SizeType> findBestPermutation(Graph<SizeType>& myGraph) {
  std::vector<SizeType> bestPermutation = myGraph.freeNodes;
  int minCrossings = std::numeric_limits<int>::max();

  // Generate permutations of edges

  auto crossings = crossGrader(myGraph);
  if (crossings < minCrossings) {
    minCrossings = crossings;
    bestPermutation = myGraph.edges;
  }
  // Here update the currentpermutation to the next one
  // With the current interface it looks kinda dubious
  return bestPermutation;
};