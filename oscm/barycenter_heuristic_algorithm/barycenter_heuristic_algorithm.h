#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

#include "oscm/barycenter_algorithm/barycenter_algorithm.h"
#include "oscm/heuristic_algorithm/heuristic_algorithm.h"
namespace barycenterHeuristic_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void barycenterHeuristicAlgorithm(GraphType& myGraph) {
  barycenter_algorithm::barycenterAlgorithm(myGraph);
  heuristic_algorithm::heuristicAlgorithm(myGraph, true, true, true);
}
}  // namespace barycenterHeuristic_algorithm