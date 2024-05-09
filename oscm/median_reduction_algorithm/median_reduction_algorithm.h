#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

#include "oscm/median_algorithm/median_algorithm.h"
#include "oscm/reduction_algorithm/reduction_algorithm.h"
namespace medianReduction_algorithm {

template <typename GraphType, class Undo>
std::tuple<typename Graph::CrossingCountType, std::vector<typename Graph::NodeType>>
medianReductionAlgorithm(GraphType& myGraph) {
  median_algorithm::medianAlgorithm(myGraph);
  auto [crossingSum, orderVector] =
      reductionalgorithms::reductionAlgorithm(myGraph, true, true, true);
  return std::make_tuple(crossingSum, orderVector);
}
}  // namespace medianReduction_algorithm