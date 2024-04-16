#pragma once

#include "ds/bipartite_graph.h"
#include "ds/helper/empty_problem.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/barycenter_algorithm/barycenter_algorithm.h"
#include "toolkit/algorithms/algorithm_impl.h"
#include "toolkit/app/app_io.pb.h"
#include "toolkit/ds/empty_problem.h"
namespace oscm::experiments::algorithms {
struct BipartiteGraphProblem : public oscm::ds::EmptyProblem<BipartiteGraph<int>> {
  const static constexpr std::string_view ds_name = "bipartite_graph";
};
}  // namespace oscm::experiments::algorithms