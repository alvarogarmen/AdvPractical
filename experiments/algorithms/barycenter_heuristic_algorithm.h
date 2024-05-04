#pragma once

#include "ds/helper/empty_problem.h"
#include "ds/heuristic_graph/heuristic_graph.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/barycenter_heuristic_algorithm/barycenter_heuristic_algorithm.h"
#include "problems.h"
#include "toolkit/algorithms/algorithm_impl.h"
#include "toolkit/app/app_io.pb.h"
#include "toolkit/ds/empty_problem.h"

namespace oscm::experiments::algorithms {

namespace {
using henrixapp::app::app_io::AlgorithmConfig;
using henrixapp::app::app_io::AlgorithmRunInformation;
using henrixapp::app::app_io::Hypergraph;
using henrixapp::app::app_io::Result;
using henrixapp::app::app_io::RunConfig;
using BG = HeuristicGraphEmptyProblem;
}  // namespace
class BarycenterHeuristicAlgorithm : public henrixapp::algorithms::AlgorithmImpl<BG> {
 public:
  static constexpr absl::string_view AlgorithmName = "barycenter_heuristic_algorithm";

 protected:
  absl::StatusOr<std::unique_ptr<BG>> Execute(const AlgorithmConfig& config,
                                              std::unique_ptr<BG> problem) override {
    barycenterHeuristic_algorithm::barycenterHeuristicAlgorithm<HeuristicGraph<int, int>>(
        problem->instance());

    return problem;
  }
  absl::Status ValidateConfig(const AlgorithmConfig& config) override { return absl::OkStatus(); }
  absl::StatusOr<std::unique_ptr<BG>> Load(const RunConfig& run_config,
                                           const Hypergraph& hypergraph) override {
    // auto hgr = henrixapp::algorithms::loadFile<StandardIntegerHypergraph>(hypergraph);
    if (hypergraph.format() == "gr") {
      std::ifstream file(hypergraph.file_path());
      auto Graph = readGraph<HeuristicGraph<int, int>>(file);
      if (!Graph.ok()) {
        return Graph.status();
      }
      return std::make_unique<BG>(std::move(Graph.value()));
    }
    return absl::UnimplementedError("Read-in-function for this format not implemented");
  }
};
}  // namespace oscm::experiments::algorithms