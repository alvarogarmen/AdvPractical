#pragma once

#include "ds/helper/empty_problem.h"
#include "ds/heuristic_graph/heuristic_graph.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/heuristic_algorithm/heuristic_algorithm.h"
#include "toolkit/algorithms/algorithm_impl.h"
#include "toolkit/app/app_io.pb.h"
#include "toolkit/ds/empty_problem.h"

namespace oscm::experiments::algorithms {
struct HeuristicGraphEmptyProblem : public oscm::ds::EmptyProblem<HeuristicGraph<int, int>> {
  const static constexpr std::string_view ds_name = "heuristic_graph";
};
namespace {
using henrixapp::app::app_io::AlgorithmConfig;
using henrixapp::app::app_io::AlgorithmRunInformation;
using henrixapp::app::app_io::Hypergraph;
using henrixapp::app::app_io::Result;
using henrixapp::app::app_io::RunConfig;
using RG = HeuristicGraphEmptyProblem;
}  // namespace
class SimpleHeuristicAlgorithm : public henrixapp::algorithms::AlgorithmImpl<RG> {
 public:
  static constexpr absl::string_view AlgorithmName = "heuristic_algorithm";

 protected:
  absl::StatusOr<std::unique_ptr<RG>> Execute(const AlgorithmConfig& config,
                                              std::unique_ptr<RG> problem) override {
    heuristic_algorithm::algorithm<HeuristicGraph<int, int>>(problem->instance(), true, true, true);

    return problem;
  }
  absl::Status ValidateConfig(const AlgorithmConfig& config) override { return absl::OkStatus(); }
  absl::StatusOr<std::unique_ptr<RG>> Load(const RunConfig& run_config,
                                           const Hypergraph& hypergraph) override {
    // auto hgr = henrixapp::algorithms::loadFile<StandardIntegerHypergraph>(hypergraph);
    if (hypergraph.format() == "gr") {
      std::ifstream file(hypergraph.file_path());
      auto Graph = readGraph<HeuristicGraph<int, int>>(file);
      if (!Graph.ok()) {
        return Graph.status();
      }
      return std::make_unique<RG>(std::move(Graph.value()));
    }
    return absl::UnimplementedError("Read-in-function for this format not implemented");
  }
};
}  // namespace oscm::experiments::algorithms