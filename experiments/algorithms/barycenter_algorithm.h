#pragma once

#include "ds/bipartite_graph.h"
#include "ds/helper/empty_problem.h"
#include "experiments/algorithms/problems.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/barycenter_algorithm/barycenter_algorithm.h"
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
}  // namespace
class BarycenterAlgorithm : public henrixapp::algorithms::AlgorithmImpl<BipartiteGraphProblem> {
 public:
  using BG = BipartiteGraphProblem;
  static constexpr absl::string_view AlgorithmName = "barycenter_algorithm";

 protected:
  absl::StatusOr<std::unique_ptr<BG>> Execute(const AlgorithmConfig& config,
                                              std::unique_ptr<BG> problem) override {
    barycenter_algorithm::barycenterAlgorithm<BipartiteGraph<int>>(problem->instance());

    return problem;
  }
  absl::Status ValidateConfig(const AlgorithmConfig& config) override { return absl::OkStatus(); }
  absl::StatusOr<std::unique_ptr<BG>> Load(const RunConfig& run_config,
                                           const Hypergraph& hypergraph) override {
    // auto hgr = henrixapp::algorithms::loadFile<StandardIntegerHypergraph>(hypergraph);
    if (hypergraph.format() == "gr") {
      std::ifstream file(hypergraph.file_path());
      auto Graph = readGraph<BipartiteGraph<int>>(file);
      if (!Graph.ok()) {
        return Graph.status();
      }
      return std::make_unique<BG>(std::move(Graph.value()));
    }
    return absl::UnimplementedError("Read-in-function for this format not implemented");
  }
};
}  // namespace oscm::experiments::algorithms