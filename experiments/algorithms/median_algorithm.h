#pragma once

#include "ds/bipartite_graph.h"
#include "ds/helper/empty_problem.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/median_algorithm/median_algorithm.h"
#include "toolkit/algorithms/algorithm_impl.h"
#include "toolkit/app/app_io.pb.h"
#include "toolkit/ds/empty_problem.h"

namespace oscm::experiments::algorithms {
struct MedianGraphEmptyProblem : public oscm::ds::EmptyProblem<BipartiteGraph<int>> {
  const static constexpr std::string_view ds_name = "bipartite_graph";
};
namespace {
using henrixapp::app::app_io::AlgorithmConfig;
using henrixapp::app::app_io::AlgorithmRunInformation;
using henrixapp::app::app_io::Hypergraph;
using henrixapp::app::app_io::Result;
using henrixapp::app::app_io::RunConfig;
using BG = MedianGraphEmptyProblem;
}  // namespace
class MedianAlgorithm : public henrixapp::algorithms::AlgorithmImpl<BG> {
 public:
  static constexpr absl::string_view AlgorithmName = "median_algorithm";

 protected:
  absl::StatusOr<std::unique_ptr<BG>> Execute(const AlgorithmConfig& config,
                                              std::unique_ptr<BG> problem) override {
    median_algorithm::medianAlgorithm<BipartiteGraph<int>>(problem->instance());

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