
#pragma once
#include <memory>

#include "absl/strings/string_view.h"

namespace oscm {
namespace ds {

template <class Graph>
class EmptyProblem {
 public:
  std::unique_ptr<Graph> _g;
  Graph& instance() { return *_g; }
  using Graph_t = Graph;
  using NodeID_t = typename Graph::NodeType;
  using EdgeID_t = typename Graph::EdgeType;
  using NodeType = typename Graph::NodeType;
  using EdgeType = typename Graph::EdgeType;
  using WeightType = typename Graph::WeightType;

  EmptyProblem(std::unique_ptr<Graph> g) : _g(std::move(g)) {}

  size_t crossings = 0;

  bool valid() { return true; }

  size_t size() const { return _g->getCrossings(); }

  int weight() const { return _g->getCrossings(); }
  int free_edges_size() { return 0; }
  double quality() const { return 0; }

  void save(std::string filename) {}
  const static constexpr std::string_view ds_name = "empty_problem";
};
}  // namespace ds
}  // namespace oscm