template <class Graph>
class EmptyProblem {
 public:
  std::unique_ptr<Graph> _g;
  Graph& instance() { return *_g; }
  using Graph_t = Graph;
  using NodeID_t = typename Graph::NodeType;
  using EdgeID_t = typename Graph::EdgeType;
  using NodeType = typename Graph::NodeType;

  // We don´t use these types. Not sure if I can delete them without braking anything
  // using EdgeType = typename Graph::EdgeType;
  // using WeightType = typename Graph::WeightType;

  EmptyProblem(std::unique_ptr<Graph> g) : _g(std::move(g)) {}

  bool valid() { return true; }

  size_t size() const { return 0; }

  long weight() const { return 0; }
  long free_edges_size() { return 0; }
  double quality() const { return 0; }

  const static constexpr std::string_view ds_name = "empty_problem";
};