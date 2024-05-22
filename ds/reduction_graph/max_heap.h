#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

template <typename NT, typename CCT>
class MaxHeap {
  std::vector<std::tuple<NT, NT, CCT>> heap;

  static bool compare(const std::tuple<NT, NT, CCT>& a, const std::tuple<NT, NT, CCT>& b) {
    return std::get<2>(a) < std::get<2>(b);
  }

 public:
  using NodeType = NT;
  using CrossingCountType = CCT;

  MaxHeap() = default;

  inline void push(NodeType key1, NodeType key2, CrossingCountType value) {
    heap.emplace_back(key1, key2, value);
    std::push_heap(heap.begin(), heap.end(), compare);
  }

  inline void pop() {
    std::pop_heap(heap.begin(), heap.end(), compare);
    auto maxNode = heap.back();
    heap.pop_back();
    return;
  }

  inline std::tuple<NodeType, NodeType, CrossingCountType> top() const {
    if (heap.empty()) {
      throw std::runtime_error("Heap is empty");
    }
    return heap.front();
  }
  const auto getSize() { return heap.size(); }
  inline bool empty() const { return heap.empty(); }
};