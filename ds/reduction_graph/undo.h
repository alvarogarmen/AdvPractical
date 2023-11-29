#pragma once
#include <vector>

template <typename NodeType, typename CrossingCountType>
struct Operation {
  NodeType leftNode;
  NodeType rightNode;
  CrossingCountType leftRightCrossing;
  CrossingCountType rightLeftCrossing;

  Operation(NodeType leftNode, NodeType rightNode, CrossingCountType leftRightCrossing,
            CrossingCountType rightLeftCrossing)
      : leftNode(leftNode),
        rightNode(rightNode),
        leftRightCrossing(leftRightCrossing),
        rightLeftCrossing(rightLeftCrossing) {}
};

template <typename NodeType, typename CrossingCountType>
class Undo {
  std::vector<Operation<NodeType, CrossingCountType>> parameterAccountingUndo;
  std::vector<NodeType> setPositionUndo;

 public:
  Undo() {}

  void addParameterAccountingUndo(NodeType leftNode, NodeType rightNode,
                                  CrossingCountType leftRightCrossing,
                                  CrossingCountType rightLeftCrossing) {
    Operation opUndo = Operation(leftNode, rightNode, leftRightCrossing, rightLeftCrossing);
    parameterAccountingUndo.push_back(opUndo);
  }

  void addSetPositionUndo(NodeType position) { setPositionUndo.push_back(position); }
  const auto& getParameterAccountingUndo() const { return parameterAccountingUndo; }
  const auto& getSetPositionUndo() const { return setPositionUndo; }
};