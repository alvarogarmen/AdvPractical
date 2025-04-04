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

template <typename NT, typename CCT>
class UndoAlgorithmStep {
  std::vector<Operation<NT, CCT>> parameterAccountingUndo;
  std::vector<NT> setPositionUndo;

 public:
  using NodeType = NT;
  using CrossingCountType = CCT;

  UndoAlgorithmStep() {}
  // Undo every ParameterAccounting for a specific step
  void addParameterAccountingUndo(NodeType leftNode, NodeType rightNode,
                                  CrossingCountType leftRightCrossing,
                                  CrossingCountType rightLeftCrossing) {
    parameterAccountingUndo.push_back(
        Operation(leftNode, rightNode, leftRightCrossing, rightLeftCrossing));
  }
  // Undo every SetPosition node for a specific step
  void addSetPositionUndo(NodeType position) { setPositionUndo.push_back(position); }
  const auto& getParameterAccountingUndo() const { return parameterAccountingUndo; }
  const auto& getSetPositionUndo() const { return setPositionUndo; }
};