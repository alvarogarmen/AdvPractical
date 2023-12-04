#pragma once
#include <fstream>
#include <iostream>
#include <string>

#include "absl/status/statusor.h"
#include "absl/strings/str_split.h"

template <typename BipartiteGraph>
absl::StatusOr<std::unique_ptr<BipartiteGraph>> readGraph(std::istream& stream,
                                                          bool verbose = false,
                                                          bool readAll = true) {
  if (!stream.good()) {
    return absl::NotFoundError("stream not good error");
  }
  std::string header;
  while (header == "") {
    if (stream.eof()) {
      return absl::NotFoundError("File does not contain enough lines. Could not read header.");
    }
    std::getline(stream, header);
  }

  std::string line;
  typename BipartiteGraph::NodeType n0 = 0;  // Number of fixed Nodes
  typename BipartiteGraph::NodeType n1 = 0;  // Number of free Nodes
  typename BipartiteGraph::NodeType m = 0;   // Number of edges

  while (std::getline(stream, line) && (line.empty() || line[0] == 'c')) {
    // Comment line, skip
  }
  // Now read in p-line
  std::vector<std::string> parts = absl::StrSplit(line, ' ');
  if (parts.size() != 5 || parts[0] != "p") {
    return absl::InvalidArgumentError("Invalid first p-line format, part size: " +
                                      std::to_string(parts.size()));
  }
  if (!(absl::SimpleAtoi(parts[2], &n0) &&   // Tries to turn the 2nd and 3rd "words" into
        absl::SimpleAtoi(parts[3], &n1))) {  // integers. If successfull, carry on

    return absl::InvalidArgumentError("Invalid n0 or n1 in the p-line");
  }
  if (!(n0 > 0 && n1 > 0)) {  // Check for empty graph
    return absl::InvalidArgumentError("n0 and n1 must be greater than 0");
  }
  if (!absl::SimpleAtoi(parts[4], &m) ||
      m < 0) {  // Try to turn the 4th word into integer. If successfull, carry on
    return absl::InvalidArgumentError("Invalid m in the p-line");
  }
  // Now create Data Structure
  auto bipartiteGraph =
      std::make_unique<BipartiteGraph>(n0, n1, m);  // BipartiteGraph needs a constructor with
                                                    // n0=numFixedNodes, n1=numFreeNodes, m=numEdges
  typename BipartiteGraph::NodeType source, target;
  const typename BipartiteGraph::NodeType size = n0 + n1;
  while (stream && stream >> target &&
         stream >> source) {  // Fixed Node appears first, free Node second
    if ((target < 1 || target > n0) ||
        (source <= n0 || source > size)) {  // Check if nodes out of bounds
      return absl::InvalidArgumentError("Invalid edge vertex indices");
    }

    bipartiteGraph->addEdge(source - 1, target - 1);

    // adds an edge going from Source to Target
    // internally we want 0-indexed arrays but the input is 1-indexed
  }
  // TODO: implement finish on data structures
  // bipartiteGraph->finish();
  return bipartiteGraph;
}