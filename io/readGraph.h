#include <absl/status/statusor.h>
#include <absl/strings/str_split.h>

#include <fstream>
#include <iostream>
#include <string>

template <typename BipartiteGraph>
absl::StatusOr<std::unique_ptr<BipartiteGraph>> readGraph(const std::string& filename,
                                                          bool verbose = false,
                                                          bool readAll = true) {
  std::ifstream file(filename);
  if (!stream.good()) {
    return absl::NotFoundError("stream not good error");
  }
  std::string header;
  while (header == "") {
    if (stream.eof()) {
      return absl::NotFoundError("File does not contain enough lines. Could not read header.");
    }
    std::getline(stream, header);
    if (header.length() > 0 && header[0] == '%') {
      header = "";
    }
  }

  std::string line;
  size_t n0 = 0;
  size_t n1 = 0;
  size_t m = 0;
  bool pLineEncountered = false;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == 'c') {
      // Comment line, skip
      continue;
    }

    if (!pLineEncountered) {
      std::vector<std::string> parts = absl::StrSplit(line, ' ');
      if (parts.size() != 5 || parts[0] != "p") {
        return absl::InvalidArgumentError("Invalid first p-line format");
      }
      pLineEncountered = true;
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
      auto bipartiteGraph = std::make_unique<BipartiteGraph>(
          n0, n1, m);  // BipartiteGraph needs a constructor with n0=numFixedNodes, n1=numFreeNodes,
                       // m=numEdges
    } else {
      int source, target;
      if (sscanf(line.c_str(), "%d %d", &source, &target) !=
          2) {  // Tries to parse string as integer. Returns true if successfull
        return absl::InvalidArgumentError("Invalid edge format");
      }
      if ((source < 1 || source > n0) ||
          (target < n0 + 1 || target > n0 + n1)) {  // Check if nodes out of bounds
        return absl::InvalidArgumentError("Invalid edge vertex indices");
      }
      bipartiteGraph->addEdge(source, target);  // adds an edge going from Source to Target
    }
  }

  return absl::OkStatus();
}