//
// Created by alvar on 26/04/2023.
//

#include <chrono>
#include <iostream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ds/bipartite_graph.h"
#include "ds/reduction_graph/reduction_graph.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/reduction_algorithm/reduction_algorithm.h"

ABSL_FLAG(std::string, example, "Default value", "Helpful text");
//  bazel run app -- --example="bazel?"; ./bazel-bin/app/app
int main(int argc, char* argv[]) {
  absl::ParseCommandLine(argc, argv);
  auto example = absl::GetFlag(FLAGS_example);  // Get the variable and store it
  std::cout << example << std::endl;
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input filename> <output filename>" << std::endl;
    return 1;  // Exit with an error code
  }

  // Open the file using std::ifstream
  std::ifstream input(argv[1]);

  // Check if the file is successfully opened
  if (!input.is_open()) {
    std::cerr << "Error opening file: " << argv[1] << std::endl;
    return 1;  // Exit with an error code
  }
  auto result = readGraph<ReductionGraph<int, int>>(input);

  if (result.status() != absl::OkStatus()) {
    std::cerr << "result status in not ok: " << std::endl;
    return 1;  // Exit with an error code
  }

  auto graph = std::move(result.value());
  auto [crossingSum, orderVector] =
      algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(*graph);
  // Create an output file stream
  std::ofstream outputFile(argv[2]);

  // Check if the file is open
  if (!outputFile.is_open()) {
    std::cerr << "Error opening file: output.txt" << std::endl;
    return 1;  // Return an error code
  }
  graph->writeResultsToFile(outputFile, orderVector);
  // Close the file stream
  outputFile.close();

  std::cout << "the solution has " << crossingSum << " crossings" << std::endl;

  // auto myGraph = inputGraphManually<int>();
  // std::cout << "Crossings: " << crossGrader<decltype(myGraph)>(myGraph) << std::endl;

  return 0;
}