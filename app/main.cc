//
// Created by alvar on 26/04/2023.
//

#include <atomic>
#include <chrono>
#include <csignal>
#include <iostream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ds/bipartite_graph.h"
#include "ds/heuristic_graph/heuristic_graph.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/heuristic_algorithm/heuristic_algorithm.h"
ABSL_FLAG(std::string, example, "Default value", "Helpful text");
//  bazel run app -- --example="bazel?"; ./bazel-bin/app/app

std::atomic<bool> terminationRequested(false);

void signalHandler(int signum) {
  if (signum == SIGTERM) {
    terminationRequested.store(true);
  }
}
int main(int argc, char* argv[]) {
  std::signal(SIGTERM, signalHandler);
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
  auto result = readGraph<HeuristicGraph<int, int>>(input);

  if (result.status() != absl::OkStatus()) {
    std::cerr << "result status in not ok: " << std::endl;
    return 1;  // Exit with an error code
  }

  auto graph = std::move(result.value());
  HeuristicGraph<int, int>& myGraph = *graph;
  // Create an output file stream
  std::ofstream outputFile(argv[2]);
  heuristic_algorithm::heuristicAlgorithm(myGraph, true, true, true, terminationRequested);
  auto solution = myGraph.getPermutation();
  // Check if the file is open
  for (size_t i = 0; i < solution.size(); i++) {
    outputFile << solution[i] << std::endl;
  }
  // Close the file stream
  outputFile.close();

  // auto myGraph = inputGraphManually<int>();
  // std::cout << "Crossings: " << crossGrader<decltype(myGraph)>(myGraph) << std::endl;

  return 0;
}