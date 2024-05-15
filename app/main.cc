//
// Created by alvar on 26/04/2023.
//

#include <chrono>
#include <iostream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ds/bipartite_graph.h"
#include "ds/heuristic_graph/heuristic_graph.h"
#include "ds/reduction_graph/reduction_graph.h"
#include "io/ocsm_graph_reader.h"
#include "oscm/heuristic_algorithm/heuristic_algorithm.h"
#include "oscm/reduction_algorithm/reduction_algorithm.h"
ABSL_FLAG(std::string, example, "Default value", "Helpful text");
//  bazel run app -- --example="bazel?"; ./bazel-bin/app/app
int main(int argc, char* argv[]) {
  absl::ParseCommandLine(argc, argv);
  auto example = absl::GetFlag(FLAGS_example);  // Get the variable and store it
  std::cout << example << std::endl;
  if (argc != 2) {
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
  // auto result = readGraph<HeuristicGraph<int, int>>(input);
  auto result = readGraph<ReductionGraph<int, int>>(input);

  if (result.status() != absl::OkStatus()) {
    std::cerr << "result status in not ok: " << std::endl;
    return 1;  // Exit with an error code
  }

  auto graph = std::move(result.value());
  // std::cout << "the crossing is" << graph->getCrossings() << std::endl;
  // heuristic_algorithm::algorithm<HeuristicGraph<int, int>>(*graph, true, true, true);
  auto start_time = std::chrono::steady_clock::now();
  auto [crossingSum, orderVector] =
      reductionalgorithms::algorithm<ReductionGraph<int, int>, UndoAlgorithmStep<int, int>>(*graph);
  std::cout << "the reduction solution has " << crossingSum << " crossings" << std::endl;
  auto end_time = std::chrono::steady_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

  // Convert the elapsed time to minutes
  auto elapsed_minutes = std::chrono::duration_cast<std::chrono::minutes>(elapsed_time);
  std::cout << "Execution time: " << elapsed_minutes.count() << " minutes" << std::endl;
  // Create an output file stream
  /*
        std::ofstream outputFile(argv[2]);

        // Check if the file is open
        if (!outputFile.is_open()) {
          std::cerr << "Error opening file: output.txt" << std::endl;
          return 1;  // Return an error code
        }
        if (auto status = graph->writeResultsToFile(outputFile, orderVector); !status.ok()) {
          std::cerr << "Error: " << status << std::endl;
          return 1;
        }
        */

  // Close the file stream
  // outputFile.close();

  // std::cout << "the solution has " << graph->getCrossings() << " crossings" << std::endl;

  // auto myGraph = inputGraphManually<int>();
  // std::cout << "Crossings: " << crossGrader<decltype(myGraph)>(myGraph) << std::endl;

  return 0;
}