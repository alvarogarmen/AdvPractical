//
// Created by alvar on 26/04/2023.
//

#include <chrono>
#include <iostream>
#include <string>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "ds/bipartite_graph.h"
#include "io/input_graph_manually.h"
#include "oscm/grader.h"

ABSL_FLAG(std::string, example, "Default value", "Helpful text");
// bazel run app -- --example="bazel?"; ./bazel-bin/app/app
int main(int argc, char* argv[]) {
  absl::ParseCommandLine(argc, argv);
  auto example = absl::GetFlag(FLAGS_example);  // Get the variable and store it
  std::cout << example << std::endl;
  // Call out functions
  std::cout << "Something";

  auto myGraph = inputGraphManually<int>();
  std::cout << "Crossings: " << crossGrader<decltype(myGraph)>(myGraph) << std::endl;

  return 0;
}