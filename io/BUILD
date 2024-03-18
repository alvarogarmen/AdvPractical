load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "io",
    hdrs = [
        "input_graph_manually.h",
        "ocsm_graph_reader.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "//ds/heuristic_graph",
        "//ds/reduction_graph",
        "//oscm/reduction_algorithm",
        "@abseil-cpp//absl/status",
    ],
)

cc_test(
    name = "oscm_io_test",
    srcs = ["ocsm_graph_reader_test.cc"],
    deps = [
        ":io",
        "//ds/reduction_graph",
        "@abseil-cpp//absl/status:statusor",
        "@abseil-cpp//absl/strings",
        "@com_google_googletest//:gtest_main",
    ],
)
