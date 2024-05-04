load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "barycenter_heuristic_algorithm",
    hdrs = [
        "barycenter_heuristic_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
        "//oscm/barycenter_algorithm",
        "//oscm/heuristic_algorithm",
    ],
)

cc_test(
    name = "barycenter_heuristic_algorithm_test",
    srcs = ["barycenter_heuristic_algorithm_test.cc"],
    deps = [
        ":barycenter_heuristic_algorithm",
        "//ds:bipartite_graph",
        "//ds/reduction_graph",
        "@googletest//:gtest_main",
    ],
)