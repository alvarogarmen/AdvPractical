load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "median_reduction_algorithm",
    hdrs = [
        "median_reduction_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
        "//oscm/median_algorithm",
        "//oscm/reduction_algorithm",
    ],
)

cc_test(
    name = "median_reduction_algorithm_test",
    srcs = ["median_reduction_algorithm_test.cc"],
    deps = [
        ":median_reduction_algorithm",
        "//ds:bipartite_graph",
        "//ds/reduction_graph",
        "@googletest//:gtest_main",
    ],
)