load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "median_algorithm",
    hdrs = [
        "median_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "median_algorithm_test",
    srcs = ["median_algorithm_test.cc"],
    deps = [
        ":median_algorithm",
        "//ds:bipartite_graph",
        "//ds/reduction_graph",
        "@googletest//:gtest_main",
    ],
)