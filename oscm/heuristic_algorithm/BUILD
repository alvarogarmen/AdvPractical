load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "heuristic_algorithm",
    hdrs = [
        "heuristic_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds/heuristic_graph",
    ],
)

cc_test(
    name = "heuristic_algorithm_test",
    srcs = ["heuristic_algorithm_test.cc"],
    deps = [
        ":heuristic_algorithm",
        "@googletest//:gtest_main",
    ],
)
