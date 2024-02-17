load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "heuristic_graph",
    hdrs = [
        "heuristic_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_test(
    name = "heuristic_graph_test",
    srcs = ["heuristic_graph_test.cc"],
    deps = [
        ":heuristic_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
