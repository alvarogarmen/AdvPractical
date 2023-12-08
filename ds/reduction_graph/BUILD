load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_graph",
    hdrs = [
        "UndoAlgorithmStep.h",
        "reduction_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_test(
    name = "reduction_graph_test",
    srcs = ["reduction_graph_test.cc"],
    deps = [
        ":reduction_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
