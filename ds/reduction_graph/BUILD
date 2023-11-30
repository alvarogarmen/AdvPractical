load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

# Additional functionality that depends on the core
cc_library(
    name = "UndoAlgorithmStep",
    hdrs = [
        "UndoAlgorithmStep.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_library(
    name = "reduction_graph",
    hdrs = [
        "reduction_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":UndoAlgorithmStep",
    ],
)

cc_test(
    name = "reduction_graph_test",
    srcs = ["reduction_graph_test.cc"],
    deps = [
        ":UndoAlgorithmStep",
        ":reduction_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
