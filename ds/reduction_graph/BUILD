load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_graph",
    hdrs = [
        "reduction_graph.h",
        "undo_algorithm_step.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
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
