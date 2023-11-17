load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_graph",
    hdrs = [
        "r_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_test(
    name = "r_graph_test",
    srcs = ["r_graph_test.cc"],
    deps = [
        ":reduction_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
