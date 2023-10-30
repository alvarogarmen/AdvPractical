load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "fpt_node",
    hdrs = [
        "fpt_node.h",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "fpt_graph",
    hdrs = [
        "fpt_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":fpt_node",
    ],
)

cc_test(
    name = "graph_test",
    srcs = ["graph_test.cc"],
    deps = [
        ":fpt_graph",
        ":fpt_node",
        "@com_google_googletest//:gtest_main",
    ],
)
