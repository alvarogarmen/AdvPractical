load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "bipartite_graph",
    hdrs = [
        "bipartite_graph.h",
    ],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "bipartite_graph_test",
    srcs = ["bipartite_graph_test.cc"],
    deps = [
        ":bipartite_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
