load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "io",
    hdrs = [
        "input_graph_manually.h",
        "readGraph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
    ],
)

cc_test(
    name = "io_test",
    srcs = ["readGraph_test.cc"],
    deps = [
        ":io",
        "@com_google_googletest//:gtest_main",
    ],
)
