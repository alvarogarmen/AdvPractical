load("@rules_cc//cc:defs.bzl", "cc_library")

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
        ":oscm",
        ":readGraph",
        "//ds/fpt:fpt_graph"
        "@com_google_googletest//:gtest_main",
    ],
)