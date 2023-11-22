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
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "io_test",
    srcs = ["readGraph_test.cc"],
    deps = [
        ":io",
        "@com_google_absl//absl/status",
        "@com_google_googletest//:gtest_main",
    ],
)
