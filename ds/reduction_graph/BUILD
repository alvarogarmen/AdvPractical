load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

# Additional functionality that depends on the core
cc_library(
    name = "undo",
    hdrs = [
        "undo.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
    ],
)

cc_library(
    name = "reduction_graph",
    hdrs = [
        "r_graph.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":undo",
    ],
)

cc_test(
    name = "r_graph_test",
    srcs = ["r_graph_test.cc"],
    deps = [
        ":reduction_graph",
        ":undo",
        "@com_google_googletest//:gtest_main",
    ],
)
