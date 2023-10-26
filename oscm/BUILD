load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "oscm",
    hdrs = [
        "grader.h",
    ],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "grader_test",
    srcs = ["grader_test.cc"],
    deps = [
        ":oscm",
        "//ds:bipartite_graph",
        "@com_google_googletest//:gtest_main",
    ],
)
