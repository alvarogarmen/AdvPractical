load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_algorithm",
    hdrs = [
        "r_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds/reduction_graph",
    ],
)

cc_test(
    name = "r_algorithm_test",
    srcs = ["r_algorithm_test.cc"],
    deps = [
        ":reduction_algorithm",
        "@com_google_googletest//:gtest_main",
    ],
)
