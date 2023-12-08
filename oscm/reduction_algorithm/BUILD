load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_algorithm",
    hdrs = [
        "compute_crossings.h",
        "reduction_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds/reduction_graph",
    ],
)

cc_test(
    name = "reduction_algorithm_test",
    srcs = ["reduction_algorithm_test.cc"],
    deps = [
        ":reduction_algorithm",
        "@com_google_googletest//:gtest_main",
    ],
)
