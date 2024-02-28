load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "barycenter_algorithm",
    hdrs = [
        "barycenter_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "barycenter_algorithm_test",
    srcs = ["barycenter_algorithm_test.cc"],
    deps = [
        ":barycenter_algorithm",
        "//ds:bipartite_graph",
        "//ds/reduction_graph",
        "@com_google_googletest//:gtest_main",
    ],
)