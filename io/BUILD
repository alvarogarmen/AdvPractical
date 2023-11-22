load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "io",
    hdrs = [
        "input_graph_manually.h",
        "ocsm_graph_reader.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "ocsm_io_test",
    srcs = ["ocsm_graph_reader_test.cc"],
    deps = [
        ":io",
        "//ds/fpt:fpt_graph",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
    ],
)
