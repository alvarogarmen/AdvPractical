load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "app",
    srcs = ["main.cc"],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "//io",
        "//oscm",
        "//toolkit/toolkit/app",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/flags:parse",
    ],
)

cc_binary(
    name = "runner",
    deps = [
        "//toolkit/toolkit/runner",
        "//algos/sample_problem",
        #add additional families here
    ],
)
