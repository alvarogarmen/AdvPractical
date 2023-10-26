load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "app",
    srcs = ["main.cc"],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "//io",
        "//oscm",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/flags:parse",
    ],
)
