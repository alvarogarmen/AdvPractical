load("@rules_cc//cc:defs.bzl", "cc_binary")

cc_binary(
    name = "app",
    srcs = ["main.cc"],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "//ds/reduction_graph",
        "//io",
        "//oscm/reduction_algorithm",
        "@abseil-cpp//absl/flags:flag",
        "@abseil-cpp//absl/flags:parse",
        "@abseil-cpp//absl/status:statusor",
    ],
)
