load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "algorithms",
    srcs = [
        "barycenter_algorithm.cc",
        "heuristic_algorithm.cc",
        "median_algorithm.cc",
        "reduction_algorithm.cc",
        "barycenter_heuristic_algorithm.cc",

    ],
    hdrs = [
        "barycenter_algorithm.h",
        "heuristic_algorithm.h",
        "median_algorithm.h",
        "reduction_algorithm.h",
        "barycenter_heuristic_algorithm.h",
        "problems.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//ds:bipartite_graph",
        "//ds/helper:empty_problem",
        "//ds/heuristic_graph",
        "//ds/reduction_graph",
        "//io",
        "//oscm/barycenter_algorithm",
        "//oscm/heuristic_algorithm",
        "//oscm/barycenter_heuristic_algorithm",
        "//oscm/median_algorithm",
        "//oscm/reduction_algorithm",
        "@abseil-cpp//absl/status:statusor",
        "@henrixapp_toolkit//toolkit/algorithms",
        "@henrixapp_toolkit//toolkit/app:app_io_cc_proto",
        "@henrixapp_toolkit//toolkit/ds:empty_problem",
    ],
    alwayslink = 1,

)

