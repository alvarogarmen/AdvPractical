cc_library(
    name = "algorithms",
    hdrs = [
        "reduction_algorithm.h",
        "barycenter_algorithm.h",
    ],
    srcs = [
        "reduction_algorithm.cc",
        "barycenter_algorithm.cc"],
    visibility = ["//visibility:public"],
    deps = [
        "//ds/reduction_graph",
        "//ds:bipartite_graph",
        "@abseil-cpp//absl/status:statusor",
        "@henrixapp_toolkit//toolkit/algorithms",
        "@henrixapp_toolkit//toolkit/app:app_io_cc_proto",
        "@henrixapp_toolkit//toolkit/ds:empty_problem",
        "//io:io",
        "//oscm/reduction_algorithm",
        "//oscm/barycenter_algorithm",
        "//ds/helper:empty_problem",
    ],
    alwayslink=1,
)