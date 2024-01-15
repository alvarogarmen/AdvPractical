cc_library(
    name = "algorithms",
    hdrs = [
        "reduction_algorithm.h",
    ],
    srcs = [
        "reduction_algorithm.cc",],
    visibility = ["//visibility:public"],
    deps = [
        "//ds/reduction_graph",
        "@abseil-cpp//absl/status:statusor",
        "@henrixapp_toolkit//toolkit/algorithms",
        "@henrixapp_toolkit//toolkit/app:app_io_cc_proto",
        "@henrixapp_toolkit//toolkit/ds:empty_problem",
        "//io:io",
    ],
    alwayslink=1,
)