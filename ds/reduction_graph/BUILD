load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "reduction_graph",
    hdrs = [
        "reduction_graph.h",
        "undo_algorithm_step.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":bit_vector",
        #"//external_hash/unordered_dense/include/ankerl",
        "@abseil-cpp//absl/status",
        #"@com_google_absl//absl/container:flat_hash_map",
        #"@com_google_absl//absl/container:node_hash_map",
    ],
)

cc_test(
    name = "reduction_graph_test",
    srcs = ["reduction_graph_test.cc"],
    deps = [
        ":reduction_graph",
        "@googletest//:gtest_main",
    ],
)


cc_library(
    name = "bit_vector",
    hdrs = [
        "bit_vector.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@abseil-cpp//absl/status",
    ],  # <-- You were missing a closing parenthesis here
)

cc_test(
    name = "bit_vector_test",
    srcs = ["bit_vector_test.cc"],
    deps = [
        ":bit_vector",
        "@googletest//:gtest_main",
    ],
)