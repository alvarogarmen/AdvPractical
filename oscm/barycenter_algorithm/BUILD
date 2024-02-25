load("@rules_cc//cc:defs.bzl", "cc_library", "cc_test")

cc_library(
    name = "barycenter",
    hdrs = [
        "barycenter_algorithm.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "barycenter_test",
    srcs = ["barycenter_test.cc"],
    deps = [
        ":oscm/barycenter_algorithm",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
    ],
)