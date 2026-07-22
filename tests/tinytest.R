if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("stan4bart", at_home = isTRUE(as.logical(Sys.getenv("NOT_CRAN"))))
} else {
  cat("package 'tinytest' not available; cannot run unit tests\n")
}
