Sys.unsetenv("R_TESTS")
if (require(testthat, quietly = TRUE)) {
  require(stan4bart)
  test_check("stan4bart")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}

