context("stan4bart WALNUTS-era prior scope and control-arg deprecations")

# Q(b)/Q(c) of docs/plans/walnuts.md, landed at C5: the shrinkage coefficient
# priors are refused with an informative R-level error before any C++ is
# invoked (ParametricModel::finalize() would otherwise throw an uncaught
# C++ exception - a process abort, not a catchable R error), and the
# NUTS-specific stan_args controls are accepted-but-ignored with a one-time
# deprecation warning rather than silently doing nothing.

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(60, TRUE, FALSE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y))
rm(testData)

fit_call <- function(...) {
  stan4bart(y ~ bart(. - g.1 - g.2 - X4) + X4 + (1 + X4 | g.1) + (1 | g.2), df,
           cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
           bart_args = list(n.trees = 3), ...)
}

test_that("shrinkage-family coefficient priors are refused with an informative error", {
  msg <- "not supported by the gradient-based sampler"
  expect_error(fit_call(stan_args = list(prior = stan4bart:::hs())), msg)
  expect_error(fit_call(stan_args = list(prior = stan4bart:::hs_plus())), msg)
  expect_error(fit_call(stan_args = list(prior = stan4bart:::lasso())), msg)
  expect_error(fit_call(stan_args = list(prior = stan4bart:::laplace())), msg)
  expect_error(fit_call(stan_args = list(prior = stan4bart:::product_normal())), msg)
})

test_that("supported coefficient priors are unaffected", {
  expect_is(fit_call(stan_args = list(prior = stan4bart:::normal())), "stan4bartFit")
  expect_is(fit_call(stan_args = list(prior = stan4bart:::student_t())), "stan4bartFit")
  expect_is(fit_call(stan_args = list(prior = stan4bart:::cauchy())), "stan4bartFit")
})

test_that("NUTS control args are accepted with a deprecation warning, not ignored silently", {
  expect_warning(
    fit_call(stan_args = list(adapt_delta = 0.99)),
    "'adapt_delta' is ignored by the gradient-based sampler and is deprecated"
  )
  expect_warning(
    fit_call(stan_args = list(max_treedepth = 12, stepsize = 0.05)),
    "'max_treedepth' is ignored.*'stepsize' is ignored"
  )
})

test_that("loop-level and init_r args do not warn", {
  expect_silent(fit_call(stan_args = list(init_r = 1.5)))
})

test_that("a fit with no deprecated args produces no warning", {
  expect_silent(fit_call())
})
