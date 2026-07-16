context("stan4bart print/summary and the mis-tuning diagnostic")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, FALSE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y))
rm(testData)

# default warmup ratio (iter %/% 2), comfortably above the documented floor,
# so this fit's mean_leapfrog should sit well under the warning threshold;
# seeded for a reproducible (non-flaky) diagnostic value
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4) + X4 + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = -1L, chains = 2, iter = 300,
                 bart_args = list(n.trees = 10), seed = 7)

test_that("a well-warmed fit's print/summary do not warn", {
  expect_true(max(fit$adaptation$mean_leapfrog) < 30)
  expect_message(print(fit), NA)
  expect_message(summary(fit), NA)
})

test_that("print and summary warn when mean_leapfrog is inflated", {
  bad_fit <- fit
  bad_fit$adaptation$mean_leapfrog[1] <- 50
  msg <- "parametric sampler tuning looks poor"
  expect_message(print(bad_fit), msg)
  expect_message(summary(bad_fit), msg)
})

test_that("a fit with no adaptation record does not error or warn", {
  no_adapt_fit <- fit
  no_adapt_fit$adaptation <- NULL
  expect_message(print(no_adapt_fit), NA)
  expect_message(summary(no_adapt_fit), NA)
})

test_that("summary returns a summary.stan4bartFit carrying the per-chain diagnostics", {
  s <- summary(fit)
  expect_true(inherits(s, "summary.stan4bartFit"))
  expect_equal(s$adaptation$mean_leapfrog, fit$adaptation$mean_leapfrog)
  expect_equal(s$adaptation$mean_leapfrog_warmup, fit$adaptation$mean_leapfrog_warmup)
})
