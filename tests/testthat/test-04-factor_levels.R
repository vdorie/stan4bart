context("stan4bart factor levels")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

df$X1 <- cut(df$X1, quantile(df$X1, seq(0.0, 1.0, 0.2)), include.lowest = TRUE)
levels(df$X1) <- letters[seq_len(nlevels(df$X1))]

df.train <- df[seq_len(floor(0.8 * nrow(df))),]
df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]

levels(df.train$X1) <- c(levels(df.train$X1), letters[1L:2L + nlevels(df.train$X1)])
levels(df.test$X1) <- c(levels(df.test$X1), letters[c(1L, 3L) + nlevels(df.test$X1)])

# keepTrees so predict() on new data can be exercised below; storage-only flag,
# does not alter sampling. X1 is a factor with unused ("empty") levels in both
# train and test, encoded under categorical splits as a single bart column.
fit <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) +
                       X4 + z + (1 + X4 | g.1) + (1 | g.2),
                 df.train, test = df.test,
                 cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                 bart_args = list(n.trees = 11, keepTrees = TRUE))

test_that("model fits with empty factor levels", {
  types <- c("ev", "indiv.fixef", "indiv.ranef", "indiv.bart")
  for (type in types) {
    value <- fitted(fit, type = type)

    expect_true(length(value) == nrow(df.train))
    expect_true(!anyNA(value))

    value <- fitted(fit, type = type, sample = "test")

    expect_true(length(value) == nrow(df.test))
    expect_true(!anyNA(value))
  }
})

test_that("varcount is named per factor, not per level", {
  # Under categorical splits a factor is one bart column (colnames(bartData@x)
  # is bare "X1", not one "X1.<level>" indicator per level), so varcount has
  # one row per factor instead of one row per level.
  varcount_names <- rownames(fit$bart_varcount)

  expect_identical(varcount_names, colnames(fit$bartData@x))
  expect_true("X1" %in% varcount_names)
  expect_false(any(grepl("^X1\\.", varcount_names)))
  expect_equal(nrow(fit$bart_varcount), 9L)
})

test_that("predict on fit-time test rows agrees with the fit-time test data", {
  # Both routes build the test design through the same makeTestModelMatrix
  # seam (getTestDataFrames / stan4bart.R), so a categorical factor column
  # should encode identically whether it arrives via `test =` at fit time or
  # via predict() afterwards on the same rows.
  samples.pred <- predict(fit, df.test, type = "indiv.bart")
  samples.ev   <- extract(fit, "indiv.bart", sample = "test")
  expect_equal(samples.pred, samples.ev)

  samples.pred <- predict(fit, df.test, type = "ev")
  samples.ev   <- extract(fit, "ev", sample = "test")
  expect_equal(samples.pred, samples.ev)
})

test_that("a new BART factor level errors informatively at predict", {
  df.newlevel <- df.test
  levels(df.newlevel$X1) <- c(levels(df.newlevel$X1), "z")
  df.newlevel$X1[1L] <- "z"

  expect_error(predict(fit, df.newlevel, type = "ev"), "new levels")
})

test_that("a new random-effects grouping level still predicts under both sample_new_levels settings", {
  # Unlike an unknown BART factor level (an error, tested above), an unknown
  # grouping level in the random-effects structure is not an error: lme4's
  # new-level machinery (levelfun, R/lme4_functions.R) is untouched by the
  # BART categorical-factor change and still draws from the prior or returns
  # zero-ev random effects, per sample_new_levels. This triggers the same
  # known warning as test-01's "predict samples new levels of a random-slope
  # grouping factor" test (levelfun slicing the per-draw covariance array for
  # the (1 + X4 | g.1) block) - expected, unrelated to this change, and
  # suppressed here since test-01 owns asserting it; this test only cares
  # that prediction succeeds.
  df.newgroup <- df.test
  df.newgroup$g.1[1L] <- max(df$g.1) + 1L

  for (snl in c(TRUE, FALSE)) {
    predictions <- suppressWarnings(
      predict(fit, df.newgroup, type = "indiv.ranef", sample_new_levels = snl)
    )
    expect_equal(dim(predictions), c(nrow(df.newgroup), (13L - 7L) * 3L))
    expect_true(all(is.finite(predictions)))
  }
})
