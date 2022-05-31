context("stan4bart with no random effects")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

fixef.true <- with(testData, mu.fixef.1 * z + mu.fixef.0 * (1 - z))
ranef.true <- with(testData, mu.ranef.1 * z + mu.ranef.0 * (1 - z))
bart.true  <- with(testData, mu.bart.1  * z + mu.bart.0  * (1 - z))

rm(testData)

test_that("fit issues warning if overparameterized", {
  expect_warning(stan4bart(y ~ X1 + bart(X1), df,
                           cores = 1L, verbose = 0L, chains = 1L, warmup = 0L, iter = 1L,
                           bart_args = list(n.trees = 1L)))
})


fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z, df,
                 cores = 1, verbose = -1L, chains = 3, warmup = 2, iter = 3,
                 bart_args = list(n.trees = 2),
                 treatment = z)

test_that("extract matches fitted in causal setting model for fixef only model", {
  samples.ev.train <- extract(fit)
  samples.ev.test  <- extract(fit, sample = "test")
  
  fitted.ev.train <- apply(samples.ev.train, 1, mean)
  fitted.ev.test  <- apply(samples.ev.test,  1, mean)
  
  expect_equal(fitted.ev.train, fitted(fit, "ev", "train"))
  expect_equal(fitted.ev.test,  fitted(fit, "ev", "test"))
})


test_that("extract matches as.array for fixef only model", {
  arr <- as.array(fit)
  fixef <- extract(fit, "fixef")
  
  expect_equal(c("sigma", dimnames(fixef)$predictor),
               dimnames(arr)$parameters)
})

test_that("extract include_samples works correctly for fixef only model", {
  sample <- extract(fit, "fixef", combine_chains = FALSE, include_warmup = FALSE)
  warmup <- extract(fit, "fixef", combine_chains = FALSE, include_warmup = "only")
  both   <- extract(fit, "fixef", combine_chains = FALSE, include_warmup = TRUE)
  
  expect_equal(dim(warmup)[2L] + dim(sample)[2L], dim(both)[2L])
  for (i in seq_len(dim(warmup)[1L])) {
    expect_equal(unname(rbind(warmup[i,,], sample[i,,])), unname(both[i,,]))
  }
  
  sample <- extract(fit, "sigma", combine_chains = FALSE, include_warmup = FALSE)
  warmup <- extract(fit, "sigma", combine_chains = FALSE, include_warmup = "only")
  both   <- extract(fit, "sigma", combine_chains = FALSE, include_warmup = TRUE)
  
  expect_equal(dim(warmup)[1L] + dim(sample)[1L], dim(both)[1L])
  expect_equal(unname(rbind(warmup, sample)), unname(both))
})



test_that("predict matches supplied data for fixef only model", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z,
                   df.train,
                   test = df.test,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11, keepTrees = TRUE))

  samples.pred <- predict(fit, df.train, type = "indiv.bart")
  samples.ev   <- extract(fit, "indiv.bart", "train")
  expect_equal(samples.pred, samples.ev)
  
  samples.pred <- predict(fit, df.train, type = "indiv.fixef")
  samples.ev   <- extract(fit, "indiv.fixef", "train")
  expect_equal(samples.pred, samples.ev)
    
  samples.pred <- predict(fit, df.train, type = "ev")
  samples.ev   <- extract(fit, "ev", "train")
  expect_equal(samples.pred, samples.ev)
  
  
  samples.pred <- predict(fit, df.test, type = "indiv.bart")
  samples.ev   <- extract(fit, "indiv.bart", "test")
  expect_equal(samples.pred, samples.ev)
  
  samples.pred <- predict(fit, df.test, type = "indiv.fixef")
  samples.ev   <- extract(fit, "indiv.fixef", "test")
  expect_equal(samples.pred, samples.ev)
  
  samples.pred <- predict(fit, df.test, type = "ev")
  samples.ev   <- extract(fit, "ev", "test")
  expect_equal(samples.pred, samples.ev)
})

