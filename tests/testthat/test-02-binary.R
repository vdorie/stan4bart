context("mstan4bart binary response")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, TRUE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)


test_that("extract matches fitted in causal setting model", {
  df.train <- df
  
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                    cores = 2, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                    bart_args = list(n.trees = 11),
                    treatment = z)
  
  
  samples.ev.train <- extract(fit)
  samples.ev.test  <- extract(fit, sample = "test")
  
  fitted.ev.train <- apply(samples.ev.train, 1, mean)
  fitted.ev.test  <- apply(samples.ev.test,  1, mean)
  
  expect_equal(fitted.ev.train, fitted(fit, "ev", "train"))
  expect_equal(fitted.ev.test,  fitted(fit, "ev", "test"))
  
})

test_that("nonlinearities are estimated well", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  
  rows.train <- seq_len(floor(0.8 * nrow(df)))
  rows.test  <- seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df))

  df.train <- df[rows.train,]
  df.test  <- df[rows.test,]
  
  mstan4bart_fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                               df.train,
                               verbose = -1L,
                               test = df.test)
  
  mstan4bart_fitted <- fitted(mstan4bart_fit, sample = "test")
  mstan4bart_dev <- -2 * mean(log(ifelse(df.test$y == 1, mstan4bart_fitted, 1 - mstan4bart_fitted)))
  
  glmer_control <- lme4::glmerControl(check.conv.grad     = "ignore",
                                      check.conv.singular = "ignore",
                                      check.conv.hess     = "ignore")
  glmer_fit <- lme4::glmer(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + z + (1 + X4 | g.1) + (1 | g.2),
                           df.train, family = binomial(link = "probit"),
                           control = glmer_control)
  
  glmer_fitted <- predict(glmer_fit, newdata = df.test, type = "response")
  glmer_dev <- -2 * mean(log(ifelse(df.test$y == 1, glmer_fitted, 1 - glmer_fitted)))
  
  
  base_bart_fit <- bart2(y ~ ., df.train, test = df.test, verbose = FALSE,
                         n.samples = 1000, n.burn = 1000, rngSeed = 0)
  
  base_bart_fitted <- fitted(base_bart_fit, sample = "test")
  base_bart_dev <- -2 * mean(log(ifelse(df.test$y == 1, base_bart_fitted, 1 - base_bart_fitted)))
  
  
  rbart_fit <- rbart_vi(y ~ . - g.2, df.train, test = df.test, group.by = g.2,
                        group.by.test = df.test$g.2, verbose = FALSE,
                        n.samples = 1000, n.burn = 1000)
  
  rbart_fitted <- fitted(rbart_fit, sample = "test")
  rbart_dev <- -2 * mean(log(ifelse(df.test$y == 1, rbart_fitted, 1 - rbart_fitted)))
 
  expect_true(mstan4bart_dev <= glmer_dev)
  # low sample size, so we cut ourselves some slack
  expect_true(mstan4bart_dev <= 1.35 * base_bart_dev)
  expect_true(mstan4bart_dev <= 1.35 * rbart_dev)
})

test_that("predict matches supplied data", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
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
  
  samples.pred <- predict(fit, df.train, type = "indiv.ranef")
  samples.ev   <- extract(fit, "indiv.ranef", "train")
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
  
  samples.pred <- predict(fit, df.test, type = "indiv.ranef")
  samples.ev   <- extract(fit, "indiv.ranef", "test")
  expect_equal(samples.pred, samples.ev)
   
  samples.pred <- predict(fit, df.test, type = "ev")
  samples.ev   <- extract(fit, "ev", "test")
  expect_equal(samples.pred, samples.ev)
})

test_that("ppd has approximately right amount of noise", {
  df.train <- df
  
  set.seed(0)
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                    cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                    bart_args = list(n.trees = 11))
  
  samples.ev  <- extract(fit)
  samples.ppd <- extract(fit, type = "ppd")
  
  expect_true(mean((apply(samples.ev, 2, mean) - apply(samples.ppd, 2, mean))^2) <= 1e-1)
  
  
  r <- (samples.ev - samples.ppd) / sqrt(samples.ev * (1 - samples.ev))
  r <- mean(r[!is.nan(r)])
  expect_true(abs(r) <= 0.01)
})

