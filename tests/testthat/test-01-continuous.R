context("mstan4bart continuous response")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

test_that("extract matches fitted in causal setting model", {
  df.train <- df
  
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                    cores = 2, verbose = 0, chains = 3, warmup = 7, iter = 13,
                    bart_args = list(n.trees = 11),
                    treatment = z)
  
  
  samples.ev.train <- extract(fit)
  samples.ev.test  <- extract(fit, sample = "test")
  
  # samples.icate <- (samples.ev.train - samples.ev.test) * (2 * df.train$z - 1)
  # samples.cate <- apply(samples.icate, 2, mean)
  # cate <- mean(samples.cate)
  
  # samples.ppd.cf <- extract(fit, type = "ppd", sample = "test")
  
  # samples.ite <- (df.train$y - samples.ppd.cf) * (2 * df.train$z - 1)
  # samples.sate <- apply(samples.ite, 2, mean)
  # sate <- mean(samples.sate)
  
  # samples.ppd.fac <- extract(fit, type = "ppd", sample = "train")
  # samples.pate <- apply((samples.ppd.fac - samples.ppd.cf) * (2 * df.train$z - 1), 2, mean)
  # pate <- mean(samples.pate)
  
  fitted.ev.train <- apply(samples.ev.train, 1, mean)
  fitted.ev.test  <- apply(samples.ev.test,  1, mean)
  
  expect_equal(fitted.ev.train, fitted(fit, "ev", "train"))
  expect_equal(fitted.ev.test,  fitted(fit, "ev", "test"))
  
})

test_that("nonlinearities are estimated well", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  bart_fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                         df.train,
                         test = df.test)
  
  bart_fitted <- fitted(bart_fit, sample = "test")
  bart_rmse <- sqrt(mean((df.test$y - bart_fitted)^2)) / sd(df.train$y)
  
  lmer_control <- lme4::lmerControl(check.conv.grad     = "ignore",
                                    check.conv.singular = "ignore",
                                    check.conv.hess     = "ignore")
  lmer_fit <- lme4::lmer(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + z + (1 + X4 | g.1) + (1 | g.2),
                         df.train,
                         control = lmer_control)
  
  lmer_fitted <- predict(lmer_fit, newdata = df.test, type = "response")
  # predict doesn't like the .
  lmer_rmse <- sqrt(mean((df.test$y - lmer_fitted)^2)) / sd(df.train$y)
  
  expect_true(bart_rmse <= lmer_rmse)
})

test_that("predict matches supplied data", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                    df.train,
                    test = df.test,
                    cores = 1, verbose = 0, chains = 3, warmup = 7, iter = 13,
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
                    cores = 1, verbose = 0, chains = 3, warmup = 7, iter = 13,
                    bart_args = list(n.trees = 11))
  
  
  samples.ev  <- extract(fit)
  samples.ppd <- extract(fit, type = "ppd")
  
  expect_true(mean((apply(samples.ev, 2, mean) - apply(samples.ppd, 2, mean))^2) <= 1e-1)
  r <- mean(apply(samples.ev - samples.ppd, 2, sd)) / fitted(fit, "sigma")
  r <- max(r, 1 / r)
  expect_true(r <= 1.1)
})


