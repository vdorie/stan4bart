context("stan4bart continuous response")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

fixef.true <- with(testData, mu.fixef.1 * z + mu.fixef.0 * (1 - z))
ranef.true <- with(testData, mu.ranef.1 * z + mu.ranef.0 * (1 - z))
bart.true  <- with(testData, mu.bart.1  * z + mu.bart.0  * (1 - z))

rm(testData)

fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 2, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                 bart_args = list(n.trees = 11),
                 treatment = z)

test_that("no intercept parameter was included", {
  expect_true(!("(Intercept)" %in% colnames(fit$X)))
  expect_true(!any(startsWith("gamma", dimnames(fit$stan)$parameters)))
})

test_that("extract matches fitted in causal setting model", {
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

test_that("extract matches as.array", {
  arr <- as.array(fit)
  ranef <- extract(fit, "ranef")
  
  grouping_factors <- names(fit$reTrms$cnms)
  for (grouping_factor in grouping_factors) {
    predictors <- fit$reTrms$cnms[[grouping_factor]]
    groups <- levels(fit$reTrms$flist[[grouping_factor]])
    for (predictor in predictors) {
      for (group in groups) {
        expect_true(all(arr[,,paste0("b[", predictor, " ", grouping_factor, ":", group, "]")] ==
                        ranef[[grouping_factor]][predictor,group,]))
      }
    }
  }
})

test_that("extract include_samples works correctly", {
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

test_that("nonlinearities are estimated well", {
  # Because this is not documented, to enable this test execute from R
  #   Sys.setenv(NOT_CRAN = "true")
  # or from shell
  #   export NOT_CRAN=true
  skip_on_cran()
  skip_if_not_installed("lme4")
  
  ind.train <- seq_len(floor(0.8 * nrow(df)))
  ind.test <- seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df))
  df.train <- df[ind.train,]
  df.test  <- df[ind.test,]
  
  bart_fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                        seed = 0,
                        df.train,
                        verbose = -1L,
                        test = df.test)
  
  bart_fitted <- fitted(bart_fit, sample = "test")
  bart_rmse <- sqrt(mean((df.test$y - bart_fitted)^2)) / sd(df.train$y)
  
  lmer_control <- lme4::lmerControl(check.conv.grad     = "ignore",
                                    check.conv.singular = "ignore",
                                    check.conv.hess     = "ignore")
  # predict doesn't like the .
  lmer_fit <- lme4::lmer(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + z + (1 + X4 | g.1) + (1 | g.2),
                         df.train,
                         control = lmer_control)
  
  lmer_fitted <- predict(lmer_fit, newdata = df.test, type = "response")
  lmer_rmse <- sqrt(mean((df.test$y - lmer_fitted)^2)) / sd(df.train$y)
  
  expect_true(bart_rmse <= lmer_rmse)
  
  indiv.bart <- fitted(bart_fit, type = "indiv.bart", sample = "train") 
  expect_true(cor(indiv.bart, bart.true[ind.train]) >= 0.95)
  indiv.ranef <- fitted(bart_fit, type = "indiv.ranef", sample = "train")
  expect_true(cor(indiv.ranef, ranef.true[ind.train]) >= 0.8)
  indiv.fixef <- fitted(bart_fit, type = "indiv.fixef", sample = "train")
  expect_true(cor(indiv.fixef, fixef.true[ind.train]) >= 0.99)
})

test_that("works when QR is true", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  
  ind.train <- seq_len(floor(0.8 * nrow(df)))
  ind.test <- seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df))
  df.train <- df[ind.train,]
  df.test  <- df[ind.test,]
  
  bart_fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                        seed = 0,
                        df.train,
                        verbose = -1L,
                        test = df.test,
                        cores = 4,
                        stan_args = list(QR = TRUE))
  
  bart_fitted <- fitted(bart_fit, sample = "test")
  bart_rmse <- sqrt(mean((df.test$y - bart_fitted)^2)) / sd(df.train$y)
  
  lmer_control <- lme4::lmerControl(check.conv.grad     = "ignore",
                                    check.conv.singular = "ignore",
                                    check.conv.hess     = "ignore")
  # predict doesn't like the .
  lmer_fit <- lme4::lmer(y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + z + (1 + X4 | g.1) + (1 | g.2),
                         df.train,
                         control = lmer_control)
  
  lmer_fitted <- predict(lmer_fit, newdata = df.test, type = "response")
  lmer_rmse <- sqrt(mean((df.test$y - lmer_fitted)^2)) / sd(df.train$y)
  
  expect_true(bart_rmse <= lmer_rmse)
  
  indiv.bart <- fitted(bart_fit, type = "indiv.bart", sample = "train") 
  expect_true(cor(indiv.bart, bart.true[ind.train]) >= 0.95)
  indiv.ranef <- fitted(bart_fit, type = "indiv.ranef", sample = "train")
  expect_true(cor(indiv.ranef, ranef.true[ind.train]) >= 0.8)
  indiv.fixef <- fitted(bart_fit, type = "indiv.fixef", sample = "train")
  expect_true(cor(indiv.fixef, fixef.true[ind.train]) >= 0.99)
})


test_that("predict matches supplied data", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
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
  
  set.seed(15)
  
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11))
  
  
  samples.ev  <- extract(fit)
  samples.ppd <- extract(fit, type = "ppd")
  
  expect_true(mean((apply(samples.ev, 2, mean) - apply(samples.ppd, 2, mean))^2) <= 0.1)
  r <- mean(apply(samples.ev - samples.ppd, 2, sd)) / fitted(fit, "sigma")
  r <- max(r, 1 / r)
  expect_true(r <= 1.1)
})


