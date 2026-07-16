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

test_that("default fit stores no warmup draws but keeps tuning summaries", {
  # slim parametric store: only the two LIVE diagnostics (lp__, stepsize__) and
  # the transformed block (aux, beta, b, theta_L); no raw rows, no placeholders
  stan_rows <- dimnames(fit$stan)$parameters
  expect_equal(fit$par_names$diagnostic, c("lp__", "stepsize__"))
  expect_length(fit$par_names$upar, 0L)
  expect_false(any(grepl("^(z_beta|z_b|z_T|rho|zeta|tau|aux_unscaled)\\.", stan_rows)))
  expect_false(any(c("accept_stat__", "treedepth__", "n_leapfrog__", "divergent__", "energy__") %in% stan_rows))
  expect_true(all(c("lp__", "stepsize__", "aux.1") %in% stan_rows))
  expect_true(any(startsWith(stan_rows, "beta.")) && any(startsWith(stan_rows, "b.")) &&
              any(startsWith(stan_rows, "theta_L.")))

  # no full warmup, but adaptation summaries + a thinned trace are present
  expect_null(fit[["warmup"]])
  expect_false(is.null(fit$adaptation))
  n_chains <- dim(fit$stan)[3L]
  expect_length(fit$adaptation$step_size, n_chains)
  expect_equal(dim(fit$adaptation$inv_mass)[2L], n_chains)
  expect_equal(dim(fit$adaptation$snapshot), c(length(stan_rows), n_chains))
  expect_equal(dim(fit$adaptation$trace$stan), c(length(stan_rows), 7L, n_chains))
  expect_equal(dim(fit$adaptation$trace$sigma), c(7L, n_chains))
})

test_that("include_warmup errors informatively when warmup was not saved", {
  expect_error(extract(fit, "fixef", include_warmup = TRUE), "save_warmup = FALSE")
  expect_error(extract(fit, "fixef", include_warmup = "only"), "save_warmup = FALSE")
  expect_error(as.array(fit, include_warmup = TRUE), "save_warmup = FALSE")
})

test_that("extract include_samples works correctly with save_warmup = TRUE", {
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                   cores = 2, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11), save_warmup = TRUE,
                   treatment = z)
  expect_false(is.null(fit$warmup$stan))

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

test_that("save_raw_parameters restores the raw unconstrained rows", {
  fit_raw <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                       cores = 1, verbose = -1L, chains = 2, warmup = 7, iter = 13,
                       bart_args = list(n.trees = 11), seed = 0,
                       stan_args = list(save_raw_parameters = TRUE))
  raw_rows <- dimnames(fit_raw$stan)$parameters
  expect_true(any(startsWith(raw_rows, "z_beta.")))
  expect_true(all(c("accept_stat__", "divergent__") %in% raw_rows))
  expect_gt(length(fit_raw$par_names$upar), 0L)
  # extract(fit, "stan") surfaces the raw rows behind the opt-in
  expect_true(any(startsWith(rownames(extract(fit_raw, "stan")), "z_beta.")))
  # the transformed draws are identical to the slim store (storage is neutral)
  expect_equal(unname(extract(fit_raw, "fixef", combine_chains = FALSE)),
               unname(extract(fit_raw, "stan", combine_chains = FALSE)[startsWith(raw_rows, "beta."),,,drop=FALSE]))
})

test_that("extract varcount works", {
  varcounts <- extract(fit, "varcount")
  expect_true(!is.null(varcounts))
  expect_equal(dim(varcounts), c(9L, (13L - 7L) * 3L))
  varcounts <- extract(fit, "varcount", combine_chains = FALSE)
  expect_true(!is.null(varcounts))
  expect_equal(dim(varcounts), c(9L, (13L - 7L), 3L))
})

test_that("verbose works with different model specifications", {
  invisible(capture.output(capture.output(
    fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                     cores = 1, verbose = TRUE, chains = 1, warmup = 7, iter = 13,
                     bart_args = list(n.trees = 11)),
    type = "message"
  )))
  expect_is(fit, "stan4bartFit")
  invisible(capture.output(capture.output(
    fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + (1 + X4 | g.1) + (1 | g.2), df,
                     cores = 1, verbose = TRUE, chains = 1, warmup = 7, iter = 13,
                     bart_args = list(n.trees = 11)),
    type = "message"
  )))
  expect_is(fit, "stan4bartFit")
  invisible(capture.output(capture.output(
    fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z, df,
                     cores = 1, verbose = TRUE, chains = 1, warmup = 7, iter = 13,
                     bart_args = list(n.trees = 11)),
    type = "message"
  )))
  expect_is(fit, "stan4bartFit")
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
  expect_true(cor(indiv.ranef, ranef.true[ind.train]) >= 0.68)
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

test_that("predict accepts 'type' positionally on a fit with random effects", {
  # 'offset' used to sit ahead of 'type' in predict.stan4bartFit's formals;
  # a positionally-supplied 'type' (the conventional third argument, as in
  # predict(fit, newdata, "ev")) was silently captured by 'offset' instead.
  # For type %in% c("ev", "ppd") that character string then hit the
  # indiv.bart + indiv.fixef + indiv.ranef + offset sum and threw "non-numeric
  # argument to binary operator"; for the indiv.* types (which never use
  # offset) 'type' quietly reverted to its default "ev", so the wrong
  # component was returned with no error at all.
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]

  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   df.train,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11, keepTrees = TRUE))

  expect_equal(predict(fit, df.test, "ev"), predict(fit, df.test, type = "ev"))
  expect_equal(predict(fit, df.test, "indiv.bart"), predict(fit, df.test, type = "indiv.bart"))
  expect_equal(predict(fit, df.test, "indiv.fixef"), predict(fit, df.test, type = "indiv.fixef"))
  expect_equal(predict(fit, df.test, "indiv.ranef"), predict(fit, df.test, type = "indiv.ranef"))

  # a genuine offset, only reachable by name, still lands where documented
  off <- rep(0.5, nrow(df.test))
  expect_equal(predict(fit, df.test, type = "ev", offset = off),
              predict(fit, df.test, type = "ev") + off)
})

test_that("predict works with one chain", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   df.train,
                   test = df.test,
                   cores = 1, verbose = -1L, chains = 1, warmup = 7L, iter = 13L,
                   bart_args = list(n.trees = 11, keepTrees = TRUE))
  expect_is(fit, "stan4bartFit")
  predictions <- predict(fit, df.test)
  expect_equal(dim(predictions), c(nrow(df.test), 13L - 7L))
  expect_equal(names(dimnames(predictions)), c("observation", "iterations:chains"))
})

test_that("predict samples new levels of a random-slope grouping factor", {
  # levelfun (R/lme4_functions.R) sampled new-level coefficients by slicing
  # the per-draw covariance array with drop = FALSE, e.g. L[,,i,drop=FALSE];
  # for a random-slope block (more than one predictor, as with (1+X4|g.1))
  # that slice is a p x p x 1 array, not the p x p matrix chol() requires,
  # and base::chol() errored "'a' must be a square matrix". Intercept-only
  # blocks (p == 1) did not trigger it: chol()'s as.matrix() coercion
  # happens to collapse a 1 x 1 x 1 array into a valid 1 x 1 matrix.
  # Trigger: predict(..., newdata) with the default sample_new_levels = TRUE
  # when newdata has a grouping-factor level absent from training.
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  df.test$g.1[1L] <- max(df$g.1) + 1L # a g.1 level absent from df.train

  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   df.train,
                   cores = 1, verbose = -1L, chains = 1, warmup = 7L, iter = 13L,
                   bart_args = list(n.trees = 11, keepTrees = TRUE))

  predictions <- predict(fit, df.test, type = "indiv.ranef")
  expect_equal(dim(predictions), c(nrow(df.test), 13L - 7L))
  expect_true(all(is.finite(predictions)))
})

test_that("ppd has approximately right amount of noise", {
  skip_if_not_installed("lme4")
  df.train <- df
  
  set.seed(15)
  
  # the ev/ppd mean gap is sigma^2 / numSamples by construction, so the
  # first check needs sigma near its converged value; run long enough that
  # convergence is not trajectory luck for the pinned seed
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                   cores = 1, verbose = -1L, chains = 2, warmup = 100, iter = 200,
                   bart_args = list(n.trees = 25))
  
  
  samples.ev  <- extract(fit)
  samples.ppd <- extract(fit, type = "ppd")
  
  expect_true(mean((apply(samples.ev, 2, mean) - apply(samples.ppd, 2, mean))^2) <= 0.1)
  r <- mean(apply(samples.ev - samples.ppd, 2, sd)) / fitted(fit, "sigma")
  r <- max(r, 1 / r)
  expect_true(r <= 1.1)
})


test_that("data split in call itself based on variable in data frame", {
  df$train <- sample(c(rep(TRUE, floor(0.8 * nrow(df))), rep(FALSE, nrow(df) - floor(0.8 * nrow(df)))))
  df$X10 <- as.factor(LETTERS[1 + round(4 * df$X10, 0)])

  fit <- stan4bart(y ~ bart(. - y - X10) + X10, df[df$train,],
                   cores = 1, verbose = -1L, chains = 3, warmup = 1, iter = 3,
                   bart_args = list(n.trees = 4),
                   test = df[!df$train,])

  terms <- attr(fit$frame, "terms")

  expect_equal(attr(terms, "varnames.fixed"), "X10")
  expect_setequal(attr(terms, "varnames.bart"), c(
    "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "g.1", "g.2",
    "z", "train"))
  expect_equal(attr(terms, "varnames.random"), character())
  expect_equal(nrow(fit$frame), sum(df$train == TRUE))

  
  expect_setequal(colnames(fit$test$frame), colnames(df))
  expect_equal(nrow(fit$test$frame), sum(df$train == FALSE))
})

