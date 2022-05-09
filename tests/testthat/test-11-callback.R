context("stan4bart callback argument")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1 = as.factor(g.1), g.2 = as.factor(g.2), y, z))
df.train <- df[1:80,]
df.test  <- df[81:100,]

callback <- function(yhat.train, yhat.test, stan_pars) {
  rho <- parent.frame()
  if (is.null(rho$fixef.ind)) {
    rho$fixef.ind <- which(grepl("^beta|gamma", names(stan_pars)))
    rho$ranef.ind <- which(startsWith(names(stan_pars), "b."))
  }
  y.test <- rho$frame$y
  fixef <- stan_pars[rho$fixef.ind]
  ranef <- stan_pars[rho$ranef.ind]

  x_means <- rho$X_means
  keep_cols <- names(x_means) != "(Intercept)"
    
  intercept_delta <- sum(fixef[keep_cols] * x_means[keep_cols])
  
  fit.fixef <- as.vector(rho$test$X %*% fixef) - intercept_delta
  fit.ranef <- as.vector(Matrix::crossprod(rho$test$reTrms$Zt, ranef))
  yhat.test.full <- yhat.test + fit.fixef + fit.ranef
  
  c(yhat.test,
    fit.fixef,
    fit.ranef)
}
fn_env <- new.env(parent = baseenv())
environment(callback) <- fn_env

test_that("callback is passed sample correctly", {
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   data = df.train,
                   test = df.test,
                   cores = 1, verbose = -1L, chains = 2, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   seed = 0,
                   callback = callback)
  expect_is(fit, "stan4bartFit")

  indiv.bart <- fit$callback[seq_len(20L),,]
  indiv.fixef <- fit$callback[seq_len(20L) + 20L,,]
  indiv.ranef <- fit$callback[seq_len(20L) + 40L,,]

  expect_equal(unname(indiv.bart),
               unname(extract(fit, type = "indiv.bart",  sample = "test", combine_chains = FALSE)))
  expect_equal(unname(indiv.fixef),
               unname(extract(fit, type = "indiv.fixef", sample = "test", combine_chains = FALSE)))
  expect_equal(unname(indiv.ranef),
               unname(extract(fit, type = "indiv.ranef", sample = "test", combine_chains = FALSE)))
  
  ev <- indiv.fixef + indiv.ranef + indiv.bart
  expect_equal(unname(ev),
               unname(extract(fit, "ev", sample = "test", combine_chains = FALSE)))

  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   data = df.train,
                   test = df.test,
                   cores = 1, verbose = -1L, chains = 2, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   seed = 0, keep_fits = FALSE,
                   callback = callback)
  expect_is(fit, "stan4bartFit")
  
  indiv.bart <- fit$callback[seq_len(20L),,]
  indiv.fixef <- fit$callback[seq_len(20L) + 20L,,]
  indiv.ranef <- fit$callback[seq_len(20L) + 40L,,]
  
  expect_equal(unname(ev), unname(indiv.fixef + indiv.ranef + indiv.bart))
})


test_that("callback works with multiple threads", {
  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   data = df.train,
                   test = df.test,
                   cores = 2, verbose = -1L, chains = 2, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   seed = 0,
                   callback = callback)
  expect_is(fit, "stan4bartFit")

  
  indiv.bart <- fit$callback[seq_len(20L),,]
  indiv.fixef <- fit$callback[seq_len(20L) + 20L,,]
  indiv.ranef <- fit$callback[seq_len(20L) + 40L,,]

  expect_equal(unname(indiv.bart),
               unname(extract(fit, type = "indiv.bart",  sample = "test", combine_chains = FALSE)))
  expect_equal(unname(indiv.fixef),
               unname(extract(fit, type = "indiv.fixef", sample = "test", combine_chains = FALSE)))
  expect_equal(unname(indiv.ranef),
               unname(extract(fit, type = "indiv.ranef", sample = "test", combine_chains = FALSE)))
  
  ev <- indiv.fixef + indiv.ranef + indiv.bart
  expect_equal(unname(ev),
               unname(extract(fit, "ev", sample = "test", combine_chains = FALSE)))
})

