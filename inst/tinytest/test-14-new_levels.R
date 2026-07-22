# predicting/extracting random effects for grouping-factor levels not seen
# during training (sample_new_levels), including levels that are entirely
# new for a multi-predictor (random slope + intercept) grouping factor

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(120, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

df.train <- df[1:90,]
df.test  <- df[91:120,]

# g.1 has a random slope and intercept (p = 2, correlated), g.2 has a random
# intercept only (p = 1); mark half the test rows with brand-new levels for
# both grouping factors, disjoint from anything in training
set.seed(11)
new.rows <- sort(sample(nrow(df.test), 15))
df.test$g.1[new.rows] <- df.test$g.1[new.rows] + 1000L
df.test$g.2[new.rows] <- df.test$g.2[new.rows] + 1000L

fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                 df.train, test = df.test,
                 cores = 1, verbose = -1L, chains = 2, warmup = 3, iter = 8,
                 bart_args = list(n.trees = 5, keepTrees = TRUE))

# sample_new_levels = FALSE: brand-new groups get exactly zero random effect,
# levels seen during training keep their normal, non-zero contribution
indiv.ranef.false <- extract(fit, "indiv.ranef", sample = "test", sample_new_levels = FALSE)
expect_equal(dim(indiv.ranef.false), c(nrow(df.test), (8L - 3L) * 2L))
expect_true(all(indiv.ranef.false[new.rows,] == 0))
expect_true(any(indiv.ranef.false[-new.rows,] != 0))

# sample_new_levels = TRUE: brand-new groups get draws from the posterior
# predictive distribution (non-zero, reproducible under a fixed seed), while
# previously-seen levels are unaffected
set.seed(3)
indiv.ranef.true1 <- extract(fit, "indiv.ranef", sample = "test", sample_new_levels = TRUE)
set.seed(3)
indiv.ranef.true2 <- extract(fit, "indiv.ranef", sample = "test", sample_new_levels = TRUE)

expect_equal(indiv.ranef.true1, indiv.ranef.true2)
expect_true(var(as.vector(indiv.ranef.true1[new.rows,])) > 0)
expect_equal(indiv.ranef.true1[-new.rows,], indiv.ranef.false[-new.rows,])

# predict() agrees with extract() for both settings
predict.false <- predict(fit, df.test, type = "indiv.ranef", sample_new_levels = FALSE)
expect_equal(unname(predict.false), unname(indiv.ranef.false))

set.seed(3)
predict.true <- predict(fit, df.test, type = "indiv.ranef", sample_new_levels = TRUE)
expect_equal(unname(predict.true), unname(indiv.ranef.true1))

# ev/ppd still compose correctly when some rows draw new-level random effects
ev <- extract(fit, "ev", sample = "test", sample_new_levels = FALSE)
indiv.fixef <- extract(fit, "indiv.fixef", sample = "test")
indiv.bart  <- extract(fit, "indiv.bart",  sample = "test")
expect_equal(unname(ev), unname(indiv.fixef + indiv.ranef.false + indiv.bart))

# statistical check (larger n, single posterior sample): draws for brand-new
# groups should reproduce the fitted Sigma for a random-slope grouping factor,
# not merely be non-zero (regression test for the levelfun sample_new_levels
# reconstruction: a prior version of this code errored outright for p > 1
# grouping factors, and separately sampled draws with the wrong covariance
# orientation whenever the random slope and intercept were correlated)
if (at_home()) {
  n.big <- 400L
  set.seed(21)
  big.train <- df[sample(nrow(df), n.big, replace = TRUE),]
  big.test  <- big.train
  big.test$g.1 <- seq_len(n.big) + 10000L  # every test row is a brand-new g.1 level
  big.test$X4  <- 0                        # isolate the intercept component

  fit.big <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                       big.train, test = big.test,
                       cores = 1, verbose = -1L, chains = 1, warmup = 5, iter = 6,
                       bart_args = list(n.trees = 5, keepTrees = TRUE))

  Sigma.g1 <- extract(fit.big, "Sigma", combine_chains = FALSE)[["g.1"]][,,1L,1L]
  ranef.new <- extract(fit.big, "indiv.ranef", sample = "test", sample_new_levels = TRUE,
                       combine_chains = FALSE)[,1L,1L]

  empirical.var <- var(ranef.new)
  target.var <- Sigma.g1["(Intercept)", "(Intercept)"]
  # loose tolerance: this is a Monte Carlo check with n.big draws, not a
  # bitwise one, but should be well within a factor of 2 of the target
  r <- empirical.var / target.var
  expect_true(r > 0.5 && r < 2)
}
