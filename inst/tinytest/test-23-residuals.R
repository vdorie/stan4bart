# residuals(): the observed response less the posterior mean expected value,
# on the rows fitted() reports

set.seed(23)
n <- 60L
df <- data.frame(x = rnorm(n), g = factor(sample(4L, n, TRUE)))
df$y <- df$x + as.integer(df$g) / 2 + rnorm(n)
df$yb <- as.integer(df$y > median(df$y))

fit <- stan4bart(y ~ bart(x) + (1 | g), df, chains = 1, iter = 40,
                 warmup = 20, verbose = -1L, seed = 1)
expect_equal(residuals(fit), df$y - fitted(fit, type = "ev"))
expect_warning(residuals(fit, type = "ppd"), "unused arguments ignored")

# a binary fit's expected value is the probability of a 1
fitBinary <- stan4bart(yb ~ bart(x) + (1 | g), df, chains = 1, iter = 40,
                       warmup = 20, verbose = -1L, seed = 1)
residualsBinary <- residuals(fitBinary)
expect_equal(residualsBinary, df$yb - fitted(fitBinary, type = "ev"))
expect_true(all(abs(residualsBinary) < 1))

# rows na.omit dropped are absent, as they are from fitted()
dfMissing <- df
dfMissing$y[3L] <- NA
fitMissing <- stan4bart(y ~ bart(x) + (1 | g), dfMissing, chains = 1,
                        iter = 40, warmup = 20, verbose = -1L, seed = 1)
expect_equal(length(residuals(fitMissing)), n - 1L)
