# stan4bart honors dbarts's probit active-row mask: 0/1 weights on a binary
# response name rows in and out of the likelihood rather than a precision, so
# they must reach the BART sampler's active-row mask instead of being dropped.

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(60, TRUE, TRUE, TRUE)
df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

a <- as.double(seq_len(nrow(df)) %% 4L != 1L)

# the sampler stan4bart drives reports the mask its weights installed
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 weights = a,
                 cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6, seed = 11,
                 bart_args = list(n.trees = 5, keepTrees = TRUE))
expect_identical(stan4bart:::getBartSampler(fit)$activeRows, a)

# a masked row leaves the likelihood: substituting arbitrary labels at the
# inactive rows leaves every active row's draw bitwise, matching dbarts's own
# front-door test for the same rule
fitWithResponse <- function(y) {
  df2 <- df
  df2$y <- y
  stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df2,
           weights = a,
           cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6, seed = 11,
           bart_args = list(n.trees = 5))
}
fit.orig <- fitWithResponse(df$y)
y.flipped <- df$y
y.flipped[a == 0] <- 1 - y.flipped[a == 0]
fit.flipped <- fitWithResponse(y.flipped)

ev.orig    <- extract(fit.orig, "ev")
ev.flipped <- extract(fit.flipped, "ev")
expect_identical(ev.orig[a == 1, , drop = FALSE], ev.flipped[a == 1, , drop = FALSE])

# the mask is not a no-op: an unmasked (all-ones weight) fit draws differently
fit.unmasked <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                          weights = rep(1, nrow(df)),
                          cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6, seed = 11,
                          bart_args = list(n.trees = 5))
ev.unmasked <- extract(fit.unmasked, "ev")
expect_false(identical(ev.orig[a == 1, , drop = FALSE], ev.unmasked[a == 1, , drop = FALSE]))

# a weight vector that is not all 0/1 has no coherent probit likelihood and is
# refused with dbarts's own message, the same way dbarts() refuses it
expect_error(
  stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
           weights = runif(nrow(df), 0.5, 2),
           cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
           bart_args = list(n.trees = 5)),
  "probit models do not support weights other than 0 and 1"
)

rm(a, fit, fitWithResponse, fit.orig, fit.flipped, ev.orig, ev.flipped,
   fit.unmasked, ev.unmasked, y.flipped, df)
