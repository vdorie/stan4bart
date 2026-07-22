# error and warning paths for extract / fitted / predict argument
# validation and type-availability checks

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(60, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

# model with both fixed and random effects, but no keepTrees and no callback
fit.full <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                      cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                      bart_args = list(n.trees = 3))

# model with fixed effects and bart only, no random effects
fit.norandom <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z, df,
                          cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                          bart_args = list(n.trees = 3))

# model with random effects and bart only, no fixed effects
fit.nofixed <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + (1 + X4 | g.1) + (1 | g.2), df,
                         cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                         bart_args = list(n.trees = 3))

# model with a binary response, for the 'sigma' unavailability check
testDataBinary <- local({
  source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)
  generateFriedmanData(60, TRUE, TRUE, TRUE)
})
dfb <- with(testDataBinary, data.frame(x, g.1, g.2, y, z))

fit.binary <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), dfb,
                        cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                        bart_args = list(n.trees = 3))

# model with the end-node sensitivity parameter modeled
fit.k <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z, df,
                   cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                   bart_args = list(n.trees = 3, k = chi(1.25, Inf)))

# extract type availability -----------------------------------------------

expect_error(extract(fit.norandom, "ranef"), "no modeled parameters")
expect_error(extract(fit.norandom, "Sigma"), "no modeled parameters")
expect_error(extract(fit.nofixed,  "fixef"), "no unmodeled parameters")

expect_error(extract(fit.binary, "sigma"), "binary outcome model")

expect_error(extract(fit.norandom, "k"), "not fit with end-node sensitivity")
k <- extract(fit.k, "k")
expect_equal(length(k), (4L - 2L) * 1L)
k2 <- extract(fit.k, "k", combine_chains = FALSE)
expect_equal(dim(k2), c(4L - 2L, 1L))

expect_error(extract(fit.norandom, "trees"), "keepTrees")

expect_error(extract(fit.norandom, "callback"), "without callback function")

# extract argument validation ----------------------------------------------

expect_error(extract(fit.full, include_warmup = "bogus"))
expect_error(extract(fit.full, include_warmup = c(TRUE, FALSE)))
expect_error(extract(fit.full, include_warmup = 1L))

expect_warning(extract(fit.full, "ev", foo = 1), "unused arguments ignored")

# predict argument validation and type availability ------------------------

expect_error(predict(fit.full, df, combine_chains = "bogus"), "combine_chains")
expect_error(predict(fit.full, df, combine_chains = c(TRUE, FALSE)), "combine_chains")
expect_error(predict(fit.full, df, combine_chains = NA), "combine_chains")

expect_error(predict(fit.nofixed, df, type = "indiv.fixef"), "does not include fixed effect terms")
expect_error(predict(fit.norandom, df, type = "indiv.ranef"), "does not include random effect terms")

# fit.full was not built with keepTrees, so predict cannot reconstruct the
# bart component for new data
expect_error(predict(fit.full, df), "keepTrees")
expect_error(predict(fit.full, df, type = "ev"), "keepTrees")

expect_warning(predict(fit.norandom, df, foo = 1), "unused arguments ignored")

# a model built with keepTrees can predict without error
fit.trees <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z, df,
                       cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                       bart_args = list(n.trees = 3, keepTrees = TRUE))
predictions <- predict(fit.trees, df)
expect_equal(dim(predictions), c(nrow(df), 4L - 2L))
