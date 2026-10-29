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

expect_error(extract(fit.norandom, "k"), "not fit with leaf sensitivity")
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

# an ignored argument is named, and a dbarts spelling is said to be one ------

named_warning <- function(expr) {
  msg <- NULL
  withCallingHandlers(expr, warning = function(w) {
    msg <<- c(msg, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  msg
}

msg <- named_warning(extract(fit.full, "ev", foo = 1, bar = 2))
expect_true(length(msg) == 1L && grepl("unused arguments ignored: 'foo', 'bar'", msg, fixed = TRUE))
expect_false(grepl("spelling", msg, fixed = TRUE))

msg <- named_warning(extract(fit.full, "ev", combineChains = FALSE))
expect_true(length(msg) == 1L && grepl("'combineChains'", msg, fixed = TRUE))
expect_true(grepl("'combineChains' is dbarts's spelling of 'combine_chains'", msg, fixed = TRUE))

msg <- named_warning(residuals(fit.full, 1))
expect_true(grepl("1 unnamed", msg, fixed = TRUE))

msg <- named_warning(predict(fit.trees, df, combineChains = FALSE, foo = 1))
expect_true(grepl("'combineChains', 'foo'", msg, fixed = TRUE))
expect_true(grepl("spelling of 'combine_chains'", msg, fixed = TRUE))

expect_true(grepl("'foo'", named_warning(fitted(fit.full, foo = 1)), fixed = TRUE))
expect_true(grepl("'foo'", named_warning(residuals(fit.full, foo = 1)), fixed = TRUE))

# fitted() and residuals() have no 'combine_chains' to point a dbarts spelling at
for (msg in list(named_warning(fitted(fit.full, combineChains = FALSE)),
                 named_warning(residuals(fit.full, combineChains = FALSE))))
{
  expect_equal(msg, "unused arguments ignored: 'combineChains'")
}

# a vector 'type' names the types that are allowed ---------------------------

predict_types <- c("ev", "ppd", "indiv.fixef", "indiv.ranef", "indiv.bart")

err <- tryCatch(predict(fit.trees, df, type = c("ev", "ppd")), error = function(e) conditionMessage(e))
expect_true(grepl("'type' must be a single string", err, fixed = TRUE))
expect_true(all(vapply(predict_types, function(t) grepl(paste0("'", t, "'"), err, fixed = TRUE), TRUE)))
expect_false(grepl("'arg'", err, fixed = TRUE))
expect_error(predict(fit.trees, df, type = character(0)), "'type' must be a single string")

# the whole default, as a wrapper forwards it, still means the first type
expect_equal(predict(fit.trees, df, type = predict_types), predict(fit.trees, df, type = "ev"))
expect_equal(dim(predict(fit.trees, df, type = "ppd")), dim(predictions))

# the same for extract() and fitted(), whose default is the first type as well

extract_types <- c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef", "indiv.bart",
                   "sigma", "Sigma", "k", "varcount", "stan", "trees", "callback")
fitted_types <- setdiff(extract_types, "trees")

for (case in list(list(f = extract, types = extract_types), list(f = fitted, types = fitted_types)))
{
  err <- tryCatch(case$f(fit.full, type = c("ev", "ppd")), error = function(e) conditionMessage(e))
  expect_true(grepl("'type' must be a single string", err, fixed = TRUE))
  expect_true(all(vapply(case$types, function(t) grepl(paste0("'", t, "'"), err, fixed = TRUE), TRUE)))
  expect_false(grepl("'arg'", err, fixed = TRUE))
  expect_error(case$f(fit.full, type = character(0)), "'type' must be a single string")
}
expect_false(grepl("'trees'", tryCatch(fitted(fit.full, type = c("ev", "ppd")), error = conditionMessage), fixed = TRUE))

# the type left out, or the whole default forwarded, is the first type
expect_equal(extract(fit.full), extract(fit.full, type = "ev"))
expect_equal(extract(fit.full, type = extract_types), extract(fit.full, type = "ev"))
expect_equal(fitted(fit.full), fitted(fit.full, type = "ev"))
expect_equal(fitted(fit.full, type = fitted_types), fitted(fit.full, type = "ev"))
