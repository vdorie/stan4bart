# weights, offset/offset_type/offset_test, the end-node sensitivity ('k')
# parameter, and subset

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(80, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

# weights ------------------------------------------------------------------

set.seed(7)
w <- runif(nrow(df), 0.5, 2)

fit.w <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                   weights = w,
                   cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                   bart_args = list(n.trees = 3))
expect_equal(as.vector(fit.w$weights), w)
expect_true("(weights)" %in% colnames(fit.w$frame))

expect_error(stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                       weights = rep(-1, nrow(df)),
                       cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                       bart_args = list(n.trees = 3)),
             "Negative weights")

expect_error(stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                       weights = as.character(rep(1, nrow(df))),
                       cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                       bart_args = list(n.trees = 3)),
             "numeric vector")

# heterogeneous weights rescale the ppd noise: observation-level ppd
# variance should be sigma^2 / weight (regression coverage for the weighted
# noise branch of extract(..., "ppd"))
if (at_home()) {
  set.seed(7)
  w.hetero <- c(rep(4, floor(nrow(df) / 2)), rep(0.25, ceiling(nrow(df) / 2)))

  fit.hetero <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                          weights = w.hetero,
                          cores = 1, verbose = -1L, chains = 1, warmup = 50, iter = 350,
                          bart_args = list(n.trees = 15))
  ev  <- extract(fit.hetero, "ev")
  ppd <- extract(fit.hetero, "ppd")
  resid.var <- apply(ppd - ev, 1, var)
  sigma <- fitted(fit.hetero, "sigma")

  ratio <- mean(resid.var * w.hetero / sigma^2)
  expect_true(ratio > 0.7 && ratio < 1.4)
}

# offset / offset_type -------------------------------------------------------

off <- rnorm(nrow(df))

# offset_type = "fixef"/"ranef" substitute that component's contribution
# with the offset in ev/ppd (regression coverage: these two used to be
# silently dropped due to a string mismatch against the internally-checked
# values "fixed"/"random")
for (ot in c("bart", "fixef", "ranef", "parametric", "default")) {
  fit.o <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                     offset = off, offset_type = ot,
                     cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                     bart_args = list(n.trees = 3))

  ev <- extract(fit.o, "ev")
  indiv.fixef <- extract(fit.o, "indiv.fixef")
  indiv.ranef <- extract(fit.o, "indiv.ranef")
  indiv.bart  <- extract(fit.o, "indiv.bart")

  expected <- switch(ot,
    bart       = sweep(indiv.fixef + indiv.ranef, 1, off, "+"),
    fixef      = sweep(indiv.ranef + indiv.bart,  1, off, "+"),
    ranef      = sweep(indiv.fixef + indiv.bart,  1, off, "+"),
    parametric = sweep(indiv.bart, 1, off, "+"),
    default    = sweep(indiv.fixef + indiv.ranef + indiv.bart, 1, off, "+")
  )
  expect_equal(unname(ev), unname(expected), info = paste("offset_type =", ot))
}

# offset_test supplies a separate offset applied only to test predictions
df.train <- df[seq_len(floor(0.75 * nrow(df))),]
df.test  <- df[seq.int(floor(0.75 * nrow(df)) + 1L, nrow(df)),]
off.test <- rnorm(nrow(df.test))

fit.ot <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                    df.train, test = df.test, offset_test = off.test,
                    cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                    bart_args = list(n.trees = 3))
expect_equal(as.vector(fit.ot$test$offset), off.test)

ev.test <- extract(fit.ot, "ev", sample = "test")
ev.test.components <- extract(fit.ot, "indiv.fixef", sample = "test") +
  extract(fit.ot, "indiv.ranef", sample = "test") + extract(fit.ot, "indiv.bart", sample = "test")
expect_equal(unname(ev.test - ev.test.components), unname(matrix(off.test, nrow(df.test), ncol(ev.test))))

# end-node sensitivity parameter 'k' -----------------------------------------

n.chains <- 2L; n.warmup <- 2L; n.iter <- 5L
fit.k <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                   cores = 1, verbose = -1L, chains = n.chains, warmup = n.warmup, iter = n.iter,
                   bart_args = list(n.trees = 3, k = chi(1.25, Inf)))

k.combined <- extract(fit.k, "k")
expect_equal(length(k.combined), (n.iter - n.warmup) * n.chains)

k.uncombined <- extract(fit.k, "k", combine_chains = FALSE)
expect_equal(dim(k.uncombined), c(n.iter - n.warmup, n.chains))
expect_true(all(k.uncombined > 0))

# subset ----------------------------------------------------------------

fit.subset <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                        subset = X1 > median(X1),
                        cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                        bart_args = list(n.trees = 3))
expect_equal(nrow(fit.subset$frame), sum(df$X1 > median(df$X1)))
