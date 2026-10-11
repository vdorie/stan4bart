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
w[c(2L, 11L)] <- 0

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

# all-zero weights leave no observation in the likelihood; refused up front on
# both response families rather than failing inside the lme4/glm init fits
expect_error(stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                       weights = rep(0, nrow(df)),
                       cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                       bart_args = list(n.trees = 3)),
             "all zero")
df.bin <- df
df.bin$y <- as.double(df$y > median(df$y))
expect_error(stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + (1 | g.1), df.bin,
                       weights = rep(0, nrow(df)),
                       cores = 1, verbose = -1L, chains = 1, warmup = 2, iter = 4,
                       bart_args = list(n.trees = 3)),
             "all zero")
rm(df.bin)

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
                     weights = w, offset = off, offset_type = ot,
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

# rows dropped for a missing value ---------------------------------------

# the weights and the offset are cut to the rows the model keeps, as the
# response and the designs are
local({
  set.seed(3)
  n <- 60L
  d <- data.frame(g = factor(sample(6L, n, TRUE)), x = rnorm(n), u = runif(n), v = runif(n),
                  w = runif(n, 0.5, 2), o = rnorm(n))
  d$y <- d$x + 3 * d$u + rnorm(6L)[d$g] + rnorm(n)
  d$u[c(4L, 17L)] <- NA # the forest's
  d$x[9L] <- NA         # a fixed column's
  d$y[30L] <- NA        # the response's
  d$w[2L] <- 0
  d$s <- !(seq_len(n) %in% c(5L, 11L, 23L, 41L, 55L))
  dropped <- c(4L, 9L, 17L, 30L)

  fit_rows <- function(use_subset, na.action, ...) {
    if (use_subset)
      stan4bart(y ~ x + bart(u + v) + (1 | g), d, weights = w, offset = o, subset = s,
                na.action = na.action, ...)
    else
      stan4bart(y ~ x + bart(u + v) + (1 | g), d, weights = w, offset = o,
                na.action = na.action, ...)
  }

  for (use_subset in c(FALSE, TRUE)) {
    rows <- if (use_subset) which(d$s) else seq_len(n)
    kept <- setdiff(rows, dropped)
    info <- paste("subset:", use_subset)

    fit <- fit_rows(use_subset, na.omit, iter = 0L)
    expect_identical(c(length(fit$weights), length(fit$offset), nrow(fit$X), ncol(fit$reTrms$Zt),
                       nrow(fit$bartData@x), length(fit$bartData@y)),
                     rep(length(kept), 6L), info = info)
    expect_identical(fit$y, d$y[kept], info = info)
    expect_identical(as.vector(fit$weights), d$w[kept], info = info)
    expect_identical(as.vector(fit$offset), d$o[kept], info = info)

    # the fit runs, and what it reports is on the kept rows: the expected value
    # holds each row's own offset, and the predictive draws each row's own weight
    fit <- fit_rows(use_subset, na.omit, cores = 1L, verbose = -1L, chains = 1L, warmup = 3L,
                    iter = 7L, bart_args = list(n.trees = 5L))
    ev <- extract(fit, "ev")
    parts <- extract(fit, "indiv.fixef") + extract(fit, "indiv.ranef") + extract(fit, "indiv.bart")
    expect_equal(ev - parts, matrix(d$o[kept], length(kept), ncol(ev)), check.attributes = FALSE,
                 info = info)
    sigma <- extract(fit, "sigma")
    warned <- character()
    set.seed(11)
    ppd <- withCallingHandlers(extract(fit, "ppd"), warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    set.seed(11)
    noise <- suppressWarnings(rnorm(length(kept) * length(sigma), 0,
                                    rep(sigma, length(kept)) * rep(sqrt(1 / d$w[kept]), each = length(sigma))))
    # a row of weight zero has no predictive draw: its one warning is rnorm's
    expect_identical(warned, "NAs produced", info = info)
    expect_equal(ppd - ev, t(matrix(noise, length(sigma))), check.attributes = FALSE, info = info)

    # under na.exclude the dropped rows are reported as missing and the kept
    # ones as above
    fit <- fit_rows(use_subset, na.exclude, cores = 1L, verbose = -1L, chains = 1L, warmup = 3L,
                    iter = 7L, bart_args = list(n.trees = 5L))
    ev <- extract(fit, "ev")
    parts <- suppressWarnings(extract(fit, "indiv.fixef") + extract(fit, "indiv.ranef") +
                              extract(fit, "indiv.bart"))
    expect_identical(which(is.na(ev[, 1L])), match(intersect(dropped, rows), rows), info = info)
    expect_equal((ev - parts)[match(kept, rows), ], matrix(d$o[kept], length(kept), ncol(ev)),
                 check.attributes = FALSE, info = info)
    warned <- character()
    ppd <- withCallingHandlers(extract(fit, "ppd"), warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_identical(warned, "NAs produced", info = info)
    expect_identical(which(is.na(ppd[, 1L])),
                     sort(match(c(2L, intersect(dropped, rows)), rows)), info = info)
  }
})
