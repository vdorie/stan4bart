#!/usr/bin/env Rscript
# bart-train-recompute C0 EXACTNESS + LATENCY harness (docs/plans/
# bart-train-recompute.md C0). Fits with keep_fits = TRUE AND keepTrees =
# TRUE, then reconstructs the full stored `bart_train` block by replaying the
# kept trees over the TRAINING design matrix through dbarts's predict path
# (C_stan4bart_predictBART on fit$bartData@x - the tree-replay route the plan
# identifies at init.cpp:332-383). The gate is elementwise:
#
#     max |bart_train_stored - bart_train_recomputed| < TOL   (TOL = 1e-10)
#
# The plan's measured agreement is ~1e-13 (measurement 4); the two differ only
# by the stored value's offset round-trip (fl((offset+treeSum)-offset) vs
# treeSum, init.cpp:688-694), never by tree-replay error, so this is a TIGHT
# tolerance, not a loose one. NEVER loosen past 1e-10: a larger divergence
# means the replay path itself has drifted and must be re-derived (plan's C0
# abort clause).
#
# It also records the recompute wall time and the stored/kept-tree sizes per
# scale. One scale is the plan's past-the-old-threshold BLOCKER case
# (n = 50000 x 200 draws) that createStoredBARTSampler -> dbarts setState
# rejected before the C1 empty-leaf-veto fix (plan measurement 6); its
# assembling + recomputing at all is the functional proof C1 is in the linked
# dbarts.
#
# Usage (one group per invocation to stay inside a foreground time budget):
#   Rscript benchmarks/R/recompute-exactness.R          # fast tiers (small/medium/binary)
#   Rscript benchmarks/R/recompute-exactness.R big      # the n=50000 past-threshold case
#   Rscript benchmarks/R/recompute-exactness.R all      # everything

suppressMessages(library(stan4bart))

TOL <- 1e-10   # exactness gate; see header - NEVER loosen beyond this

# Friedman-style BART signal over X1..X5, matching inst/common/friedmanData.R
# and the posterior-baseline harness so BART's side is well identified.
friedman_f <- function(X)
  10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]

simulate <- function(seed, n, n_levels, binary = FALSE) {
  set.seed(seed)
  X <- matrix(runif(n * 5), n, 5, dimnames = list(NULL, paste0("X", 1:5)))
  Xfix <- runif(n, -1, 1)
  g <- factor(sample.int(n_levels, n, replace = TRUE))
  b <- rnorm(n_levels, 0, 1.2)
  mu <- friedman_f(X) + 2 * Xfix + b[g]
  if (binary) {
    loc <- mean(mu); scale <- sd(mu) / qnorm(0.15)
    mu <- (mu - loc) / scale
    y <- rbinom(n, 1L, pnorm(mu))
  } else {
    y <- mu + rnorm(n, 0, 1)
  }
  data.frame(X, Xfix = Xfix, g = g, y = y)
}

# Literal call sites: stan4bart's bart()/bar-term detection walks the
# UNEVALUATED formula from match.call(), so the formula must be written out at
# the call site (a formula held in a variable silently fails to strip bart()).
fit_scale <- function(df, binary, warmup, iter, chains, n.trees, seed) {
  if (binary)
    stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 | g), df,
              cores = 1L, verbose = -1L, chains = chains,
              warmup = warmup, iter = iter, seed = seed, keep_fits = TRUE,
              bart_args = list(n.trees = n.trees, keepTrees = TRUE))
  else
    stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 | g), df,
              cores = 1L, verbose = -1L, chains = chains,
              warmup = warmup, iter = iter, seed = seed, keep_fits = TRUE,
              bart_args = list(n.trees = n.trees, keepTrees = TRUE))
}

fmt_mb <- function(bytes) sprintf("%.1f MB", as.numeric(bytes) / 1024^2)

run_scale <- function(label, n, n_levels, warmup, iter, chains, n.trees,
                      seed = 99L, binary = FALSE) {
  cat(sprintf("\n=== %s: n=%d, %d draws x %d chain(s), n.trees=%d%s ===\n",
              label, n, iter - warmup, chains, n.trees,
              if (binary) ", binary" else ""))
  df <- simulate(seed, n, n_levels, binary)

  t_fit <- system.time(
    fit <- fit_scale(df, binary, warmup, iter, chains, n.trees, seed)
  )[["elapsed"]]

  # C1 proof: the stored-tree sampler assembled (this is the createStoredBART-
  # Sampler -> dbarts setState step that rejected the 50000x200 case pre-fix).
  if (is.null(fit$sampler.bart))
    stop(sprintf("[%s] FAIL: keepTrees fit did not assemble a stored BART sampler", label))

  stored <- fit$bart_train                       # n x S x C, the block we replace
  X.bart <- fit$bartData@x                        # the training design matrix

  t_rec <- system.time(
    recomputed <- .Call(stan4bart:::C_stan4bart_predictBART, fit$sampler.bart, X.bart, NULL)
  )[["elapsed"]]
  if (length(dim(recomputed)) == 2L)
    dim(recomputed) <- c(dim(recomputed), 1L)     # chains==1 arrives 2-D

  if (!all(dim(stored) == dim(recomputed)))
    stop(sprintf("[%s] FAIL: dim(stored)=%s vs dim(recomputed)=%s", label,
                 paste(dim(stored), collapse = "x"), paste(dim(recomputed), collapse = "x")))

  maxdiff <- max(abs(stored - recomputed))
  pass <- is.finite(maxdiff) && maxdiff < TOL

  size_train <- object.size(stored)
  size_state <- if (!is.null(fit$state.bart)) object.size(fit$state.bart) else NA

  cat(sprintf("  fit wall: %.2f s | recompute wall: %.2f s\n", t_fit, t_rec))
  cat(sprintf("  bart_train stored: %s | retained state.bart: %s\n",
              fmt_mb(size_train),
              if (is.na(size_state)) "(absent)" else fmt_mb(size_state)))
  cat(sprintf("  max |stored - recomputed| = %.3e  (tol %.0e)  ->  %s\n",
              maxdiff, TOL, if (pass) "PASS" else "FAIL"))

  invisible(list(label = label, n = n, draws = iter - warmup, chains = chains,
                 fit_s = t_fit, recompute_s = t_rec, maxdiff = maxdiff,
                 pass = pass, size_train = size_train, size_state = size_state))
}

# Scale definitions. `small` sits below the plan's n~=550 crossover; `medium`
# above it; `binary` a small probit tier for gaussian+binary coverage; `big`
# is the past-the-old-threshold blocker (n=50000 x 200 draws), sized to run in
# minutes (chains=1, warmup=50, iter=250, n.trees=50).
SCALES <- list(
  small  = function() run_scale("small (gaussian)",  n =   300L, n_levels =  8L,
                                warmup =  50L, iter = 250L, chains = 2L, n.trees = 75L),
  medium = function() run_scale("medium (gaussian)", n =  2000L, n_levels = 20L,
                                warmup =  50L, iter = 250L, chains = 2L, n.trees = 75L),
  binary = function() run_scale("small (binary)",    n =   400L, n_levels =  8L,
                                warmup =  50L, iter = 250L, chains = 2L, n.trees = 75L,
                                binary = TRUE),
  big    = function() run_scale("PAST-THRESHOLD (gaussian, C1 blocker)",
                                n = 50000L, n_levels = 40L,
                                warmup =  50L, iter = 250L, chains = 1L, n.trees = 50L)
)

.args <- commandArgs(trailingOnly = TRUE)
which <- if (length(.args) >= 1L) .args[[1L]] else "fast"
groups <- switch(which,
                 fast = c("small", "medium", "binary"),
                 big  = c("big"),
                 all  = c("small", "medium", "binary", "big"),
                 if (which %in% names(SCALES)) which else stop("unknown scale: ", which))

results <- lapply(groups, function(g) SCALES[[g]]())
ok <- all(vapply(results, function(r) isTRUE(r$pass), logical(1)))
cat(sprintf("\n==== EXACTNESS: %s (%d/%d tiers within %.0e) ====\n",
            if (ok) "ALL PASS" else "SOME FAIL",
            sum(vapply(results, function(r) isTRUE(r$pass), logical(1))),
            length(results), TOL))
if (!ok) quit(status = 1L)
