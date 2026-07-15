#!/usr/bin/env Rscript
# C0(b) POSTERIOR baselines (docs/plans/walnuts.md): fixef/ranef/sigma/Sigma
# posterior means, sds, and MCSEs, recorded from the Stan-era sampler on five
# reference datasets spanning the WALNUTS gradient's tiers, while Stan is
# still the only sampler. These are the oracle the C2/C3 distributional-
# equivalence gates compare a future WALNUTS fit against.
#
# Usage:
#   Rscript benchmarks/R/record-posterior-baselines.R record [outfile.rds]
#   Rscript benchmarks/R/record-posterior-baselines.R rerun-check <tier> <new_seed> [baseline.rds]
#
# "record" regenerates all five reference datasets from their fixed seeds
# (see TIERS below) and fits each with stan4bart, deterministically. "rerun-
# check" is the C0 stability gate: refits one tier's SAME data at a
# different MCMC seed and reports whether every gated summary's mean shift
# is within ~2x the combined MCSE of the two runs (the plan's abort
# criterion: "if a reference summary is not stable across reruns ... lengthen
# before freezing").
#
# PREREQUISITE - THE structSize FIX (C0(c)): these baselines are only valid
# against a build carrying the bart_util.hpp fix for IterableBartResults:
#   dbarts_results current = {};                             // VD's one-liner
#   current.structSize = sizeof(dbarts_results);             // in the ctor -
#                                                            // ALSO REQUIRED
# Without it, `current.structSize` is uninitialized (in practice 0 on fresh
# heap), dbarts_sampler_run's DBARTS_RESULTS_HAS skips EVERY output field,
# the train-fit buffers stay zero, and init.cpp then feeds Stan
# stanOffset = 0 - parametricMean: the Gibbs coupling decoheres (beta
# random-walks hundreds of prior sds out, stored "BART" fits are exactly
# -X*beta, predict-vs-extract testthat tests fail, and chains progressively
# lock into a permanently-divergent zero-acceptance state as the runaway
# target outruns the frozen post-warmup stepsize). Diagnosed at C0 recording
# time; VD's uncommitted one-liner alone (the `= {}`) makes the skip
# deterministic but does NOT restore output - the structSize assignment is
# the load-bearing half. On the fixed build all chains mix cleanly at
# DEFAULT stan_args across every tier; chain_health() below is retained as a
# cheap guard (record aborts if chains die), not as a survival filter.

suppressMessages(library(stan4bart))

# ---- shared simulation pieces ----------------------------------------------

## Friedman-style nonlinear signal for the BART covariates (X1..X5); the same
## functional form as inst/common/friedmanData.R, so BART's side is
## well-identified and the setup matches the package's own test conventions.
friedman_f <- function(X)
  10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]

## One continuous-response reference dataset. `nc` sets the random-effect
## block dimension (1, 2, or 3+ - the decov "onion" tier); `weighted` adds an
## observation-weight column exercising the weighted-likelihood log-
## normalizer (continuous.stan:365's N-not-sum(w) quirk). Fully seeded, so
## the same `seed` always regenerates the same data frame (no external data
## files - regeneration is the reproducibility contract for these
## baselines).
simulate_continuous <- function(seed, n, nc, n_levels, weighted = FALSE) {
  set.seed(seed)

  X <- matrix(runif(n * 5), n, 5, dimnames = list(NULL, paste0("X", 1:5)))
  Xfix <- runif(n, -1, 1)
  Xr1 <- if (nc >= 2) runif(n, -1, 1) else NULL
  Xr2 <- if (nc >= 3) runif(n, -1, 1) else NULL
  g <- factor(sample.int(n_levels, n, replace = TRUE))

  if (nc == 1) {
    true_sd <- 1.2
    true_corr <- NULL
    b <- rnorm(n_levels, 0, true_sd)
    ranef <- b[g]
  } else if (nc == 2) {
    true_corr <- matrix(c(1, 0.4, 0.4, 1), 2)
    true_sd <- c(1.5, 1.0)
    Sigma_b <- diag(true_sd) %*% true_corr %*% diag(true_sd)
    b <- matrix(rnorm(n_levels * 2), n_levels) %*% chol(Sigma_b)
    ranef <- b[g, 1] + Xr1 * b[g, 2]
  } else {
    true_corr <- matrix(c(1, 0.35, -0.25, 0.35, 1, 0.3, -0.25, 0.3, 1), 3)
    true_sd <- c(1.4, 1.1, 0.9)
    Sigma_b <- diag(true_sd) %*% true_corr %*% diag(true_sd)
    b <- matrix(rnorm(n_levels * 3), n_levels) %*% chol(Sigma_b)
    ranef <- b[g, 1] + Xr1 * b[g, 2] + Xr2 * b[g, 3]
  }

  mu <- friedman_f(X) + 2 * Xfix + ranef

  weights <- if (weighted) runif(n, 0.4, 2.5) else NULL
  noise_sd <- if (weighted) 1 / sqrt(weights) else 1
  y <- mu + rnorm(n, 0, noise_sd)

  df <- data.frame(X, Xfix = Xfix, g = g, y = y)
  if (!is.null(Xr1)) df$Xr1 <- Xr1
  if (!is.null(Xr2)) df$Xr2 <- Xr2
  if (!is.null(weights)) df$w <- weights

  # NOTE: the formula is NOT stored as an object here (see fit_tier()'s
  # docstring) - stan4bart's bart()/lme4-bar term detection walks the
  # UNEVALUATED call captured by match.call(), so it only fires when the
  # formula is written literally at the stan4bart() call site; a formula
  # held in a variable (this list's $shape is a proxy for it) silently
  # fails to strip the bart() term, and model.frame() then tries to call
  # dbarts::bart() for real and errors on a missing y.train. `shape`
  # records which of fit_tier()'s literal call sites to use.
  shape <- paste0("nc", nc)

  list(df = df, shape = shape, weights_col = if (weighted) "w" else NULL,
       is_binary = FALSE,
       sim = list(seed = seed, n = n, nc = nc, n_levels = n_levels,
                  weighted = weighted, true_sd = true_sd, true_corr = true_corr))
}

## The binary reference dataset: probit response over TWO random-effect
## blocks - an nc=2 correlated intercept+slope block (g) plus an nc=1
## random-intercept block (g2) - mirroring inst/common/friedmanData.R's
## g.1/g.2 pair and exercising multi-block geometry, which the continuous
## tiers (single block each) do not. The latent scaling via
## sd(mu)/qnorm(0.15) follows friedmanData.R's convention.
simulate_binary <- function(seed, n, n_levels, n_levels2) {
  set.seed(seed)

  X <- matrix(runif(n * 5), n, 5, dimnames = list(NULL, paste0("X", 1:5)))
  Xfix <- runif(n, -1, 1)
  Xr1 <- runif(n, -1, 1)
  g <- factor(sample.int(n_levels, n, replace = TRUE))
  g2 <- factor(sample.int(n_levels2, n, replace = TRUE))

  true_corr <- matrix(c(1, 0.4, 0.4, 1), 2)
  true_sd <- c(1.5, 1.0)
  Sigma_b <- diag(true_sd) %*% true_corr %*% diag(true_sd)
  b <- matrix(rnorm(n_levels * 2), n_levels) %*% chol(Sigma_b)
  true_sd2 <- 1.0
  b2 <- rnorm(n_levels2, 0, true_sd2)
  ranef <- b[g, 1] + Xr1 * b[g, 2] + b2[g2]

  mu <- friedman_f(X) + 2 * Xfix + ranef
  loc <- mean(mu); scale <- sd(mu) / qnorm(0.15)
  mu <- (mu - loc) / scale
  y <- rbinom(n, 1L, pnorm(mu))

  df <- data.frame(X, Xfix = Xfix, Xr1 = Xr1, g = g, g2 = g2, y = y)

  list(df = df, shape = "binary", weights_col = NULL, is_binary = TRUE,
       sim = list(seed = seed, n = n, n_levels = n_levels, n_levels2 = n_levels2,
                  true_sd = true_sd, true_sd2 = true_sd2, true_corr = true_corr))
}

# ---- reference tiers --------------------------------------------------------

TIERS <- list(
  continuous_nc1      = list(kind = "continuous", seed = 20260715L, n = 350L, nc = 1L, n_levels = 20L),
  continuous_nc2      = list(kind = "continuous", seed = 20260716L, n = 400L, nc = 2L, n_levels = 20L),
  continuous_nc3      = list(kind = "continuous", seed = 20260717L, n = 500L, nc = 3L, n_levels = 30L),
  weighted_continuous = list(kind = "continuous", seed = 20260718L, n = 400L, nc = 2L, n_levels = 20L, weighted = TRUE),
  binary_multilevel   = list(kind = "binary",     seed = 20260719L, n = 400L, n_levels = 20L, n_levels2 = 10L)
)

build_tier_data <- function(cfg) {
  if (cfg$kind == "continuous") {
    simulate_continuous(cfg$seed, cfg$n, cfg$nc, cfg$n_levels, isTRUE(cfg$weighted))
  } else {
    simulate_binary(cfg$seed, cfg$n, cfg$n_levels, cfg$n_levels2)
  }
}

# ---- MCMC control -----------------------------------------------------------

## One control for every tier: DEFAULT stan_args (on the structSize-fixed
## build all tiers mix cleanly without tuning - see the module docstring),
## 8 chains x 3000 kept draws = 24000 retained draws per summary, long
## enough that the batch-means MCSEs make the C2/C3 k=4 tolerance sharp.
CONTROL <- list(chains = 8L, cores = 8L, warmup = 1000L, iter = 4000L,
                n.trees = 50L, stan_args = NULL)

control_for <- function(cfg) CONTROL

## Dispatches to one of four LITERAL stan4bart() call sites, keyed by
## tier_data$shape. This duplication is required, not stylistic: see the
## comment in simulate_continuous() - stan4bart's bart()/lme4-bar term
## detection needs the formula written out at the call site, so a single
## generic call built from a stored formula object silently mis-fits.
## `weights` is always passed (possibly NULL, which stan4bart treats as
## "no weights") so the weighted and unweighted paths share one branch.
fit_tier <- function(tier_data, control, mcmc_seed) {
  data <- tier_data$df
  weights_vec <- if (!is.null(tier_data$weights_col)) data[[tier_data$weights_col]] else NULL

  if (tier_data$shape == "nc1") {
    stan4bart::stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 | g), data,
                         weights = weights_vec,
                         cores = control$cores, verbose = -1L,
                         chains = control$chains, warmup = control$warmup,
                         iter = control$iter, seed = mcmc_seed,
                         bart_args = list(n.trees = control$n.trees),
                         stan_args = control$stan_args, keep_fits = TRUE)
  } else if (tier_data$shape == "nc2") {
    stan4bart::stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 + Xr1 | g), data,
                         weights = weights_vec,
                         cores = control$cores, verbose = -1L,
                         chains = control$chains, warmup = control$warmup,
                         iter = control$iter, seed = mcmc_seed,
                         bart_args = list(n.trees = control$n.trees),
                         stan_args = control$stan_args, keep_fits = TRUE)
  } else if (tier_data$shape == "nc3") {
    stan4bart::stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 + Xr1 + Xr2 | g), data,
                         weights = weights_vec,
                         cores = control$cores, verbose = -1L,
                         chains = control$chains, warmup = control$warmup,
                         iter = control$iter, seed = mcmc_seed,
                         bart_args = list(n.trees = control$n.trees),
                         stan_args = control$stan_args, keep_fits = TRUE)
  } else if (tier_data$shape == "binary") {
    stan4bart::stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 + Xr1 | g) + (1 | g2), data,
                         cores = control$cores, verbose = -1L,
                         chains = control$chains, warmup = control$warmup,
                         iter = control$iter, seed = mcmc_seed,
                         bart_args = list(n.trees = control$n.trees),
                         stan_args = control$stan_args, keep_fits = TRUE)
  } else {
    stop("unknown tier shape: ", tier_data$shape)
  }
}

# ---- chain health guard ------------------------------------------------------

## Flags chains whose HMC transitions stopped mixing in their trailing
## window. On a broken build (missing the structSize fix - module docstring)
## chains progressively lock into a permanently-divergent zero-acceptance
## state as the decohered target runs away; on the fixed build every chain
## passes, so this is a cheap tripwire, not a survival filter. A chain is
## unhealthy only if its trailing 20% shows BOTH collapsed acceptance AND
## pervasive divergences - a momentarily-low accept_stat with no divergences
## is ordinary NUTS noise, not this failure.
chain_health <- function(fit, tail_frac = 0.2, accept_floor = 0.3, divergent_ceiling = 0.5) {
  nms <- dimnames(fit$stan)[[1L]]
  accept_idx <- match("accept_stat__", nms)
  div_idx <- match("divergent__", nms)
  n_iter <- dim(fit$stan)[2L]
  n_chains <- dim(fit$stan)[3L]
  tail_idx <- seq.int(max(1L, floor(n_iter * (1 - tail_frac)) + 1L), n_iter)
  vapply(seq_len(n_chains), function(cc) {
    accept_tail <- fit$stan[accept_idx, tail_idx, cc]
    div_tail <- fit$stan[div_idx, tail_idx, cc]
    mean(accept_tail) >= accept_floor && mean(div_tail) <= divergent_ceiling
  }, logical(1L))
}

# ---- batch-means MCSE --------------------------------------------------------

## Posterior mean, sd, and MCSE of `draws` (an [iterations, chains] matrix,
## one column per KEPT chain). Two hand-rolled estimators (no coda
## dependency), reported separately, with the HONEST mcse = their max:
##   mcse_bm: non-overlapping batch means, batch b = n^(2/3) per chain
##     (Flegal-Jones sizing; sqrt(n) batches proved far too short here),
##     pooled across chains (independent chains -> independent batches).
##   mcse_bc: sd of the per-chain means / sqrt(n_chains). Chains are
##     independent, so this captures ALL within-chain autocorrelation by
##     construction - including the slow BART-vs-parametric tradeoff that
##     inflates autocorrelation far past any within-chain batch (observed at
##     C0 recording: mcse_bc 3-8x the sqrt(n)-batch estimate on the ranef
##     summaries). Only n_chains - 1 df, hence the max() rather than using
##     it alone.
## The C0 rerun gate and the C2/C3 tolerance both consume `mcse`.
batch_means_summary <- function(draws) {
  draws <- as.matrix(draws)
  n <- nrow(draws); n_chains <- ncol(draws)

  mcse_bc <- if (n_chains >= 2L) sd(colMeans(draws)) / sqrt(n_chains) else NA_real_

  b <- max(2L, floor(n^(2 / 3)))
  m <- floor(n / b)
  if (m >= 2L) {
    batch_means <- vapply(seq_len(n_chains), function(cc) {
      x <- draws[seq_len(m * b), cc]
      colMeans(matrix(x, nrow = b))
    }, numeric(m))
    mcse_bm <- sd(as.vector(batch_means)) / sqrt(m * n_chains)
  } else {
    # too few draws per chain to batch - the naive (autocorrelation-blind)
    # standard error, better than nothing when mcse_bc is also unavailable.
    mcse_bm <- sd(as.vector(draws)) / sqrt(length(draws))
    m <- NA_integer_
  }

  list(mean = mean(draws), sd = sd(as.vector(draws)),
       mcse = max(mcse_bm, mcse_bc, na.rm = TRUE),
       mcse_bm = mcse_bm, mcse_bc = mcse_bc,
       n_draws = n * n_chains,
       n_batches = if (is.na(m)) NA_integer_ else m * n_chains)
}

# ---- posterior summaries: fixef/ranef/sigma/Sigma (the C2/C3 gated set) ----

## Every summary the plan's C2/C3 gate gives a tolerance for: fixef, ranef
## (b), sigma (continuous only), and Sigma (reduced to per-block sds and
## pairwise correlations, the natural scale-free summary of a covariance
## draw - see docs/plans/walnuts.md C2's "for every fixef/ranef/sigma/Sigma
## summary"). `keep_chains` is the logical vector from chain_health().
summarize_fit <- function(fit, is_binary, keep_chains) {
  rows <- list()
  add_row <- function(quantity, name, s)
    rows[[length(rows) + 1L]] <<- data.frame(
      quantity = quantity, name = name,
      mean = s$mean, sd = s$sd, mcse = s$mcse,
      mcse_bm = s$mcse_bm, mcse_bc = s$mcse_bc, n_draws = s$n_draws,
      stringsAsFactors = FALSE)

  fixef <- extract(fit, "fixef", combine_chains = FALSE)[, , keep_chains, drop = FALSE]
  for (p in dimnames(fixef)$predictor)
    add_row("fixef", p, batch_means_summary(fixef[p, , ]))

  ranef_full <- extract(fit, "ranef", combine_chains = FALSE)
  for (gname in names(ranef_full)) {
    arr <- ranef_full[[gname]][, , , keep_chains, drop = FALSE]
    for (p in dimnames(arr)$predictor)
      for (grp in dimnames(arr)$group)
        add_row("ranef", sprintf("%s:%s:%s", gname, p, grp),
                batch_means_summary(arr[p, grp, , ]))
  }

  if (!is_binary) {
    sigma <- extract(fit, "sigma", combine_chains = FALSE)[, keep_chains, drop = FALSE]
    add_row("sigma", "sigma", batch_means_summary(sigma))
  }

  Sigma_full <- extract(fit, "Sigma", combine_chains = FALSE)
  for (gname in names(Sigma_full)) {
    arr <- Sigma_full[[gname]][, , , keep_chains, drop = FALSE]
    p <- dim(arr)[1L]; nms <- dimnames(arr)[[1L]]
    n_iter <- dim(arr)[3L]; n_kept <- dim(arr)[4L]
    sd_draws <- array(NA_real_, c(p, n_iter, n_kept))
    for (i in seq_len(p)) sd_draws[i, , ] <- sqrt(arr[i, i, , ])
    for (i in seq_len(p))
      add_row("Sigma_sd", sprintf("%s:%s", gname, nms[i]),
              batch_means_summary(sd_draws[i, , ]))
    if (p >= 2L) {
      for (i in seq_len(p - 1L)) for (j in seq.int(i + 1L, p)) {
        corr_draws <- arr[i, j, , ] / (sd_draws[i, , ] * sd_draws[j, , ])
        add_row("Sigma_corr", sprintf("%s:%s,%s", gname, nms[i], nms[j]),
                batch_means_summary(corr_draws))
      }
    }
  }

  do.call(rbind, rows)
}

# ---- one tier: simulate, fit, filter, summarize ----------------------------

record_one_tier <- function(tier_name, cfg, mcmc_seed = cfg$seed) {
  tier_data <- build_tier_data(cfg)
  control <- control_for(cfg)
  fit <- fit_tier(tier_data, control, mcmc_seed)
  healthy <- chain_health(fit)
  n_healthy <- sum(healthy)
  cat(sprintf("[%s] mcmc_seed=%d healthy chains: %d/%d\n",
              tier_name, mcmc_seed, n_healthy, length(healthy)))
  # On the structSize-fixed build every chain mixes (see module docstring);
  # any unhealthy chain now warrants a look, and more than one aborts.
  if (n_healthy < length(healthy) - 1L)
    stop(sprintf("[%s] only %d/%d healthy chains at mcmc_seed=%d - the build is likely missing the structSize fix; investigate before freezing",
                 tier_name, n_healthy, length(healthy), mcmc_seed))
  summary_df <- summarize_fit(fit, is_binary = identical(cfg$kind, "binary"), keep_chains = healthy)
  list(tier = tier_name, cfg = cfg, control = control, mcmc_seed = mcmc_seed,
       n_chains = length(healthy), n_healthy = n_healthy, healthy = healthy,
       data = tier_data$df, shape = tier_data$shape,
       weights_col = tier_data$weights_col, sim = tier_data$sim,
       summary = summary_df)
}

# ---- record mode: all five tiers -> one combined .rds ----------------------

## Saves incrementally after every tier (a late-tier failure must not lose
## completed tiers); re-recording a subset merges into an existing outfile.
record_baselines <- function(outfile, tier_names = names(TIERS)) {
  existing <- if (file.exists(outfile)) readRDS(outfile)$tiers else list()
  manifest <- list(commit = "967a1c6",
                   package_state = "967a1c6 + bart_util.hpp structSize fix (C0(c); see module docstring)",
                   date = as.character(Sys.Date()),
                   package_version = as.character(utils::packageVersion("stan4bart")),
                   r_version = R.version.string,
                   control = CONTROL,
                   chain_health_params = list(tail_frac = 0.2, accept_floor = 0.3,
                                               divergent_ceiling = 0.5))
  results <- existing
  for (nm in tier_names) {
    results[[nm]] <- record_one_tier(nm, TIERS[[nm]])
    saveRDS(list(manifest = manifest, tiers = results), outfile)
    cat("saved through", nm, "\n")
  }
  cat("wrote", outfile, "\n")
  invisible(results)
}

# ---- rerun-check mode: the C0 stability gate --------------------------------

## Refits `tier_name`'s SAME reference data at a different MCMC seed and
## compares every gated summary's mean shift against k times the combined
## MCSE of the two runs. This is the plan's C0 abort criterion made concrete
## - not yet the C2/C3 equivalence gate itself (that also compares against a
## WALNUTS fit, which does not exist yet), but the same comparison machinery.
## READING THE RESULT: with ~100 gated summaries, a handful of ratios in
## (2, 3) is the EXPECTED tail of well-calibrated MCSEs (P(|Z| > 2) ~ 4.6%
## per summary), not instability; what aborts is a systematic excess or
## ratios approaching C2's k = 4. The verdict line prints the expected
## exceedance count for calibration.
rerun_stability_check <- function(tier_name, new_seed, baseline_file, k = 2) {
  baseline <- readRDS(baseline_file)
  base_tier <- baseline$tiers[[tier_name]]
  if (is.null(base_tier)) stop("tier not found in baseline: ", tier_name)

  rerun <- record_one_tier(tier_name, base_tier$cfg, mcmc_seed = new_seed)

  merged <- merge(base_tier$summary, rerun$summary, by = c("quantity", "name"),
                  suffixes = c("_base", "_rerun"))
  merged$combined_mcse <- sqrt(merged$mcse_base^2 + merged$mcse_rerun^2)
  merged$shift <- abs(merged$mean_base - merged$mean_rerun)
  merged$ratio <- merged$shift / merged$combined_mcse
  merged$pass <- merged$ratio <= k

  cat(sprintf("[%s] stability check (k=%g): %d/%d summaries within tolerance (expected exceedances if calibrated: ~%.1f); max ratio %.2f\n",
              tier_name, k, sum(merged$pass), nrow(merged),
              2 * pnorm(-k) * nrow(merged), max(merged$ratio)))
  print(merged[order(-merged$ratio),
              c("quantity", "name", "mean_base", "mean_rerun", "combined_mcse", "ratio", "pass")],
       row.names = FALSE)
  invisible(merged)
}

# ---- CLI --------------------------------------------------------------------

## Guarded on THIS file being the one Rscript was invoked on (not merely
## source()d - bench-perf.R sources this file for its simulate_*/TIERS
## definitions, and commandArgs() reflects the outer process's real argv
## regardless of which file called source(), so without this guard sourcing
## from bench-perf.R would re-trigger whatever mode the OUTER script was
## invoked with).
.self_file_rpb <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
.is_main_rpb <- !is.na(.self_file_rpb) && basename(.self_file_rpb) == "record-posterior-baselines.R"

.default_outfile <- "benchmarks/baselines/posterior-baselines-967a1c6.rds"
.args <- commandArgs(trailingOnly = TRUE)
if (.is_main_rpb && length(.args) >= 1L) {
  mode <- .args[[1L]]
  if (mode == "record") {
    outfile <- if (length(.args) >= 2L) .args[[2L]] else .default_outfile
    record_baselines(outfile)
  } else if (mode == "rerun-check") {
    tier_name <- .args[[2L]]
    new_seed <- as.integer(.args[[3L]])
    baseline_file <- if (length(.args) >= 4L) .args[[4L]] else .default_outfile
    rerun_stability_check(tier_name, new_seed, baseline_file)
  } else {
    stop("unknown mode: ", mode, " (expected 'record' or 'rerun-check')")
  }
}
