#!/usr/bin/env Rscript
# Mixing of the group standard deviation (the random-intercept sd).
#
# dbarts 1.0-0 removed rbart_vi and points users at stan4bart for multilevel
# structure. The condition attached to that removal is a bar on this package's
# group-sd chain, stated absolutely in dbarts's
# docs/design/retire-grouped-random-effects.md: lag-1 autocorrelation below 0.8
# and an effective sample size of at least 100 per 1000 kept draws. The bar was
# set on a gaussian design with n = 2000, twenty groups, a Friedman mean
# function and a true group sd of 1, where rbart_vi measured lag-1 0.141-0.150
# and stan4bart 0.967, with group-sd ESS of 20.4/154.6/636.7 against 2.1/13.2/
# 17.3 over three seeds.
#
# This script measures where stan4bart stands, on four designs covering the
# shapes the dbarts tau review used: the bar's own reference design, many small
# groups, few large groups, and a weak group signal (the corner where that
# review found mixing worst). Three seeds each. The BART-side residual sd and
# the one linear fixed effect are measured alongside, so the reader can see
# which parameter is the hard one.
#
# Usage:
#   Rscript benchmarks/mixing-group-sd.R [outfile.rds] [case ...]
#
# With no arguments every case runs and the table goes to stdout. Naming cases
# restricts the run; naming an outfile saves the per-draw summaries as well.
# Runs against the INSTALLED stan4bart; install first.

suppressMessages(library(stan4bart))

# ---- diagnostics (hand-rolled; the package has none and this adds no
# ---- dependency) ------------------------------------------------------------

## Lag-1 autocorrelation of a single chain.
lag1_acf <- function(x) {
  x <- x - mean(x)
  n <- length(x)
  sum(x[-n] * x[-1L]) / sum(x * x)
}

## Autocovariance of one chain at lags 0..n-1, by FFT (Wiener-Khinchin), with
## the 1/n normalization Stan's estimator assumes.
acov <- function(x) {
  n <- length(x)
  n_pad <- 2^ceiling(log2(2 * n))
  xc <- c(x - mean(x), rep(0, n_pad - n))
  Re(fft(Mod(fft(xc))^2, inverse = TRUE))[seq_len(n)] / (n_pad * n)
}

## Halve every chain, so that a chain drifting within itself counts against
## the diagnostic instead of hiding inside it. Both the ESS and the R-hat
## below are the split forms, which is what Stan and posterior report.
split_chains <- function(draws) {
  n <- nrow(draws)
  half <- n %/% 2L
  if (half < 2L) return(draws)
  cbind(draws[seq_len(half), , drop = FALSE],
        draws[seq.int(n - half + 1L, n), , drop = FALSE])
}

## Effective sample size over all chains of `draws` ([iterations, chains]),
## by the multi-chain autocorrelation of Gelman et al. combined with Geyer's
## initial positive sequence truncation - the same construction Stan and
## posterior use, written out here rather than depended upon.
##
## rho_t is formed from the between/within decomposition so that chains stuck
## in different places are penalized, not just slow within-chain movement; the
## initial positive sequence then truncates the sum at the first negative pair,
## which keeps the estimator from accumulating noise in the tail lags.
ess <- function(draws) {
  draws <- split_chains(as.matrix(draws))
  n <- nrow(draws); m <- ncol(draws)
  if (n < 8L) return(NA_real_)

  acovs <- vapply(seq_len(m), function(cc) acov(draws[, cc]), numeric(n))
  mean_acov <- rowMeans(acovs)
  # W, the unbiased within-chain variance; var_plus, the usual overestimate
  # that mixes in the disagreement between chain means, so chains sitting in
  # different places are penalized and not merely slow ones.
  W <- mean_acov[1L] * n / (n - 1)
  if (!is.finite(W) || W <= 0) return(NA_real_)
  var_plus <- mean_acov[1L] + if (m > 1L) var(colMeans(draws)) else 0

  # Geyer's initial positive sequence, taken over the pairs (rho_0, rho_1),
  # (rho_2, rho_3), ... and truncated at the first pair that sums to a
  # negative number; pairing from lag ZERO is what makes the sum come out at
  # 1 + 2 * sum_{t >= 1} rho_t, and pairing from lag one understates the
  # autocorrelation time by 2.
  rho_of <- function(lag) 1 - (W - mean_acov[lag + 1L]) / var_plus
  rho <- numeric(n)
  rho[1L] <- 1
  rho[2L] <- rho_of(1L)
  rho_even <- 1; rho_odd <- rho[2L]
  t <- 0L; max_t <- 0L
  while (t < n - 5L && is.finite(rho_even + rho_odd) && rho_even + rho_odd > 0) {
    t <- t + 2L
    rho_even <- rho_of(t)
    rho_odd <- rho_of(t + 1L)
    if (rho_even + rho_odd >= 0) {
      rho[t + 1L] <- rho_even
      rho[t + 2L] <- rho_odd
    }
    max_t <- t
  }
  if (rho_even > 0) rho[max_t + 1L] <- rho_even

  # Geyer's initial monotone sequence: the pair sums are non-increasing in
  # theory, so flattening any rise removes estimator noise from the tail lags.
  t <- 0L
  while (t <= max_t - 4L) {
    t <- t + 2L
    if (rho[t + 1L] + rho[t + 2L] > rho[t - 1L] + rho[t]) {
      rho[t + 1L] <- (rho[t - 1L] + rho[t]) / 2
      rho[t + 2L] <- rho[t + 1L]
    }
  }

  n_total <- n * m
  tau_hat <- -1 + 2 * sum(rho[seq_len(max_t + 1L)]) + rho[max_t + 2L]
  tau_hat <- max(tau_hat, 1 / log10(n_total))
  n_total / tau_hat
}

## Split R-hat: each chain halved before the between/within comparison, so a
## chain that drifts within itself is caught.
rhat <- function(draws) {
  split <- split_chains(as.matrix(draws))
  n <- nrow(split); m <- ncol(split)
  if (n < 2L) return(NA_real_)
  chain_means <- colMeans(split)
  W <- mean(apply(split, 2L, var))
  if (W <= 0) return(NA_real_)
  B <- n * var(chain_means)
  sqrt((((n - 1) * W + B) / n) / W)
}

## The three numbers the bar is stated in, plus the per-chain lag-1 values it
## is stated on. `draws` is [iterations, chains].
diagnose <- function(draws) {
  draws <- as.matrix(draws)
  per_chain_lag1 <- apply(draws, 2L, lag1_acf)
  n_draws <- nrow(draws) * ncol(draws)
  e <- ess(draws)
  list(lag1_mean = mean(per_chain_lag1),
       lag1_max = max(per_chain_lag1),
       lag1_per_chain = per_chain_lag1,
       ess = e,
       ess_per_1000 = 1000 * e / n_draws,
       rhat = rhat(draws),
       post_mean = mean(draws),
       n_draws = n_draws)
}

# ---- data -------------------------------------------------------------------

## The Friedman mean function, as in benchmarks/R/record-posterior-baselines.R
## and inst/common/friedmanData.R: the nonlinear part BART has to fit, so the
## group intercepts are competing with a real forest rather than with noise.
friedman_f <- function(X)
  10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]

## y = f(X) + 2 * Xfix + b_{g(i)} + eps, b_j ~ N(0, group_sd^2), eps ~ N(0, 1).
## Groups are balanced by construction (n divided evenly among n_groups), which
## keeps the design's group sizes exactly what the case name claims.
simulate <- function(seed, n, n_groups, group_sd) {
  set.seed(seed)
  X <- matrix(runif(n * 5), n, 5, dimnames = list(NULL, paste0("X", 1:5)))
  Xfix <- runif(n, -1, 1)
  g <- factor(rep_len(seq_len(n_groups), n))
  b <- rnorm(n_groups, 0, group_sd)
  y <- friedman_f(X) + 2 * Xfix + b[g] + rnorm(n)
  data.frame(X, Xfix = Xfix, g = g, y = y)
}

# ---- cases ------------------------------------------------------------------

## `bar_reference` is the design the dbarts bar was set on. The other three are
## the shapes its tau review varied over: group count at fixed n, and the
## weak-signal corner (group sd small against a residual sd of 1) where that
## review measured the worst mixing.
CASES <- list(
  bar_reference = list(n = 2000L, n_groups =  20L, group_sd = 1.0,
                       note = "20 groups of 100, group sd 1"),
  many_small    = list(n = 2000L, n_groups =  50L, group_sd = 1.0,
                       note = "50 groups of 40, group sd 1"),
  few_large     = list(n = 2000L, n_groups =   5L, group_sd = 1.0,
                       note = "5 groups of 400, group sd 1"),
  weak_signal   = list(n = 2000L, n_groups =  20L, group_sd = 0.2,
                       note = "20 groups of 100, group sd 0.2")
)

SEEDS <- c(20260913L, 20260914L, 20260915L)

## stan4bart's own defaults: chains = 4, iter = 2000 with warmup = iter %/% 2,
## so 1000 warmup and 1000 kept per chain, 4000 kept draws in all; bart_args
## and stan_args left alone (n.trees is dbarts's default 75). `cores` is the
## one argument set away from its default (getOption("mc.cores", 1)) - it
## divides the same chains across cores and changes wall time only.
CHAINS <- 4L
ITER <- 2000L
CORES <- 4L

# ---- one fit ----------------------------------------------------------------

## The formula is written literally here rather than built and passed as an
## object: stan4bart detects the bart() and lme4-bar terms by walking the
## UNEVALUATED call, so a formula held in a variable silently fails to strip
## the bart() term and model.frame() then tries to call dbarts::bart() for real.
fit_case <- function(data, mcmc_seed) {
  stan4bart::stan4bart(y ~ bart(X1 + X2 + X3 + X4 + X5) + Xfix + (1 | g), data,
                       cores = CORES, chains = CHAINS, iter = ITER,
                       verbose = -1L, seed = mcmc_seed)
}

## The three monitored scalars. The group sd is the square root of the single
## diagonal entry of the (1 | g) block's covariance - for a random intercept
## the block is 1x1, so this is the whole of it, and it is rbart_vi's tau.
monitored_draws <- function(fit) {
  Sigma <- extract(fit, "Sigma", combine_chains = FALSE)[[1L]]
  list(group_sd = sqrt(Sigma[1L, 1L, , , drop = TRUE]),
       sigma = extract(fit, "sigma", combine_chains = FALSE),
       fixef_Xfix = extract(fit, "fixef", combine_chains = FALSE)["Xfix", , ])
}

run_one <- function(case_name, seed_index) {
  cfg <- CASES[[case_name]]
  seed <- SEEDS[seed_index]
  data <- simulate(seed, cfg$n, cfg$n_groups, cfg$group_sd)

  t0 <- proc.time()[["elapsed"]]
  fit <- fit_case(data, seed + 1L)
  elapsed <- proc.time()[["elapsed"]] - t0

  draws <- monitored_draws(fit)
  rows <- lapply(names(draws), function(par) {
    d <- diagnose(draws[[par]])
    data.frame(case = case_name, seed = seed, parameter = par,
               lag1 = d$lag1_mean, lag1_max = d$lag1_max,
               ess = d$ess, ess_per_1000 = d$ess_per_1000, rhat = d$rhat,
               post_mean = d$post_mean, n_draws = d$n_draws,
               elapsed = elapsed, stringsAsFactors = FALSE)
  })
  list(summary = do.call(rbind, rows),
       lag1_per_chain = lapply(draws, function(x) apply(as.matrix(x), 2L, lag1_acf)))
}

# ---- the bar ----------------------------------------------------------------

LAG1_BAR <- 0.8
ESS_BAR <- 100

verdict <- function(row) {
  if (is.na(row$lag1) || is.na(row$ess_per_1000)) return("unknown")
  if (row$lag1 < LAG1_BAR && row$ess_per_1000 >= ESS_BAR) "pass" else "FAIL"
}

# ---- main -------------------------------------------------------------------

main <- function(args) {
  outfile <- NULL
  if (length(args) > 0L && grepl("\\.rds$", args[1L])) {
    outfile <- args[1L]
    args <- args[-1L]
  }
  case_names <- if (length(args) > 0L) args else names(CASES)
  unknown <- setdiff(case_names, names(CASES))
  if (length(unknown) > 0L)
    stop("unknown case(s): ", paste(unknown, collapse = ", "))

  results <- list()
  for (case_name in case_names) {
    for (seed_index in seq_along(SEEDS)) {
      cat(sprintf("fitting %s, seed %d ...\n", case_name, SEEDS[seed_index]))
      flush(stdout())
      results[[length(results) + 1L]] <- run_one(case_name, seed_index)
    }
  }
  summaries <- do.call(rbind, lapply(results, `[[`, "summary"))

  cat("\n")
  cat(sprintf("stan4bart defaults: %d chains, %d warmup / %d kept per chain, %d kept draws\n",
              CHAINS, ITER %/% 2L, ITER %/% 2L, CHAINS * (ITER %/% 2L)))
  cat(sprintf("bar: lag-1 below %.1f AND ESS per 1000 draws at least %d\n\n",
              LAG1_BAR, ESS_BAR))
  for (i in seq_len(nrow(summaries))) {
    row <- summaries[i, ]
    cat(sprintf("%-14s %9d %-11s lag1 %6.3f (max %6.3f)  ESS/1000 %8.1f  Rhat %6.3f  mean %8.3f  %s\n",
                row$case, row$seed, row$parameter, row$lag1, row$lag1_max,
                row$ess_per_1000, row$rhat, row$post_mean,
                if (row$parameter == "group_sd") verdict(row) else ""))
  }

  cat("\nper case, group sd only:\n")
  gs <- summaries[summaries$parameter == "group_sd", ]
  for (case_name in unique(gs$case)) {
    sub <- gs[gs$case == case_name, ]
    n_pass <- sum(vapply(seq_len(nrow(sub)), function(i) verdict(sub[i, ]) == "pass", logical(1L)))
    cat(sprintf("  %-14s %d/%d seeds pass; lag1 %.3f-%.3f, ESS/1000 %.1f-%.1f\n",
                case_name, n_pass, nrow(sub), min(sub$lag1), max(sub$lag1),
                min(sub$ess_per_1000), max(sub$ess_per_1000)))
  }
  cat(sprintf("\ntotal sampling wall time: %.1f s\n",
              sum(summaries$elapsed[summaries$parameter == "group_sd"])))

  if (!is.null(outfile)) {
    saveRDS(list(summaries = summaries, results = results, cases = CASES,
                 seeds = SEEDS, chains = CHAINS, iter = ITER),
            outfile)
    cat(sprintf("wrote %s\n", outfile))
  }
  invisible(summaries)
}

## Run only when this file is the script Rscript was pointed at, so that
## sourcing it for the diagnostics alone does not start a fit (the guard
## benchmarks/R/record-posterior-baselines.R uses).
.self_file <- sub("--file=", "",
                  grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
if (!is.na(.self_file) && basename(.self_file) == "mixing-group-sd.R")
  main(commandArgs(trailingOnly = TRUE))
