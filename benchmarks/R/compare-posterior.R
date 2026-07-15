#!/usr/bin/env Rscript
# C2 DISTRIBUTIONAL-EQUIVALENCE gate (docs/plans/walnuts.md C2). Refits each
# CONTINUOUS reference tier on the WALNUTS build using the SAME stored data,
# seed, chains, and iterations as the C0 Stan-era baseline, recomputes the same
# fixef/ranef/sigma/Sigma summaries with the same batch-means MCSE estimator,
# and gates every gated summary against the PRE-REGISTERED tolerance:
#   means:  |mean_walnuts - mean_stan| < k * sqrt(mcse_stan^2 + mcse_walnuts^2), k = 4
#   sds:    sd_walnuts / sd_stan in [0.8, 1.25]
# As of C3 the binary tier is gated too (binary routed through WALNUTS; its
# gated set has no sigma row). Reuses fit_tier / summarize_fit / chain_health /
# batch_means_summary from record-posterior-baselines.R (sourced; its CLI is
# guarded, so sourcing is inert).
#
# Usage (one tier per invocation to stay inside the foreground time budget):
#   Rscript benchmarks/R/compare-posterior.R compare <tier> [baseline.rds]
#   Rscript benchmarks/R/compare-posterior.R compare-all [baseline.rds]
# tiers: continuous_nc1 continuous_nc2 continuous_nc3 weighted_continuous
#        binary_multilevel

suppressMessages(library(stan4bart))

.here <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
.rpb <- file.path(dirname(.here), "record-posterior-baselines.R")
if (!file.exists(.rpb)) .rpb <- "benchmarks/R/record-posterior-baselines.R"
source(.rpb, local = FALSE)   # fit_tier, summarize_fit, chain_health, batch_means_summary

.default_baseline <- "benchmarks/baselines/posterior-baselines-75d7970.rds"
GATED_TIERS <- c("continuous_nc1", "continuous_nc2", "continuous_nc3",
                 "weighted_continuous", "binary_multilevel")

K_MEAN <- 4          # pre-registered mean-tolerance multiplier
SD_LO  <- 0.80       # pre-registered sd-ratio band
SD_HI  <- 1.25

# Refit one tier's stored data on the WALNUTS build and summarize, reusing the
# baseline's control + mcmc_seed so the comparison is like-for-like.
walnuts_fit_summary <- function(base_tier) {
  is_binary <- identical(base_tier$shape, "binary")
  tier_data <- list(df = base_tier$data, shape = base_tier$shape,
                    weights_col = base_tier$weights_col, is_binary = is_binary)
  fit <- fit_tier(tier_data, base_tier$control, base_tier$mcmc_seed)
  healthy <- chain_health(fit)
  cat(sprintf("[%s] WALNUTS healthy chains: %d/%d\n",
              base_tier$tier, sum(healthy), length(healthy)))
  summarize_fit(fit, is_binary = is_binary, keep_chains = healthy)
}

gate_tier <- function(tier_name, baseline_file) {
  baseline <- readRDS(baseline_file)
  base_tier <- baseline$tiers[[tier_name]]
  if (is.null(base_tier)) stop("tier not found in baseline: ", tier_name)

  wal <- walnuts_fit_summary(base_tier)

  m <- merge(base_tier$summary, wal, by = c("quantity", "name"),
             suffixes = c("_stan", "_wal"))
  m$combined_mcse <- sqrt(m$mcse_stan^2 + m$mcse_wal^2)
  m$dmean <- abs(m$mean_stan - m$mean_wal)
  m$band <- K_MEAN * m$combined_mcse
  m$mean_ratio <- m$dmean / m$combined_mcse           # pass if < K_MEAN
  m$sd_ratio <- m$sd_wal / m$sd_stan                  # pass if in [SD_LO, SD_HI]
  m$mean_pass <- m$mean_ratio < K_MEAN
  m$sd_pass <- m$sd_ratio >= SD_LO & m$sd_ratio <= SD_HI

  worst_mean <- m[which.max(m$mean_ratio), ]
  sd_dev <- pmax(m$sd_ratio, 1 / m$sd_ratio)          # distance from 1 either way
  worst_sd <- m[which.max(sd_dev), ]

  cat(sprintf("\n=== TIER %s: %d gated summaries ===\n", tier_name, nrow(m)))
  cat(sprintf("  MEAN gate (k=%g): %d/%d pass; worst mean_ratio=%.2f (%s: |dmean|=%.4g, band=%.4g)\n",
              K_MEAN, sum(m$mean_pass), nrow(m), max(m$mean_ratio),
              paste(worst_mean$quantity, worst_mean$name), worst_mean$dmean, worst_mean$band))
  cat(sprintf("  SD gate [%.2f,%.2f]: %d/%d pass; worst sd_ratio=%.3f (%s)\n",
              SD_LO, SD_HI, sum(m$sd_pass), nrow(m), worst_sd$sd_ratio,
              paste(worst_sd$quantity, worst_sd$name)))
  fails <- m[!m$mean_pass | !m$sd_pass,
             c("quantity", "name", "mean_stan", "mean_wal", "combined_mcse",
               "mean_ratio", "sd_ratio")]
  if (nrow(fails) > 0) {
    cat("  FAILURES:\n"); print(fails[order(-fails$mean_ratio), ], row.names = FALSE)
  }
  verdict <- all(m$mean_pass) && all(m$sd_pass)
  cat(sprintf("  VERDICT [%s]: %s\n", tier_name, if (verdict) "PASS" else "FAIL"))
  invisible(list(tier = tier_name, pass = verdict, merged = m))
}

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L) {
  mode <- .args[[1L]]
  if (mode == "compare") {
    tier_name <- .args[[2L]]
    baseline_file <- if (length(.args) >= 3L) .args[[3L]] else .default_baseline
    res <- gate_tier(tier_name, baseline_file)
    if (!res$pass) quit(status = 1L)
  } else if (mode == "compare-all") {
    baseline_file <- if (length(.args) >= 2L) .args[[2L]] else .default_baseline
    ok <- TRUE
    for (tn in GATED_TIERS) ok <- gate_tier(tn, baseline_file)$pass && ok
    cat(sprintf("\n==== OVERALL: %s ====\n", if (ok) "ALL TIERS PASS" else "SOME TIERS FAIL"))
    if (!ok) quit(status = 1L)
  } else {
    stop("unknown mode: ", mode, " (expected 'compare' or 'compare-all')")
  }
}
