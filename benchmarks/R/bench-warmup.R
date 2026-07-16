#!/usr/bin/env Rscript
# C0 WARMUP measurement harness (docs/plans/walnuts-warmup.md C0) - the
# instrument that judges every later warmup/seed candidate. For each tier x
# warmup it reports the three quantities the arc trades off:
#   - per-iteration SAMPLING wall (ms), by matched-warmup differencing (the
#     bench-perf.R pattern: same warmup + seed, a short and a long `iter`; the
#     slope (t_long - t_short)/(n_samp_long - n_samp_short) cancels the fixed
#     setup + warmup cost), so it isolates the cost of one sampling transition.
#   - ESS/sec: posterior::ess_bulk / ess_tail of the MONITORED scalar set (each
#     fixef, sigma, and the log-scale group SDs from Sigma - sample-storage.md's
#     slow-mixing canaries) divided by the differenced sampling wall. The honest
#     floor: a tuning change that cuts leapfrog steps but hurts mixing is a net
#     loss. min across the set is the canary; median the typical scalar.
#   - mean leapfrog / eval per transition, for the WARMUP and SAMPLING phases,
#     read from fit$adaptation$mean_leapfrog[_warmup] (the C0 eval counter). This
#     is the mechanistic quantity the per-iter wall tracks.
#
# Tiers: the five distributional reference tiers (record-posterior-baselines.R's
# TIERS) plus a new large-n-RE tier (n=1e5 continuous, q=1000 random-intercept
# levels) reproducing perf-review.md's 69->28 ms knee. The warmup grid
# {20,50,100,200,400,1000} traces that knee; on the small reference tiers it is
# cheap and uniform. Single chain, single core; keep_fits = TRUE (fit$stan
# carries the parametric draws only when samples are stored, and it matches how
# perf-review.md measured the knee).
#
# TIMING DISCIPLINE: run on an otherwise-idle machine, tiers sequentially, no
# concurrent load (the large-n-RE tier is memory- and wall-heavy). The reference
# tiers are seeded deterministically from record-posterior-baselines.R; the
# large tier from its own fixed seed.
#
# Usage:
#   Rscript benchmarks/R/bench-warmup.R record [outfile.csv] [tier1,tier2,...]
#   Rscript benchmarks/R/bench-warmup.R quick  [outfile.csv]   # tiny smoke test
# With no tier list, `record` runs all six tiers. `quick` runs a couple of
# tiers at a tiny grid/size to check the plumbing (NOT a baseline).

suppressMessages({
  library(stan4bart)
  library(posterior)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

## Sibling-relative source, like bench-perf.R (Rscript sets no working dir).
.self_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
SELF_DIR <- if (!is.na(.self_file) && nzchar(.self_file)) dirname(normalizePath(.self_file)) else "benchmarks/R"
source(file.path(SELF_DIR, "record-posterior-baselines.R"))  # simulate_*, TIERS, build_tier_data

COMMIT <- tryCatch(system2("git", c("-C", SELF_DIR, "rev-parse", "--short", "HEAD"),
                           stdout = TRUE, stderr = FALSE)[1L], error = function(e) NA_character_)

# ---- tiers -----------------------------------------------------------------

## The five reference tiers (reused verbatim) plus the large-n-RE perf tier.
## Each entry yields a tier_data list (df, shape, weights_col, is_binary) via
## build_one() below.
WARMUP_TIERS <- c(
  as.list(setNames(names(TIERS), names(TIERS))),  # placeholders; resolved in build_one
  list(large_n_re = "large_n_re"))

## Fixed seed for the large-n-RE tier (continuous, nc=1, q=1000 intercepts).
LARGE_N_RE_SIM <- list(seed = 20260720L, n = 100000L, nc = 1L, n_levels = 1000L)

build_one <- function(tier_name) {
  if (tier_name == "large_n_re") {
    with(LARGE_N_RE_SIM, simulate_continuous(seed, n, nc, n_levels, weighted = FALSE))
  } else {
    build_tier_data(TIERS[[tier_name]])
  }
}

## n and q (RE levels) for the CSV, read off the tier config.
tier_n_q <- function(tier_name) {
  if (tier_name == "large_n_re") return(c(n = LARGE_N_RE_SIM$n, q = LARGE_N_RE_SIM$n_levels))
  cfg <- TIERS[[tier_name]]
  q <- cfg$n_levels %||% NA_integer_
  c(n = cfg$n, q = q)
}

# ---- one fit at a given warmup/iter -----------------------------------------

## LITERAL-formula dispatch by shape (the bart()/bar-term detection walks the
## unevaluated call - see record-posterior-baselines.R:simulate_continuous - so
## each formula must be written out, not passed in a variable). Single chain,
## single core, keep_fits = TRUE (fit$stan is populated only when samples store).
fit_at <- function(tier_data, warmup, iter, seed, n.trees = 50L) {
  data <- tier_data$df
  weights_vec <- if (!is.null(tier_data$weights_col)) data[[tier_data$weights_col]] else NULL
  ba <- list(n.trees = n.trees)
  if (tier_data$shape == "nc1") {
    stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 | g), data, weights = weights_vec,
              cores = 1L, chains = 1L, verbose = -1L, warmup = warmup, iter = iter,
              seed = seed, bart_args = ba, keep_fits = TRUE)
  } else if (tier_data$shape == "nc2") {
    stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 | g), data, weights = weights_vec,
              cores = 1L, chains = 1L, verbose = -1L, warmup = warmup, iter = iter,
              seed = seed, bart_args = ba, keep_fits = TRUE)
  } else if (tier_data$shape == "nc3") {
    stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 + Xr2 | g), data, weights = weights_vec,
              cores = 1L, chains = 1L, verbose = -1L, warmup = warmup, iter = iter,
              seed = seed, bart_args = ba, keep_fits = TRUE)
  } else if (tier_data$shape == "binary") {
    stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 | g) + (1 | g2), data,
              cores = 1L, chains = 1L, verbose = -1L, warmup = warmup, iter = iter,
              seed = seed, bart_args = ba, keep_fits = TRUE)
  } else stop("unknown shape: ", tier_data$shape)
}

# ---- monitored scalars + ESS ------------------------------------------------

## The monitored set as single-chain draw vectors: each fixef, sigma (continuous),
## and every group's log-scale SDs from Sigma (the boundary-variance canaries).
monitored_draws <- function(fit, is_binary) {
  out <- list()
  fixef <- extract(fit, "fixef", combine_chains = FALSE)          # [predictor, iter, chain]
  for (p in dimnames(fixef)$predictor) out[[paste0("fixef:", p)]] <- as.numeric(fixef[p, , 1L])
  if (!is_binary) {
    sigma <- extract(fit, "sigma", combine_chains = FALSE)        # [iter, chain]
    out[["sigma"]] <- as.numeric(sigma[, 1L])
  }
  Sigma <- extract(fit, "Sigma", combine_chains = FALSE)          # list of [p,p,iter,chain]
  for (gname in names(Sigma)) {
    arr <- Sigma[[gname]]; p <- dim(arr)[1L]; nms <- dimnames(arr)[[1L]]
    for (i in seq_len(p))
      out[[sprintf("logsd:%s:%s", gname, nms[i])]] <- log(sqrt(as.numeric(arr[i, i, , 1L])))
  }
  out
}

ess_bulk_tail <- function(x) {
  c(bulk = tryCatch(posterior::ess_bulk(x), error = function(e) NA_real_),
    tail = tryCatch(posterior::ess_tail(x), error = function(e) NA_real_))
}

## min/median over the finite entries only, NA when none - so a noise-dominated
## (non-positive) differenced per-iter at small n yields a clean NA, not Inf.
safe_agg <- function(f, x) { x <- x[is.finite(x)]; if (length(x)) f(x) else NA_real_ }

# ---- one (tier, warmup) cell ------------------------------------------------

## Fits short and long at matched warmup + seed; differences for per-iter wall;
## reads ESS + mean_leapfrog from the long fit. Returns a one-row data.frame.
measure_cell <- function(tier_name, tier_data, warmup, samp_short, samp_long,
                         seed, n.trees) {
  is_binary <- isTRUE(tier_data$is_binary)
  iter_short <- warmup + samp_short
  iter_long  <- warmup + samp_long

  t_short <- system.time(fit_at(tier_data, warmup, iter_short, seed, n.trees))[["elapsed"]]
  t_long_obj <- NULL
  t_long <- system.time(t_long_obj <- fit_at(tier_data, warmup, iter_long, seed, n.trees))[["elapsed"]]

  per_iter_s <- (t_long - t_short) / (samp_long - samp_short)
  per_iter_ms <- 1000 * per_iter_s
  sampling_wall_s <- if (is.finite(per_iter_s) && per_iter_s > 0) per_iter_s * samp_long else NA_real_

  mono <- monitored_draws(t_long_obj, is_binary)
  ess <- vapply(mono, ess_bulk_tail, numeric(2L))   # 2 x n_monitored
  bulk <- ess["bulk", ]; tail <- ess["tail", ]
  ess_bulk_per_s <- if (!is.na(sampling_wall_s)) bulk / sampling_wall_s else rep(NA_real_, length(bulk))
  ess_tail_per_s <- if (!is.na(sampling_wall_s)) tail / sampling_wall_s else rep(NA_real_, length(tail))

  ml <- t_long_obj$adaptation$mean_leapfrog[[1L]]
  ml_w <- t_long_obj$adaptation$mean_leapfrog_warmup[[1L]]
  nq <- tier_n_q(tier_name)

  data.frame(
    tier = tier_name, n = nq[["n"]], q = nq[["q"]], warmup = warmup,
    samp_short = samp_short, samp_long = samp_long, seed = seed,
    wall_short_s = round(t_short, 4), wall_long_s = round(t_long, 4),
    per_iter_ms = round(per_iter_ms, 4),
    sampling_wall_s = round(sampling_wall_s, 4),
    mean_leapfrog_sampling = round(ml, 3), mean_leapfrog_warmup = round(ml_w, 3),
    n_monitored = length(mono),
    min_ess_bulk = round(safe_agg(min, bulk), 1),
    median_ess_bulk = round(safe_agg(median, bulk), 1),
    min_ess_tail = round(safe_agg(min, tail), 1),
    median_ess_tail = round(safe_agg(median, tail), 1),
    min_ess_bulk_per_s = round(safe_agg(min, ess_bulk_per_s), 2),
    median_ess_bulk_per_s = round(safe_agg(median, ess_bulk_per_s), 2),
    min_ess_tail_per_s = round(safe_agg(min, ess_tail_per_s), 2),
    commit = COMMIT %||% NA_character_, date = as.character(Sys.Date()),
    stringsAsFactors = FALSE)
}

# ---- record over tiers x warmup grid ----------------------------------------

record_warmup <- function(outfile, tier_names, warmup_grid,
                          samp_short, samp_long, n.trees = 50L, seed = 1L) {
  rows <- list()
  for (tier_name in tier_names) {
    cat(sprintf("=== tier %s (building data) ===\n", tier_name))
    tier_data <- build_one(tier_name)
    ss <- samp_short[[tier_name]] %||% samp_short[["default"]]
    sl <- samp_long[[tier_name]] %||% samp_long[["default"]]
    for (w in warmup_grid) {
      cat(sprintf("  warmup=%d (short iter=%d, long iter=%d) ...\n", w, w + ss, w + sl))
      row <- measure_cell(tier_name, tier_data, w, ss, sl, seed, n.trees)
      cat(sprintf("    per-iter %.3f ms | mean-leapfrog samp %.1f warm %.1f | min ESS/s bulk %.2f\n",
                  row$per_iter_ms, row$mean_leapfrog_sampling, row$mean_leapfrog_warmup,
                  row$min_ess_bulk_per_s))
      rows[[length(rows) + 1L]] <- row
      # write incrementally so a late-tier failure keeps completed cells
      write.csv(do.call(rbind, rows), outfile, row.names = FALSE)
    }
    rm(tier_data); gc()
  }
  result <- do.call(rbind, rows)
  cat("\nwrote", outfile, "\n\n")
  print_summary(result)
  invisible(result)
}

print_summary <- function(df) {
  cat("==== bench-warmup summary (per-iter ms | leapfrog samp/warm | min ESS/s bulk) ====\n")
  for (tier_name in unique(df$tier)) {
    sub <- df[df$tier == tier_name, ]
    cat(sprintf("\n%s (n=%d, q=%s):\n", tier_name, sub$n[1L], as.character(sub$q[1L])))
    tbl <- data.frame(
      warmup = sub$warmup,
      per_iter_ms = sub$per_iter_ms,
      lf_samp = sub$mean_leapfrog_sampling,
      lf_warm = sub$mean_leapfrog_warmup,
      minESSs_bulk = sub$min_ess_bulk_per_s,
      medESSs_bulk = sub$median_ess_bulk_per_s)
    print(tbl, row.names = FALSE)
  }
}

# ---- CLI --------------------------------------------------------------------

WARMUP_GRID <- c(20L, 50L, 100L, 200L, 400L, 1000L)
# Sampling-draw counts for the differencing + ESS. The large tier is expensive,
# so it uses fewer sampling draws; the small reference tiers can afford more for
# a steadier ESS estimate.
SAMP_SHORT <- list(default = 200L, large_n_re = 100L)
SAMP_LONG  <- list(default = 1000L, large_n_re = 400L)

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L) {
  mode <- .args[[1L]]
  if (mode == "record") {
    outfile <- if (length(.args) >= 2L) .args[[2L]] else
      file.path(SELF_DIR, "..", "baselines", "bench-warmup-BASELINE.csv")
    tier_names <- if (length(.args) >= 3L) strsplit(.args[[3L]], ",")[[1L]] else names(WARMUP_TIERS)
    record_warmup(outfile, tier_names, WARMUP_GRID, SAMP_SHORT, SAMP_LONG)
  } else if (mode == "quick") {
    outfile <- if (length(.args) >= 2L) .args[[2L]] else tempfile("bench-warmup-quick-", fileext = ".csv")
    record_warmup(outfile, c("continuous_nc1", "continuous_nc2"),
                  warmup_grid = c(20L, 100L),
                  samp_short = list(default = 60L), samp_long = list(default = 200L),
                  n.trees = 25L)
  } else stop("unknown mode: ", mode, " (expected 'record' or 'quick')")
}
