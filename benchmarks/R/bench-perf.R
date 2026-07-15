#!/usr/bin/env Rscript
# C0(a) PERF baselines (docs/plans/walnuts.md) - the arc's MOTIVATION,
# measured not asserted. NOT YET RUN as of the C0(b) posterior-baseline
# recording pass (needs a quiet machine with no concurrent load, per
# CLAUDE.local.md's bench-sampler discipline; that is a separate, later
# grant). This script is committed so only the RUN is deferred.
#
# Measures, on three reference fits shared with the C0(b) posterior
# baselines (continuous multilevel = continuous_nc2, binary multilevel =
# binary_multilevel, weighted continuous = weighted_continuous - see
# record-posterior-baselines.R's TIERS):
#   - wall-clock `R CMD INSTALL` time for the package
#   - peak sampling RSS for each reference fit (own subprocess per fit, so
#     RSS reflects that fit alone)
#   - per-iteration sampling wall-time for each reference fit, isolated from
#     fixed per-call setup cost by differencing two runs at different
#     iteration counts (same warmup, longer vs shorter `iter`); the slope
#     (time_long - time_short) / (iter_long - iter_short) is the
#     per-iteration cost.
#
# Usage:
#   Rscript benchmarks/R/bench-perf.R record <outfile.csv> [pkg_path]
#
# Compare the AFTER numbers (C4, WALNUTS-only build) against this recording
# the same way dbarts/benchmarks/R/bench-sampler.R compares record vs
# compare: same machine, no concurrent load, `record` again post-C4 and diff
# by hand (three reference fits is small enough that a bespoke `compare`
# mode did not seem worth the complexity here).
#
# Concrete before-number already in hand as a cheap proxy (see the plan):
# stan_sampler.o compiles to ~50 MB (init.o ~7 MB) - the stanc-generated
# translation unit dominates build cost; C4 must show this gone. This
# script's `R CMD INSTALL` timing is the wall-clock analog of that.

suppressMessages(library(stan4bart))

`%||%` <- function(a, b) if (is.null(a)) b else a

## Directory this script lives in, so it can find its sibling regardless of
## the caller's working directory (Rscript does not set one automatically).
.self_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
SELF_DIR <- if (!is.na(.self_file) && nzchar(.self_file)) dirname(normalizePath(.self_file)) else "benchmarks/R"

# Reuse the exact reference-dataset simulators / tier definitions from the
# C0(b) posterior baselines so the perf fits are the SAME data the
# distributional gate exercises, not a fourth ad hoc dataset.
source(file.path(SELF_DIR, "record-posterior-baselines.R"))

REFERENCE_FITS <- c("continuous_nc2", "binary_multilevel", "weighted_continuous")
FIT_LABELS <- c(continuous_nc2 = "continuous multilevel",
               binary_multilevel = "binary multilevel",
               weighted_continuous = "weighted continuous")

# ---- R CMD INSTALL wall-clock ----------------------------------------------

## Wall-clock time for `R CMD INSTALL --preclean pkg_path` (a clean rebuild,
## matching what a contributor pays after a header change - see
## CLAUDE.local.md's --preclean guidance). Runs INSTALL as a subprocess so
## the *.o/*.so from compiling this very script's session are irrelevant.
time_install <- function(pkg_path) {
  t0 <- proc.time()
  status <- system2("R", c("CMD", "INSTALL", "--preclean", shQuote(pkg_path)),
                    stdout = FALSE, stderr = FALSE)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  if (status != 0L) warning("R CMD INSTALL exited ", status, " - timing may not reflect a clean install")
  elapsed
}

# ---- peak RSS via a timed subprocess ---------------------------------------

## Runs `r_expr` (a string of R code performing one reference fit) in a
## fresh Rscript subprocess under the platform's `/usr/bin/time`, parsing
## wall-clock and peak RSS from its output. One subprocess per fit, so RSS
## reflects that fit alone (not the harness's own footprint).
run_timed_subprocess <- function(r_expr) {
  is_macos <- Sys.info()[["sysname"]] == "Darwin"
  time_flag <- if (is_macos) "-l" else "-v"
  script_file <- tempfile(fileext = ".R")
  writeLines(r_expr, script_file)
  out_file <- tempfile(fileext = ".txt")
  on.exit(unlink(c(script_file, out_file)), add = TRUE)

  t0 <- proc.time()
  system2("/usr/bin/time", c(time_flag, "Rscript", shQuote(script_file)),
         stdout = FALSE, stderr = out_file)
  elapsed <- (proc.time() - t0)[["elapsed"]]

  time_output <- readLines(out_file, warn = FALSE)
  peak_rss_bytes <- if (is_macos) {
    line <- grep("maximum resident set size", time_output, value = TRUE)
    if (length(line) == 0L) NA_real_ else as.numeric(trimws(sub("\\s*maximum.*", "", line[1L])))
  } else {
    line <- grep("Maximum resident set size", time_output, value = TRUE)
    if (length(line) == 0L) NA_real_ else as.numeric(trimws(sub(".*:\\s*", "", line[1L]))) * 1024
  }
  list(elapsed = elapsed, peak_rss_bytes = peak_rss_bytes)
}

## Builds the subprocess script for one reference fit at a given warmup/iter,
## at a single-chain, single-core setting - RSS/per-iteration cost should
## reflect one worker process, matching how stan4bart_fit.R spawns one
## single-threaded chain per cluster worker.
build_fit_script <- function(tier_name, warmup, iter) {
  cfg <- TIERS[[tier_name]]
  sprintf('
suppressMessages(library(stan4bart))
source("%s")
tier_data <- build_tier_data(TIERS[["%s"]])
weights_vec <- if (!is.null(tier_data$weights_col)) tier_data$df[[tier_data$weights_col]] else NULL
data <- tier_data$df
shape <- tier_data$shape
n.trees <- 50L
if (shape == "nc1") {
 stan4bart::stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 | g), data,
                      weights = weights_vec, cores = 1L, verbose = -1L, chains = 1L,
                      warmup = %d, iter = %d, seed = 1L, bart_args = list(n.trees = n.trees))
} else if (shape == "nc2") {
 stan4bart::stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 | g), data,
                      weights = weights_vec, cores = 1L, verbose = -1L, chains = 1L,
                      warmup = %d, iter = %d, seed = 1L, bart_args = list(n.trees = n.trees))
} else if (shape == "nc3") {
 stan4bart::stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 + Xr2 | g), data,
                      weights = weights_vec, cores = 1L, verbose = -1L, chains = 1L,
                      warmup = %d, iter = %d, seed = 1L, bart_args = list(n.trees = n.trees))
} else {
 stan4bart::stan4bart(y ~ bart(X1+X2+X3+X4+X5) + Xfix + (1 + Xr1 | g) + (1 | g2), data,
                      cores = 1L, verbose = -1L, chains = 1L,
                      warmup = %d, iter = %d, seed = 1L, bart_args = list(n.trees = n.trees))
}
',
         file.path(SELF_DIR, "record-posterior-baselines.R"),
         tier_name, warmup, iter, warmup, iter, warmup, iter, warmup, iter)
}

## Per-iteration wall-time via differencing: fit the same warmup at a short
## and a long `iter`, and divide the elapsed-time difference by the
## post-warmup-draw difference. This removes the fixed per-call setup cost
## (data marshalling, Stan model construction) from the per-iteration
## estimate. warmup is short (50) since only the SLOPE across two `iter`
## values is used, not the warmup phase's own cost.
per_iteration_time <- function(tier_name, warmup = 50L, iter_short = 100L, iter_long = 300L) {
  short <- run_timed_subprocess(build_fit_script(tier_name, warmup, iter_short))
  long <- run_timed_subprocess(build_fit_script(tier_name, warmup, iter_long))
  list(per_iteration_s = (long$elapsed - short$elapsed) / (iter_long - iter_short),
      peak_rss_bytes = max(short$peak_rss_bytes, long$peak_rss_bytes, na.rm = TRUE),
      short = short, long = long)
}

# ---- record mode ------------------------------------------------------------

record_perf <- function(outfile, pkg_path = "/Users/vdorie/Repositories/stan4bart/.claude/worktrees/walnuts") {
  install_s <- time_install(pkg_path)
  cat(sprintf("R CMD INSTALL --preclean: %.1fs\n", install_s))

  rows <- list(data.frame(metric = "install_wall_s", fit = NA_character_, value = install_s))
  for (tier_name in REFERENCE_FITS) {
    cat("timing", FIT_LABELS[[tier_name]], "...\n")
    res <- per_iteration_time(tier_name)
    cat(sprintf("  per-iteration: %.4fs  peak RSS: %.0f MB\n",
               res$per_iteration_s, res$peak_rss_bytes / 1e6))
    rows[[length(rows) + 1L]] <- data.frame(metric = "per_iteration_s", fit = tier_name,
                                            value = res$per_iteration_s)
    rows[[length(rows) + 1L]] <- data.frame(metric = "peak_rss_bytes", fit = tier_name,
                                            value = res$peak_rss_bytes)
  }
  result <- do.call(rbind, rows)
  result$commit <- "967a1c6"
  result$date <- as.character(Sys.Date())
  write.csv(result, outfile, row.names = FALSE)
  cat("wrote", outfile, "\n")
  invisible(result)
}

# ---- CLI --------------------------------------------------------------------

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L && .args[[1L]] == "record") {
  outfile <- if (length(.args) >= 2L) .args[[2L]] else "benchmarks/baselines/bench-perf-967a1c6.csv"
  pkg_path <- if (length(.args) >= 3L) .args[[3L]] else "/Users/vdorie/Repositories/stan4bart/.claude/worktrees/walnuts"
  record_perf(outfile, pkg_path)
}
