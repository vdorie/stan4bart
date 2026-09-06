#!/usr/bin/env Rscript
# Random-effect SCALE mixing harness (docs/plans/re-scale-mixing.md).
#
# The quantity under test is the random-effect standard deviation, the one
# parametric scalar the sampler does not mix: the non-centered block ties the
# log scale to the standardized effects on a product ridge, and a diagonal
# metric cannot rescale a rotated direction. For each design x skip x seed it
# reports
#   - the scale's lag-1 autocorrelation and ESS over the kept draws, which is
#     what the release bar is stated in (acf1 < 0.8, ESS >= 100 per 1000),
#   - its posterior mean, which catches a chain that walked off along the
#     shared level rather than mixing (the level and the scale move together
#     because Sigma is the spread of the effects about zero, not about their
#     own mean),
#   - the shared level's OWN autocorrelation and ESS. That is a second, separate
#     failure and no skip reaches it: the level of the group intercepts is
#     confounded with the forest's overall level, so it is the two-block
#     alternation that has to cross it, not the parametric transition. Freezing
#     the forest takes it from acf1 ~ 0.98 to ~ 0.09 while leaving the scale at
#     ~ 0.95, which is what separates the two,
#   - the group-contrast ESS, the forest RMSE, and wall time, so a skip that
#     buys scale mixing at another quantity's expense is visible.
#
# Designs are the three grouped fits the dbarts hand-off was judged on: a
# 20-group and a 5-group Gaussian random intercept and a 20-group probit. They
# regenerate from fixed seeds; no data files.
#
# TIMING DISCIPLINE: wall time is reported but is not the point of this
# harness; run it on a quiet machine only if the seconds column is to be read.
#
# Usage:
#   Rscript benchmarks/R/bench-re-scale.R record [outfile.csv] [ds1,ds2,...]
#   Rscript benchmarks/R/bench-re-scale.R quick  [outfile.csv]

suppressMessages({
  library(stan4bart)
  library(posterior)
})

.self_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L])
SELF_DIR <- if (!is.na(.self_file) && nzchar(.self_file)) dirname(normalizePath(.self_file)) else "benchmarks/R"

COMMIT <- tryCatch(system2("git", c("-C", SELF_DIR, "rev-parse", "--short", "HEAD"),
                           stdout = TRUE, stderr = FALSE)[1L], error = function(e) NA_character_)

# ---- the bar ---------------------------------------------------------------

BAR_ACF1 <- 0.8     # scale lag-1 autocorrelation must be below this
BAR_ESS  <- 100     # scale ESS per 1000 kept draws must be at least this

# ---- designs ---------------------------------------------------------------

friedman <- function(x)
  10 * sin(pi * x[, 1L] * x[, 2L]) + 20 * (x[, 3L] - 0.5)^2 +
  10 * x[, 4L] + 5 * x[, 5L]

## n must be a multiple of K; equal group sizes keep the identification
## strength - observations per group - the only thing the designs vary.
simulate_design <- function(seed, n, p, K, tau, binary) {
  set.seed(seed)
  x <- matrix(runif(n * p), n, p)
  fx <- friedman(x)
  if (binary) fx <- (fx - mean(fx)) / sd(fx)
  g <- rep(seq_len(K), each = n %/% K)
  b <- rnorm(K, 0, tau)
  eta <- fx + b[g] + rnorm(n, 0, 1)
  df <- as.data.frame(x)
  colnames(df) <- paste0("x", seq_len(p))
  df$y <- if (binary) as.numeric(eta > 0) else eta
  df$g <- factor(g)
  list(df = df, f = fx, b = b, tau = tau, K = K, binary = binary)
}

DESIGNS <- list(
  gaussian_k20 = list(seed = 101L, n = 2000L, p = 10L, K = 20L, tau = 1.0, binary = FALSE),
  gaussian_k5  = list(seed = 102L, n = 2000L, p = 10L, K =  5L, tau = 2.0, binary = FALSE),
  probit_k20   = list(seed = 103L, n = 2000L, p = 10L, K = 20L, tau = 0.5, binary = TRUE))

SKIP_GRID <- c(1L, 8L, 16L)
SEEDS     <- c(20260907L, 7L, 42L, 1234L, 99L)

# ---- one fit ---------------------------------------------------------------

acf1 <- function(x) as.numeric(acf(as.numeric(x), lag.max = 1L, plot = FALSE)$acf[2L])
ess  <- function(x) tryCatch(posterior::ess_basic(as.numeric(x)),
                             error = function(e) NA_real_)

## The bart() term has to be written out literally - a formula built by
## reformulate() is not recognized as the bart component.
fit_one <- function(d, skip_stan, seed, iter, warmup, n.trees) {
  set.seed(seed)
  tm <- system.time(fit <- stan4bart(
    y ~ bart(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10) + (1 | g),
    data = d$df, iter = iter, warmup = warmup, chains = 1L, cores = 1L,
    verbose = -1L, skip = c(bart = 1L, stan = skip_stan),
    bart_args = list(n.trees = n.trees)))
  Sigma <- extract(fit, type = "Sigma")[["g"]]
  tau <- sqrt(as.numeric(Sigma[1L, 1L, ]))
  ranef <- t(extract(fit, type = "ranef")[["g"]][1L, , , drop = TRUE])
  level <- rowMeans(ranef)
  contrasts <- ranef - level
  data.frame(
    tau_acf1 = acf1(tau),
    tau_ess  = ess(tau) * 1000 / length(tau),
    tau_mean = mean(tau),
    tau_sd   = sd(tau),
    level    = mean(level),
    level_acf1 = acf1(level),
    level_ess  = ess(level) * 1000 / length(tau),
    contrast_ess = median(apply(contrasts, 2L, ess)) * 1000 / length(tau),
    f_rmse   = sqrt(mean((as.numeric(fitted(fit, type = "indiv.bart")) - d$f)^2)),
    mean_leapfrog = as.numeric(fit$adaptation$mean_leapfrog),
    seconds  = tm[["elapsed"]])
}

# ---- driver ----------------------------------------------------------------

record_re_scale <- function(outfile, design_names, skips, seeds,
                            iter = 2000L, warmup = 1000L, n.trees = 200L) {
  rows <- list()
  for (nm in design_names) {
    spec <- DESIGNS[[nm]]
    if (is.null(spec)) stop("unknown design: ", nm)
    d <- do.call(simulate_design, spec)
    for (k in skips) for (s in seeds) {
      r <- fit_one(d, k, s, iter, warmup, n.trees)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(commit = COMMIT, design = nm, K = spec$K, true_tau = spec$tau,
                   skip_stan = k, seed = s, iter = iter, warmup = warmup,
                   n.trees = n.trees, stringsAsFactors = FALSE), r)
      cat(sprintf("%-13s skip %2d seed %8d: acf1 %.3f ess %6.1f mean %.3f | level acf1 %.3f ess %5.1f (%.1f s)\n",
                  nm, k, s, r$tau_acf1, r$tau_ess, r$tau_mean,
                  r$level_acf1, r$level_ess, r$seconds))
      utils::write.csv(do.call(rbind, rows), outfile, row.names = FALSE)
    }
  }
  res <- do.call(rbind, rows)
  cat("\nwrote ", outfile, "\n\n", sep = "")
  report(res)
  invisible(res)
}

## Verdict per design x skip, over seeds: the bar is stated on the WORST seed,
## not the median, because a single chain is what a user runs.
report <- function(res) {
  for (nm in unique(res$design)) for (k in sort(unique(res$skip_stan))) {
    z <- res[res$design == nm & res$skip_stan == k, ]
    if (nrow(z) == 0L) next
    pass <- max(z$tau_acf1) < BAR_ACF1 && min(z$tau_ess) >= BAR_ESS
    cat(sprintf("%-13s skip %2d  acf1 med %.3f worst %.3f | ess med %6.1f worst %6.1f | tau_mean %.3f | level acf1 %.3f ess %5.1f | contrast ess %6.1f | %5.2f s | %s\n",
                nm, k, median(z$tau_acf1), max(z$tau_acf1),
                median(z$tau_ess), min(z$tau_ess), mean(z$tau_mean),
                median(z$level_acf1), median(z$level_ess),
                median(z$contrast_ess), mean(z$seconds),
                if (pass) "PASS" else "fail"))
  }
  cat(sprintf("\nbar: worst-seed acf1 < %.1f and worst-seed ESS >= %d per 1000 kept draws\n",
              BAR_ACF1, BAR_ESS))
  cat("the bar is stated on the scale only; the level columns are reported\n",
      "because no skip reaches them and a reader should not mistake a passing\n",
      "scale for a well-mixed random-effect block.\n", sep = "")
}

.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1L) {
  mode <- .args[[1L]]
  if (mode == "record") {
    outfile <- if (length(.args) >= 2L) .args[[2L]] else
      file.path(SELF_DIR, "..", "baselines", "re-scale-BASELINE.csv")
    design_names <- if (length(.args) >= 3L) strsplit(.args[[3L]], ",")[[1L]] else names(DESIGNS)
    record_re_scale(outfile, design_names, SKIP_GRID, SEEDS)
  } else if (mode == "quick") {
    outfile <- if (length(.args) >= 2L) .args[[2L]] else
      tempfile("re-scale-quick-", fileext = ".csv")
    record_re_scale(outfile, "gaussian_k20", c(1L, 8L), SEEDS[1:2],
                    iter = 800L, warmup = 400L, n.trees = 50L)
  } else stop("unknown mode: ", mode, " (expected 'record' or 'quick')")
}
