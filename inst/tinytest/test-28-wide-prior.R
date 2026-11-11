# A coefficient prior far wider than the posterior fits as a flat one does, and
# a flat prior fits whatever the units of its column, with QR as without. The
# forest's covariate is unrelated to the response, so the fixed effects'
# posterior means sit near least squares. The prior scales run from three
# times the default's to 1.6e7 times; three are centered far from the
# estimate, two of them in the student_t family, whose map from the sampler's
# coordinate is not affine. One flat fit measures x1 in thousandths, which
# makes its coefficient 0.002, and the QR fits measure x2 in ten-thousandths,
# the column that sets their basis.

local({
  set.seed(31)
  n <- 200L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), xb = runif(n))
  df$y <- 1 + 2 * df$x1 - 1.5 * df$x2 + rnorm(n, 0, 0.5)
  df_milli <- within(df, x1 <- 1000 * x1)
  df_qr <- within(df, x2 <- 1e4 * x2)

  # the fit's fixed effects and the warnings it raised
  fit_fixef <- function(data, ...) {
    warned <- character()
    fit <- withCallingHandlers(
      stan4bart(y ~ x1 + x2 + bart(xb), data, ...,
                chains = 2L, cores = 1L, iter = 400L, warmup = 200L, seed = 7L, verbose = -1L,
                bart_args = list(n.trees = 20L, n.threads = 1L)),
      warning = function(w) {
        warned <<- c(warned, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(fixef = extract(fit, "fixef", combine_chains = FALSE), data = data, warned = warned)
  }
  fits <- list(
    normal_1e3    = fit_fixef(df, prior = normal(0, 1e3, autoscale = FALSE)),
    normal_1e8    = fit_fixef(df, prior = normal(0, 1e8, autoscale = FALSE)),
    normal_500_1e4 = fit_fixef(df, prior = normal(500, 1e4, autoscale = FALSE)),
    normal_500_20 = fit_fixef(df, prior = normal(500, 20, autoscale = FALSE)),
    student_t_1e4 = fit_fixef(df, prior = student_t(3, 0, 1e4, autoscale = FALSE)),
    student_t_500_1e4 = fit_fixef(df, prior = student_t(3, 500, 1e4, autoscale = FALSE)),
    cauchy_1e3    = fit_fixef(df, prior = cauchy(0, 1e3, autoscale = FALSE)),
    cauchy_500_1e3 = fit_fixef(df, prior = cauchy(500, 1e3, autoscale = FALSE)),
    flat          = fit_fixef(df, prior = NULL),
    flat_milli    = fit_fixef(df_milli, prior = NULL),
    qr_flat       = fit_fixef(df_qr, prior = NULL, QR = TRUE),
    qr_normal_1e6 = fit_fixef(df_qr, prior = normal(0, 1e6, autoscale = FALSE), QR = TRUE)
  )
  # every chain within four least-squares standard errors; the posterior sd
  # is about one, the Monte Carlo error of a chain's mean about a tenth
  for (nm in names(fits)) {
    ols <- lm(y ~ x1 + x2, fits[[nm]]$data)
    chain_means <- apply(fits[[nm]]$fixef[c("x1", "x2"), , , drop = FALSE], c(1L, 3L), mean)
    z <- abs(chain_means - coef(ols)[c("x1", "x2")]) / sqrt(diag(vcov(ols)))[c("x1", "x2")]
    expect_true(all(z < 4), info = sprintf("%s: max |mean - lm| / se = %.3g", nm, max(z)))
    # x2 in ten-thousandths draws the design's scale warning, and only that;
    # x1 in thousandths stays under the warning's threshold
    expected <- if (identical(fits[[nm]]$data, df_qr)) 1L else 0L
    expect_true(length(fits[[nm]]$warned) == expected &&
                  all(startsWith(fits[[nm]]$warned, "Some predictor variables are on very different scales")),
                info = paste(nm, "warnings:", paste(fits[[nm]]$warned, collapse = " | ")))
  }
})

# The raw z_beta rows hold the prior's own standardized coefficient, not the
# sampler's coordinate, where the two differ.
local({
  op <- options(stan4bart.save_raw_parameters = TRUE)
  on.exit(options(op))
  set.seed(32)
  n <- 100L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), xb = runif(n))
  df$y <- 1 + 2 * df$x1 - 1.5 * df$x2 + rnorm(n, 0, 0.5)
  fit <- stan4bart(y ~ x1 + x2 + bart(xb), df, prior = normal(c(500, -3), 1e4, autoscale = FALSE),
                   chains = 1L, cores = 1L, iter = 20L, warmup = 10L, seed = 7L, verbose = -1L,
                   bart_args = list(n.trees = 5L, n.threads = 1L))
  rows <- fit$stan[, , 1L]
  z_beta <- rows[c("z_beta.1", "z_beta.2"), ]
  beta <- rows[c("beta.1", "beta.2"), ]
  expect_true(max(abs(z_beta - (beta - c(500, -3)) / 1e4)) < 1e-12, info = "raw z_beta rows")
})
