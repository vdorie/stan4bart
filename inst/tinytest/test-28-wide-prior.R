# A coefficient prior far wider than the posterior fits as a flat one does, and
# a flat prior fits whatever the units of its column. The forest's covariate is
# unrelated to the response, so the fixed effects' posterior means sit near
# least squares. The prior scales run from about 150 times the default's to
# 1.6e7 times, one centered far from the estimate; the last fit measures x1 in
# thousandths, which makes its coefficient 0.002.

local({
  set.seed(31)
  n <- 200L
  df <- data.frame(x1 = rnorm(n), x2 = rnorm(n), xb = runif(n))
  df$y <- 1 + 2 * df$x1 - 1.5 * df$x2 + rnorm(n, 0, 0.5)
  df_milli <- within(df, x1 <- 1000 * x1)

  fit_fixef <- function(data, ...) {
    fit <- stan4bart(y ~ x1 + x2 + bart(xb), data, ...,
                     chains = 2L, cores = 1L, iter = 400L, warmup = 200L, seed = 7L, verbose = -1L,
                     bart_args = list(n.trees = 20L, n.threads = 1L))
    list(fixef = extract(fit, "fixef", combine_chains = FALSE), data = data)
  }
  fits <- list(
    normal_1e3    = fit_fixef(df, prior = normal(0, 1e3, autoscale = FALSE)),
    normal_1e8    = fit_fixef(df, prior = normal(0, 1e8, autoscale = FALSE)),
    normal_500_1e4 = fit_fixef(df, prior = normal(500, 1e4, autoscale = FALSE)),
    student_t_1e4 = fit_fixef(df, prior = student_t(3, 0, 1e4, autoscale = FALSE)),
    cauchy_1e3    = fit_fixef(df, prior = cauchy(0, 1e3, autoscale = FALSE)),
    flat          = fit_fixef(df, prior = NULL),
    flat_milli    = fit_fixef(df_milli, prior = NULL)
  )
  # every chain within four least-squares standard errors; the posterior sd
  # is about one, the Monte Carlo error of a chain's mean about a tenth
  for (nm in names(fits)) {
    ols <- lm(y ~ x1 + x2, fits[[nm]]$data)
    chain_means <- apply(fits[[nm]]$fixef[c("x1", "x2"), , , drop = FALSE], c(1L, 3L), mean)
    z <- abs(chain_means - coef(ols)[c("x1", "x2")]) / sqrt(diag(vcov(ols)))[c("x1", "x2")]
    expect_true(all(z < 4), info = sprintf("%s: max |mean - lm| / se = %.3g", nm, max(z)))
  }
})
