# mvbart: multivariate BART by composition (correlated outcomes / SUR)

# ---------------------------------------------------------------------------
# helpers (also used by the at_home statistical checks)
# ---------------------------------------------------------------------------
cov2rho <- function(S) S[1L, 2L] / sqrt(S[1L, 1L] * S[2L, 2L])
riwishart <- function(df, Psi) solve(rWishart(1L, df, solve(Psi))[, , 1L])
ci <- function(z) unname(quantile(z, c(0.025, 0.975)))

simSUR <- function(n, rho, s1 = 1.0, s2 = 1.3, p = 3L, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  B <- matrix(c(1.5, -1.0, 0.7, -0.6, 1.1, -0.9), p, 2L)
  Sig <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2L, 2L)
  E <- matrix(rnorm(n * 2L), n, 2L) %*% chol(Sig)
  list(X = as.data.frame(X), Y = X %*% B + E, rho = rho)
}

# ---------------------------------------------------------------------------
# (c) shape / finiteness gates on a small, fast fit  (always run)
# ---------------------------------------------------------------------------
d <- simSUR(150L, 0.7, seed = 11L)
df <- data.frame(y1 = d$Y[, 1L], y2 = d$Y[, 2L], d$X)
tr <- 1:120; te <- 121:150

fit <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df[tr, ],
              test = df[te, c("x1", "x2", "x3")],
              n.samples = 120L, n.burn = 50L, n.chains = 2L,
              n.trees = 30L, seed = 42L)

# object class and dimensions
expect_inherits(fit, "mvbartFit")
expect_equal(dim(fit$Sigma), c(2L, 2L, 120L, 2L))
expect_equal(dim(fit$train), c(120L, 2L, 120L, 2L))
expect_equal(dim(fit$test),  c(30L, 2L, 120L, 2L))
expect_equal(dim(fit$mMax),  c(120L, 2L))
expect_equal(fit$outcome.names, c("y1", "y2"))

# everything finite
expect_true(all(is.finite(fit$Sigma)))
expect_true(all(is.finite(fit$train)))
expect_true(all(is.finite(fit$test)))
expect_true(all(is.finite(fit$mMax)))

# every Sigma draw is symmetric and positive definite
S11 <- fit$Sigma[1L, 1L, , ]; S22 <- fit$Sigma[2L, 2L, , ]
S12 <- fit$Sigma[1L, 2L, , ]; S21 <- fit$Sigma[2L, 1L, , ]
expect_equal(S12, S21)
expect_true(all(S11 > 0) && all(S22 > 0))
expect_true(all(S11 * S22 - S12^2 > 0))            # positive determinant

# scale-drift diagnostic stays well anchored (< 1 means no rescale needed)
expect_true(max(fit$mMax) < 1.0)

# default weakly-informative inverse-Wishart prior: df = q + 2, prior mean of
# Sigma equals the diagonal of marginal outcome variances
expect_equal(fit$prior$df, 4)
Y.tr <- as.matrix(df[tr, c("y1", "y2")])
expect_equal(fit$prior$scale / (fit$prior$df - 2 - 1), diag(apply(Y.tr, 2L, var)))

# a univariate response is rejected
expect_error(mvbart(y1 ~ x1 + x2 + x3, data = df[tr, ],
                    n.samples = 10L, n.burn = 5L, n.chains = 1L),
             "multivariate")

# same seed reproduces bitwise; different seed does not
fit.a <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df[tr, ],
                n.samples = 40L, n.burn = 20L, n.chains = 1L, n.trees = 25L,
                seed = 99L)
fit.b <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df[tr, ],
                n.samples = 40L, n.burn = 20L, n.chains = 1L, n.trees = 25L,
                seed = 99L)
fit.c <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df[tr, ],
                n.samples = 40L, n.burn = 20L, n.chains = 1L, n.trees = 25L,
                seed = 100L)
expect_identical(fit.a$Sigma, fit.b$Sigma)
expect_false(identical(fit.a$Sigma, fit.c$Sigma))

# directional sanity: strongly positively correlated data -> positive rho
S <- fit$Sigma
rho.draws <- S[1L, 2L, , ] / sqrt(S[1L, 1L, , ] * S[2L, 2L, , ])
expect_true(mean(rho.draws) > 0.3)


# ---------------------------------------------------------------------------
# (a) SUR-oracle recovery + (b) value case  (heavier; local only)
# ---------------------------------------------------------------------------
if (at_home()) {
  # analytic MNIW oracle for the residual correlation under the same IW prior
  oracleSUR <- function(Y, X, nu0, Psi0, ndraw = 3000L, kappa = 1e-4) {
    Y <- as.matrix(Y); Xd <- cbind(1.0, as.matrix(X))
    n <- nrow(Y); pp <- ncol(Xd)
    Om <- crossprod(Xd) + kappa * diag(pp)
    Bn <- solve(Om, crossprod(Xd, Y))
    Psi_n <- Psi0 + crossprod(Y) - crossprod(Bn, Om %*% Bn)
    vapply(seq_len(ndraw), function(i) cov2rho(riwishart(nu0 + n, Psi_n)),
           numeric(1L))
  }

  # (a) rho recovered at rho in {0, 0.5, 0.9}: agrees with the analytic oracle
  #     and lands within a band of the truth
  for (rho in c(0.0, 0.5, 0.9)) {
    d <- simSUR(350L, rho, seed = 100L + round(100 * rho))
    Y <- cbind(y1 = d$Y[, 1L], y2 = d$Y[, 2L])
    df <- data.frame(Y, d$X)
    nu0 <- 4L; Psi0 <- diag(apply(Y, 2L, var)) * (nu0 - 2L - 1L)

    fit <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df,
                  prior = list(df = nu0, scale = Psi0),
                  n.samples = 350L, n.burn = 150L, n.chains = 1L,
                  n.trees = 60L, seed = 7L)
    S <- fit$Sigma
    rc <- S[1L, 2L, , ] / sqrt(S[1L, 1L, , ] * S[2L, 2L, , ])
    ro <- oracleSUR(Y, d$X, nu0, Psi0)

    # composition posterior agrees with the analytic oracle
    expect_true(abs(mean(rc) - mean(ro)) < 0.05,
                info = paste0("comp-vs-oracle rho at rho=", rho))
    # and recovers the true correlation within a tolerance band
    expect_true(abs(mean(rc) - rho) < 0.15,
                info = paste0("rho recovery at rho=", rho))
    # credible intervals overlap the oracle's
    expect_true(ci(rc)[1L] <= ci(ro)[2L] && ci(rc)[2L] >= ci(ro)[1L],
                info = paste0("CI overlap at rho=", rho))
  }

  # (b) value case: on strongly correlated nonlinear data, the composed joint
  #     log predictive score beats q independent BARTs
  f1 <- function(X) 10 * sin(pi * X[, 1L] * X[, 2L]) + 20 * (X[, 3L] - 0.5)^2 +
                    10 * X[, 4L] + 5 * X[, 5L]
  f2 <- function(X) 14 * cos(pi * X[, 1L] * X[, 3L]) + 18 * (X[, 2L] - 0.5)^2 +
                    8 * X[, 5L] + 6 * X[, 4L]
  simFried <- function(n, rho, s1, s2, p = 6L, seed = 1L) {
    set.seed(seed)
    X <- matrix(runif(n * p), n, p); colnames(X) <- paste0("x", seq_len(p))
    Fm <- cbind(f1(X), f2(X))
    Sig <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2L, 2L)
    list(X = as.data.frame(X), Y = Fm + matrix(rnorm(n * 2L), n, 2L) %*% chol(Sig))
  }
  lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
  jointFull <- function(fTest, Sig, Yt) {
    nT <- nrow(Yt); Sn <- dim(fTest)[3L]; LD <- matrix(0.0, nT, Sn)
    for (s in seq_len(Sn)) {
      Ss <- Sig[, , s]; inv <- solve(Ss); ld <- log(Ss[1L, 1L] * Ss[2L, 2L] - Ss[1L, 2L]^2)
      r1 <- Yt[, 1L] - fTest[, 1L, s]; r2 <- Yt[, 2L] - fTest[, 2L, s]
      qd <- inv[1L, 1L] * r1^2 + 2 * inv[1L, 2L] * r1 * r2 + inv[2L, 2L] * r2^2
      LD[, s] <- -log(2 * pi) - 0.5 * ld - 0.5 * qd
    }
    sum(apply(LD, 1L, lse) - log(Sn))
  }
  jointDiag <- function(fTest, sigma, Yt) {
    nT <- nrow(Yt); Sn <- dim(fTest)[3L]; LD <- matrix(0.0, nT, Sn)
    for (s in seq_len(Sn))
      LD[, s] <- dnorm(Yt[, 1L], fTest[, 1L, s], sigma[s, 1L], log = TRUE) +
                 dnorm(Yt[, 2L], fTest[, 2L, s], sigma[s, 2L], log = TRUE)
    sum(apply(LD, 1L, lse) - log(Sn))
  }
  indepBart <- function(Y, X, Xtest, n.trees, n.burn, n.save, seed = 21L) {
    q <- ncol(Y); nT <- nrow(Xtest)
    ctrl <- dbarts::dbartsControl(n.chains = 1L, n.threads = 1L, n.trees = n.trees,
                                  n.burn = 0L, n.samples = 1L, updateState = FALSE,
                                  verbose = FALSE)
    fTest <- array(0.0, c(nT, q, n.save)); sigma <- matrix(0.0, n.save, q)
    for (k in seq_len(q)) {
      dfk <- data.frame(.y = Y[, k], X)
      s <- eval(bquote(dbarts::dbarts(.y ~ ., data = dfk, test = Xtest,
                                      control = ctrl, seed = .(as.integer(seed + k)))))
      for (it in seq_len(n.burn)) s$run(0L, 1L)
      for (it in seq_len(n.save)) {
        r <- s$run(0L, 1L); fTest[, k, it] <- r$test[, 1L]; sigma[it, k] <- r$sigma[1L]
      }
    }
    list(fTest = fTest, sigma = sigma)
  }

  rho <- 0.8; s1 <- 3.0; s2 <- 3.0
  dTr <- simFried(400L, rho, s1, s2, seed = 2080L)
  dTe <- simFried(500L, rho, s1, s2, seed = 8080L)
  Y <- cbind(y1 = dTr$Y[, 1L], y2 = dTr$Y[, 2L])
  dfTr <- data.frame(Y, dTr$X)

  fit <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3 + x4 + x5 + x6, data = dfTr,
                test = dTe$X, n.samples = 300L, n.burn = 150L, n.chains = 1L,
                n.trees = 75L, seed = 21L)
  ifit <- indepBart(Y, dTr$X, dTe$X, 75L, 150L, 300L, seed = 21L)

  scC <- jointFull(fit$test[, , , 1L], fit$Sigma[, , , 1L], as.matrix(dTe$Y))
  scI <- jointDiag(ifit$fTest, ifit$sigma, as.matrix(dTe$Y))

  # composed model wins the joint predictive by a clear margin
  expect_true(scC - scI > 50,
              info = paste0("joint log-score gain = ", round(scC - scI, 1)))
  # and recovers the strong correlation (up to the small positive shared-X bias)
  rc <- fit$Sigma[1L, 2L, , 1L] / sqrt(fit$Sigma[1L, 1L, , 1L] * fit$Sigma[2L, 2L, , 1L])
  expect_true(abs(mean(rc) - rho) < 0.1,
              info = paste0("value-case rho_hat = ", round(mean(rc), 3)))
}
