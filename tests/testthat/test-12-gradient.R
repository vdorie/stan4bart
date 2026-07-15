context("WALNUTS parametric target: gradient FD + Stan log_prob cross-check")

# Correctness gates for src/parametric_model.hpp. Two proofs on
# self-contained fixtures spanning every gradient tier:
#   (a) FD GATE: the hand-written gradient vs central finite differences of the
#       hand-written value, at fixed-seed random unconstrained points.
#   (b) TARGET CROSS-CHECK: value(mine) - Stan's log_prob<false, true> must be
#       CONSTANT across points within each fixture (targets differ only by the
#       constants the hand target drops).
#
# The Stan data lists are hand-built (not from the full pipeline) so the gate
# tests the target math directly; the SAME list and point feed both the hand
# functor and Stan's own log_prob, so any valid list is a fair oracle.

C_grad <- function(d, par) .Call(stan4bart:::C_stan4bart_logdensity_grad, d, par)
C_stan <- function(d, par) .Call(stan4bart:::C_stan4bart_stan_logdensity, d, par)

# ---- build a marshaled Stan data list -------------------------------------
#
# blocks: list of list(nc =, group = <int 1..n_levels>, slopes = N x (nc-1) or
# NULL). Z columns follow make_b's convention: block bi, level j (0-based), coef
# c -> col boff[bi] + j*nc + c, so the hand functor's make_b and this Z agree.
build_stan_data <- function(N, X, y, blocks, is_binary,
                            weights = NULL,
                            prior_dist = 1L, prior_scale = NULL, prior_df = NULL,
                            prior_dist_for_aux = 3L, prior_scale_for_aux = 1.0,
                            prior_mean_for_aux = 0.0, prior_df_for_aux = 0.0) {
  K <- ncol(X)
  t <- length(blocks)
  if (is.null(prior_scale)) prior_scale <- rep(2.5, K)
  if (is.null(prior_df))    prior_df    <- rep(3, K)

  if (t > 0L) {
    p <- vapply(blocks, function(b) as.integer(b$nc), integer(1))
    l <- vapply(blocks, function(b) as.integer(max(b$group)), integer(1))
    boff <- c(0L, cumsum(p * l))
    cols <- vector("list", N); vals <- vector("list", N)
    for (i in seq_len(N)) {
      ci <- integer(0); vi <- numeric(0)
      for (bi in seq_len(t)) {
        b <- blocks[[bi]]; nc <- b$nc
        lev <- b$group[i] - 1L
        design <- c(1, if (nc > 1L) as.numeric(b$slopes[i, ]) else numeric(0))
        base <- boff[bi] + lev * nc
        ci <- c(ci, base + 0:(nc - 1L)); vi <- c(vi, design)
      }
      o <- order(ci); cols[[i]] <- ci[o]; vals[[i]] <- vi[o]
    }
    w <- unlist(vals); v <- as.integer(unlist(cols))
    u <- c(0L, cumsum(vapply(cols, length, integer(1))))
    q <- sum(p * l)
    len_theta_L <- sum(choose(p, 2) + p)
    len_concentration <- sum(p[p > 1])
    len_regularization <- sum(p > 1)
    shape <- rep(1, t); scale <- rep(1, t)
    concentration <- rep(1, len_concentration)
    regularization <- rep(1, len_regularization)
  } else {
    p <- integer(0); l <- integer(0); q <- 0L
    w <- numeric(0); v <- integer(0); u <- rep(0L, N + 1L)
    len_theta_L <- 0L; len_concentration <- 0L; len_regularization <- 0L
    shape <- numeric(0); scale <- numeric(0)
    concentration <- numeric(0); regularization <- numeric(0)
  }

  has_weights <- !is.null(weights)
  list(
    N = as.integer(N), K = as.integer(K), X = X,
    len_y = as.integer(N), lb_y = -Inf, ub_y = Inf, y = as.numeric(y),
    has_intercept = 0L, is_binary = as.integer(is_binary),
    prior_dist = as.integer(prior_dist), prior_dist_for_intercept = 0L,
    prior_dist_for_aux = as.integer(prior_dist_for_aux),
    has_weights = as.integer(has_weights),
    weights = if (has_weights) as.numeric(weights) else numeric(0),
    offset_ = numeric(N),
    prior_scale = as.numeric(prior_scale), prior_scale_for_intercept = 1.0,
    prior_scale_for_aux = as.numeric(prior_scale_for_aux),
    prior_mean = rep(0, K), prior_mean_for_intercept = 0.0,
    prior_mean_for_aux = as.numeric(prior_mean_for_aux),
    prior_df = as.numeric(prior_df), prior_df_for_intercept = 1.0,
    prior_df_for_aux = as.numeric(prior_df_for_aux),
    global_prior_df = 1.0, global_prior_scale = 1.0, slab_df = 1.0, slab_scale = 1.0,
    num_normals = integer(0),
    t = as.integer(t), p = as.integer(p), l = as.integer(l), q = as.integer(q),
    len_theta_L = as.integer(len_theta_L),
    shape = as.numeric(shape), scale = as.numeric(scale),
    len_concentration = as.integer(len_concentration),
    concentration = as.numeric(concentration),
    len_regularization = as.integer(len_regularization),
    regularization = as.numeric(regularization),
    num_non_zero = as.integer(length(w)),
    w = as.numeric(w), v = as.integer(v), u = as.integer(u)
  )
}

# unconstrained-vector dimension implied by the geometry
model_dim <- function(d) {
  p <- d$p; t <- d$t
  len_rho <- sum(p) - t
  len_z_T <- sum(vapply(p, function(pi) if (pi >= 3L) (pi - 2L) * (pi - 1L) else 0L, integer(1)))
  d$K + d$q + len_z_T + len_rho + d$len_concentration + t + (if (d$is_binary) 0L else 1L)
}

# a random unconstrained point with well-conditioned per-segment scales
gen_point <- function(d, seed) {
  set.seed(seed)
  p <- d$p; t <- d$t
  len_rho <- sum(p) - t
  len_z_T <- sum(vapply(p, function(pi) if (pi >= 3L) (pi - 2L) * (pi - 1L) else 0L, integer(1)))
  c(rnorm(d$K, 0, 0.7),                # z_beta
    rnorm(d$q, 0, 0.7),                # z_b
    rnorm(len_z_T, 0, 1.0),           # z_T
    rnorm(len_rho, 0, 0.9),           # rho_free
    rnorm(d$len_concentration, 0, 0.5),# zeta_free
    rnorm(t, 0, 0.5),                 # tau_free
    if (d$is_binary) numeric(0) else rnorm(1, 0, 0.4))  # aux_unscaled_free
}

# ---- fixtures: one per tier -----------------------------------------------
make_fixtures <- function() {
  set.seed(20260715L)
  N <- 50L
  X <- matrix(rnorm(N * 2L), N, 2L)
  yc <- rnorm(N, 0, 2)               # continuous response
  yl <- rnorm(N)                     # probit latents (binary)
  wts <- runif(N, 0.4, 2.5)
  g1 <- sample.int(6L, N, replace = TRUE)   # 6 levels
  g2 <- sample.int(4L, N, replace = TRUE)   # 4 levels
  s1 <- matrix(rnorm(N), N, 1L)             # one slope
  s2 <- matrix(rnorm(N * 2L), N, 2L)        # two slopes (nc=3)

  nc1 <- list(nc = 1L, group = g1, slopes = NULL)
  nc2 <- list(nc = 2L, group = g1, slopes = s1)
  nc3 <- list(nc = 3L, group = g1, slopes = s2)
  nc1b <- list(nc = 1L, group = g2, slopes = NULL)

  list(
    fe          = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(), FALSE)),
    nc1         = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc1), FALSE)),
    nc2         = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc2), FALSE)),
    nc3_onion   = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc3), FALSE)),
    nc2_weighted= list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc2), FALSE, weights = wts)),
    # continuous multi-block: the aux -> dispersion coupling must ACCUMULATE
    # across heterogeneous blocks (binary covers multi-block without aux)
    nc_mixed    = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc2, nc1b), FALSE)),
    binary      = list(is_binary = TRUE,  data = build_stan_data(N, X, yl, list(nc2, nc1b), TRUE)),
    studentt_beta = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc2), FALSE, prior_dist = 2L)),
    normal_aux  = list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc2), FALSE,
                                                                 prior_dist_for_aux = 1L,
                                                                 prior_scale_for_aux = 1.5,
                                                                 prior_mean_for_aux = 0.5)),
    studentt_aux= list(is_binary = FALSE, data = build_stan_data(N, X, yc, list(nc1), FALSE,
                                                                 prior_dist_for_aux = 2L,
                                                                 prior_scale_for_aux = 1.2,
                                                                 prior_df_for_aux = 4.0))
  )
}

fixtures <- make_fixtures()
n_points <- 4L
h <- 1e-5

# ---- (a) FD GATE -----------------------------------------------------------
# central differences of the hand value vs the hand gradient; require
# max |fd - grad| / max(|grad|, 1) < fd_tol per fixture.
fd_tol <- 1e-6

test_that("gradient matches central finite differences on every tier", {
  worst_overall <- 0
  for (nm in names(fixtures)) {
    d <- fixtures[[nm]]$data
    dim <- model_dim(d)
    worst <- 0
    for (pt in seq_len(n_points)) {
      par <- gen_point(d, 20260715L + pt * 101L + which(names(fixtures) == nm))
      g <- C_grad(d, par)$gradient
      fd <- numeric(dim)
      for (j in seq_len(dim)) {
        pp <- par; pp[j] <- pp[j] + h
        pm <- par; pm[j] <- pm[j] - h
        fd[j] <- (C_grad(d, pp)$value - C_grad(d, pm)$value) / (2 * h)
      }
      rel <- abs(fd - g) / pmax(abs(g), 1)
      worst <- max(worst, max(rel))
    }
    cat(sprintf("  FD [%-13s] dim=%2d worst rel err = %.2e\n", nm, dim, worst))
    expect_lt(worst, fd_tol)
    worst_overall <- max(worst_overall, worst)
  }
  cat(sprintf("  FD worst across all tiers = %.2e (tol %.0e)\n", worst_overall, fd_tol))
})

# ---- (b) TARGET CROSS-CHECK ------------------------------------------------
# value(mine) - Stan log_prob<false, true> must be constant across points.
xc_tol <- 1e-6

test_that("hand target minus Stan log_prob is constant within each tier", {
  worst_spread <- 0
  for (nm in names(fixtures)) {
    d <- fixtures[[nm]]$data
    diffs <- numeric(n_points)
    for (pt in seq_len(n_points)) {
      par <- gen_point(d, 424242L + pt * 97L + which(names(fixtures) == nm))
      mine <- C_grad(d, par)$value
      stan <- C_stan(d, par)
      diffs[pt] <- mine - stan
    }
    spread <- diff(range(diffs))
    cat(sprintf("  XC [%-13s] const=%+.6f spread=%.2e\n", nm, mean(diffs), spread))
    expect_lt(spread, xc_tol)
    worst_spread <- max(worst_spread, spread)
  }
  cat(sprintf("  XC worst spread across all tiers = %.2e (tol %.0e)\n", worst_spread, xc_tol))
})
