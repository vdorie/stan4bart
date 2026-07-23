## Multivariate BART by composition (seemingly-unrelated regressions, SUR).
##
## Jointly models q continuous outcomes that share a predictor set but have
## correlated residuals,
##
##     Y_ik = f_k(x_i) + e_ik,   e_i = (e_i1, ..., e_iq) ~ N_q(0, Sigma).
##
## Each mean surface f_k is an independent dbarts sum-of-trees forest; the
## forests are coupled ONLY across equations, through the residual covariance
## Sigma.  Conditional on the other equations' current residuals, updating one
## forest is an ordinary univariate BART step: the outcome is Gaussian with a
## per-observation offset (the conditional mean) and a scalar conditional
## variance,
##
##     Y_ik | . ~ N( f_k(x_i) + m_ik, v_k ),
##       m_ik = Sigma_{k,-k} Sigma_{-k,-k}^{-1} e_{i,-k},
##       v_k  = Sigma_kk - Sigma_{k,-k} Sigma_{-k,-k}^{-1} Sigma_{-k,k}.
##
## m and v are injected into the k-th sampler through dbarts's setOffset() and
## setSigma() verbs; each sampler is created with resid.prior = fixed() so the
## engine never draws its own nugget -- Sigma owns ALL covariance and is drawn
## conjugately each sweep from an inverse-Wishart on the stacked residual
## cross-product.  The whole routine is built on dbarts's public R interface; it
## adds no compiled code.
##
## This is a PURE-SUR sampler: it has no random-effect / multilevel component.
## Coupling it with stan4bart's WALNUTS multilevel machinery (SUR x ranef, i.e. a
## multivariate mixed model) is a documented future extension and is NOT built
## here.
##
## Load-bearing dbarts facts this relies on (validated against a conjugate SUR
## oracle):
##   * resid.prior = fixed(v) reads v as a VARIANCE; setSigma(s) takes s as a
##     STANDARD DEVIATION, and overrides the fixed() creation value every sweep.
##   * setOffset(m) injects the per-observation conditional mean; the run() train
##     fit INCLUDES that offset, so the forest is recovered as f_k = train - m
##     and the structural residual as e_k = Y_k - f_k.
##   * offset.test = 0 at creation pins the test slot to the pure forest
##     f_k(x*), so a per-sweep train offset never leaks into test predictions.
##   * The response Y_k is fixed (only the offset moves), so the creation-time
##     [-0.5, 0.5] rescale anchor stays put; setResponse()/updateScale are never
##     invoked.  The per-sweep max|offset| / range(Y_k) diagnostic (mMax) is
##     exposed so heavy-tailed cases that would need a warmup-only rescale can be
##     caught.

mvbart <-
  function(formula,
           data = NULL,
           subset,
           na.action = getOption("na.action", "na.omit"),
           test = NULL,
           n.samples = 1000L,
           n.burn = 1000L,
           n.chains = 4L,
           n.trees = 75L,
           prior = NULL,
           seed = NA_integer_,
           verbose = FALSE,
           bart_args = list())
{
  call <- match.call()

  ## ---- model frame (multivariate LHS, e.g. cbind(y1, y2) ~ x) -------------
  if (!is.null(data)) data <- as.data.frame(data)
  mfCall <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "subset", "na.action"), names(mfCall), 0L)
  mfCall <- mfCall[c(1L, keep)]
  mfCall$drop.unused.levels <- TRUE
  mfCall[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfCall, parent.frame())
  mt <- attr(mf, "terms")

  Y <- model.response(mf)
  if (is.null(dim(Y)) || ncol(Y) < 2L)
    stop("'mvbart' requires a multivariate response with at least two outcomes, ",
         "e.g. cbind(y1, y2) ~ x")
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  n <- nrow(Y)
  q <- ncol(Y)

  outcome.names <- colnames(Y)
  if (is.null(outcome.names)) outcome.names <- paste0("y", seq_len(q))

  Xdf <- mf[-attr(mt, "response")]
  if (ncol(Xdf) == 0L) stop("model has no predictors")

  ## ---- test predictors ---------------------------------------------------
  hasTest <- !is.null(test)
  Xtest <- NULL
  nTest <- 0L
  if (hasTest) {
    factorCols <- vapply(Xdf, is.factor, logical(1L))
    xlev <- lapply(Xdf[factorCols], levels)
    Xtest <- model.frame(delete.response(mt), as.data.frame(test),
                         na.action = na.action,
                         xlev = if (length(xlev)) xlev else NULL)
    nTest <- nrow(Xtest)
  }

  ## ---- inverse-Wishart prior on Sigma (weakly-informative defaults) -------
  ## Density proportional to |S|^{-(df+q+1)/2} exp(-tr(scale S^{-1})/2), so
  ## E[Sigma] = scale / (df - q - 1).  The default df = q + 2 is the smallest
  ## integer giving a finite mean; the default scale then centres E[Sigma] on the
  ## marginal outcome variances.
  nu0 <- if (!is.null(prior[["df"]])) prior[["df"]] else q + 2
  if (nu0 <= q - 1)
    stop("inverse-Wishart prior 'df' must exceed q - 1 = ", q - 1L)
  Psi0 <- if (!is.null(prior[["scale"]])) as.matrix(prior[["scale"]])
          else base::diag(apply(Y, 2L, var), q) * (nu0 - q - 1)
  if (nrow(Psi0) != q || ncol(Psi0) != q)
    stop("inverse-Wishart prior 'scale' must be a ", q, " x ", q, " matrix")

  yRange <- apply(Y, 2L, function(z) diff(range(z)))

  ## ---- shared dbarts control --------------------------------------------
  ## Each sweep advances every forest by a single Gibbs step (n.burn = 0,
  ## n.samples = 1); the OUTER loop here supplies burn-in and thinning.
  ctrlCall <- quote(dbarts::dbartsControl(n.chains = 1L, n.threads = 1L,
                                          n.trees = n.trees, n.burn = 0L,
                                          n.samples = 1L, updateState = FALSE,
                                          verbose = FALSE))
  ctrlFormals <- names(formals(dbarts::dbartsControl))
  for (nm in intersect(names(bart_args), setdiff(ctrlFormals, names(ctrlCall))))
    ctrlCall[[nm]] <- bart_args[[nm]]
  control <- eval(ctrlCall)

  nIter <- as.integer(n.burn) + as.integer(n.samples)

  ## ---- one coupled Gibbs chain ------------------------------------------
  run_chain <- function(chainSeed) {
    if (!is.na(chainSeed)) set.seed(chainSeed)

    ## One sampler per equation.  fixed(var(Y_k)) suppresses the engine nugget
    ## (setSigma overrides it before the first run); offset.test = 0 keeps the
    ## test slot equal to the pure forest fit.  Prior tokens fixed()/normal() are
    ## resolved by dbarts's non-standard evaluation, so the call is built and
    ## eval()'d rather than invoked with pre-computed prior objects.
    samplers <- vector("list", q)
    for (k in seq_len(q)) {
      seedk <- if (is.na(chainSeed)) NA_integer_
               else as.integer((chainSeed %% (.Machine$integer.max - q - 1L)) + k)
      dfk <- data.frame(.y = Y[, k], Xdf, check.names = FALSE)
      vk <- var(Y[, k])
      dbCall <- bquote(dbarts::dbarts(.y ~ ., data = dfk,
                                      resid.prior = fixed(.(vk)),
                                      control = control, seed = .(seedk)))
      if (hasTest) {
        dbCall$test <- quote(Xtest)
        dbCall$offset.test <- quote(rep(0.0, nTest))
      }
      if (!is.null(bart_args[["k"]]))
        dbCall$node.prior <- bquote(normal(.(bart_args[["k"]])))
      samplers[[k]] <- eval(dbCall)
    }

    ## Initialise structural residuals with one uncoupled sweep per equation,
    ## then seed Sigma from their empirical (prior-shrunk) cross-product.
    E <- matrix(0.0, n, q)
    for (k in seq_len(q)) {
      samplers[[k]]$setSigma(sd(Y[, k]))
      res <- samplers[[k]]$run(0L, 1L)
      E[, k] <- Y[, k] - res$train[, 1L]
    }
    Sigma <- (crossprod(E) + Psi0) / (n + nu0)

    SigmaDraws <- array(0.0, c(q, q, n.samples))
    fTrain <- array(0.0, c(n, q, n.samples))
    fTest <- if (hasTest) array(0.0, c(nTest, q, n.samples)) else NULL
    mMax <- numeric(n.samples)
    driftMax <- 0.0

    for (it in seq_len(nIter)) {
      drift <- 0.0

      ## Forest sweep: update each equation against fresh residuals.
      for (k in seq_len(q)) {
        oth <- seq_len(q)[-k]
        Sko <- Sigma[k, oth]
        beta <- solve(Sigma[oth, oth, drop = FALSE], Sko)   # (q-1) vector
        m <- as.numeric(E[, oth, drop = FALSE] %*% beta)    # per-obs offset
        v <- max(Sigma[k, k] - sum(Sko * beta), 1e-8)       # conditional var

        samplers[[k]]$setOffset(m)
        samplers[[k]]$setSigma(sqrt(v))
        res <- samplers[[k]]$run(0L, 1L)

        fk <- res$train[, 1L] - m           # train INCLUDES offset -> subtract
        E[, k] <- Y[, k] - fk               # structural residual
        drift <- max(drift, max(abs(m)) / yRange[k])

        if (it > n.burn) {
          s <- it - n.burn
          fTrain[, k, s] <- fk
          if (hasTest) fTest[, k, s] <- res$test[, 1L]      # offset.test == 0
        }
      }
      driftMax <- max(driftMax, drift)

      ## Covariance draw: inverse-Wishart on the stacked residual cross-product.
      Sigma <- .riwishart(nu0 + n, Psi0 + crossprod(E))

      if (it > n.burn) {
        SigmaDraws[, , it - n.burn] <- Sigma
        mMax[it - n.burn] <- drift
      }
      if (is.numeric(verbose) && verbose > 0 &&
          it %% max(1L, as.integer(verbose)) == 0L)
        cat(sprintf("  sweep %4d/%4d  drift = %.3f\n", it, nIter, drift))
    }

    list(Sigma = SigmaDraws, train = fTrain, test = fTest,
         mMax = mMax, driftMax = driftMax)
  }

  ## ---- chain seeds & sequential run --------------------------------------
  haveSeed <- !is.na(seed)
  oldSeed <- if (exists(".Random.seed", envir = .GlobalEnv))
               get(".Random.seed", envir = .GlobalEnv) else NULL
  if (haveSeed) {
    set.seed(seed)
    chainSeeds <- sample.int(.Machine$integer.max - q - 2L, n.chains)
  } else {
    chainSeeds <- rep(NA_integer_, n.chains)
  }

  chainResults <- vector("list", n.chains)
  for (ci in seq_len(n.chains)) {
    if (is.numeric(verbose) && verbose > 0 && n.chains > 1L)
      cat(sprintf("chain %d/%d\n", ci, n.chains))
    chainResults[[ci]] <- run_chain(chainSeeds[ci])
  }

  ## Restore the caller's RNG state (chains set the global seed for the
  ## inverse-Wishart draws).
  if (haveSeed) {
    if (!is.null(oldSeed)) assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv))
      rm(".Random.seed", envir = .GlobalEnv)
  }

  ## ---- assemble output ---------------------------------------------------
  chainNames <- paste0("chain:", seq_len(n.chains))
  Sigma <- array(unlist(lapply(chainResults, `[[`, "Sigma")),
                 dim = c(q, q, n.samples, n.chains),
                 dimnames = list(outcome.names, outcome.names,
                                 iterations = NULL, chain = chainNames))
  train <- array(unlist(lapply(chainResults, `[[`, "train")),
                 dim = c(n, q, n.samples, n.chains),
                 dimnames = list(observation = NULL, outcome = outcome.names,
                                 iterations = NULL, chain = chainNames))
  test <- NULL
  if (hasTest)
    test <- array(unlist(lapply(chainResults, `[[`, "test")),
                  dim = c(nTest, q, n.samples, n.chains),
                  dimnames = list(observation = NULL, outcome = outcome.names,
                                  iterations = NULL, chain = chainNames))
  mMax <- matrix(unlist(lapply(chainResults, `[[`, "mMax")), n.samples, n.chains,
                 dimnames = list(iterations = NULL, chain = chainNames))
  driftMax <- max(vapply(chainResults, `[[`, numeric(1L), "driftMax"))

  if (driftMax > 1.0)
    warning("per-sweep conditional-mean offset exceeded the outcome range ",
            "(max drift ", round(driftMax, 2L), "); the fixed [-0.5, 0.5] ",
            "response scale may be mis-anchored. Consider a longer warmup or a ",
            "warmup-only response rescale.", call. = FALSE)

  result <- list(Sigma = Sigma, train = train, test = test, mMax = mMax,
                 y = Y, outcome.names = outcome.names,
                 n.samples = n.samples, n.burn = n.burn, n.chains = n.chains,
                 prior = list(df = nu0, scale = Psi0),
                 call = call, formula = formula, terms = mt)
  class(result) <- "mvbartFit"
  result
}

## Draw one q x q matrix from the inverse-Wishart IW(df, scale) whose mean is
## scale / (df - q - 1); implemented as a Wishart on the inverse.
.riwishart <- function(df, scale) {
  W <- stats::rWishart(1L, df = df, Sigma = solve(scale))[, , 1L]
  solve(W)
}

print.mvbartFit <- function(x, ...) {
  cat("mvbart fit (multivariate BART by composition)\n")
  cat("  outcomes:  ", paste(x$outcome.names, collapse = ", "), "\n", sep = "")
  cat("  draws:     ", x$n.samples, " x ", x$n.chains, " chain(s), ", x$n.burn,
      " burn-in\n", sep = "")
  Sbar <- apply(x$Sigma, c(1L, 2L), mean)
  cat("  posterior mean residual covariance Sigma:\n")
  print(round(Sbar, 4L))
  invisible(x)
}
