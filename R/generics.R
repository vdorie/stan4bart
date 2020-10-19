
extract.mstan4bartFit <-
  function(object,
           value = c("mu", "ppd", "fixef", "mu.fixef", "ranef", "mu.ranef", "bart", "sigma"),
           sample = c("train", "test"),
           combine_chains = TRUE,
           ...)
{
  value  <- match.arg(value)
  sample <- match.arg(sample)
  
  n_samples <- dim(object$bart_train)[2L]
  n_obs     <- dim(object$bart_train)[1L]
  n_chains  <- dim(object$bart_train)[3L]
  n_obs_test <- if (!is.null(object$bart_test)) dim(object$bart_test)[1L] else 0L
  n_fixef <- if (!is.null(object$fixef)) dim(object$fixef)[1L] else 0L
  n_warmup <- if (!is.null(object$warmup$bart_train)) dim(object$warmup$bart_train)[1L] else 0L
  n_bart_vars <- dim(object$bart_varcount)[1L]
  
  if (!is.null(object$ranef)) {
    n_ranef_levels <- length(object$ranef)
    n_ranef_at_level  <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[1L])
    n_groups_at_level <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[2L])
  } else {
    n_ranef_levels <- n_ranef_at_level <- n_groups_at_level <- 0L
  }
  
    
  mu.fixef <- mu.ranef <- mu.bart <- 0
  if (sample == "train") {
    if (n_fixef > 0L)        X <- object$X
    if (n_ranef_levels > 0L) Zt <- object$Zt
    n_obs_inf <- n_obs
  } else {
    if (n_fixef > 0L)        X <- object$X.test
    if (n_ranef_levels > 0L) Zt <- object$Zt.test
    n_obs_inf <- n_obs_test
  }
  
  mu.fixef <- mu.ranef <- 0
  if (value %in% c("mu", "ppd", "mu.fixef") && n_fixef > 0L) {
    # b_0 + b_1 * (x_1 - mean(x_1)) + ... = 
    # b_0 - b_1 * mean(x_1) - b_2 * mean(x_2) + 
    fixef_flat <- matrix(object$fixef, dim(object$fixef)[1L], prod(dim(object$fixef)[-1L]),
                         dimnames = list(predictor = dimnames(object$fixef)[[1L]], sample = NULL))
    
    keep_cols <- names(object$X_means) != "(Intercept)"
    intercept_delta <- apply(fixef_flat[keep_cols,,drop = FALSE] * object$X_means[keep_cols], 2L, sum)
    
    mu.fixef <- array(t(t(X %*% fixef_flat) - intercept_delta),
                      c(n_obs_inf, n_samples, n_chains),
                      dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  if (value %in% c("mu", "ppd", "mu.ranef") %% n_ranef_levels > 0L) {
    # put the dimensions of the random effects at the end so that they can be combined on
    # a per-sample, per-chain basis with the other random effects
    # then permute back to original
    b <- aperm(array(unlist(lapply(object$ranef, function(x) aperm(x, c(3L, 4L, 1L, 2L)))),
               c(n_samples, n_chains, sum(n_groups_at_level * n_ranef_at_level))), c(3L, 1L, 2L))
  
    # b_rows <- grep("^b\\.", rownames(object$chain_results[[1L]]$sample$stan$raw))
    # all(b[,,1] == object$chain_results[[1L]]$sample$stan$raw[b_rows,])
    # all(b[,,2] == object$chain_results[[2L]]$sample$stan$raw[b_rows,])
    
    b_mat <- matrix(b, dim(b)[1L], prod(dim(b)[-1L]))
    Zb <- Matrix::crossprod(Zt, b_mat) # getting a memory error on this
    # Zb <- crossprod(as.matrix(Zt), b_mat)
    mu.ranef <- array(Zb,
                      c(n_obs_inf, n_samples, n_chains),
                      dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  bart <- if (sample == "train") object$bart_train else object$bart_test
  
  if (value %in% "ppd") {
    eps <- array(rnorm(n_obs_inf * n_samples * n_chains, 0, rep(as.vector(object$sigma), each = n_obs_inf)),
                 c(n_obs_inf, n_samples, n_chains))
  }
  
  result <- switch(value,
                   mu       = bart + mu.fixef + mu.ranef,
                   ppd      = bart + mu.fixef + mu.ranef + eps,
                   mu.fixef = mu.fixef,
                   mu.ranef = mu.ranef,
                   bart     = bart,
                   varcount = object$bart_varcount,
                   ranef    = object$ranef,
                   fixef    = object$fixef,
                   Sigma    = object$Sigma,
                   sigma    = object$sigma)
  
  combine_chains_f <- function(x) {
    if (is.array(x)) {
      d <- dim(x)
      l_d <- length(d)
      t_d <- seq.int(l_d - 1L, l_d)
      dn <- dimnames(x)
      if (!is.null(dn)) {
        dn <- dn[seq.int(1L, l_d - 2L)]
        dn[l_d - 1L] <- list(NULL)
        names(dn)[l_d - 1L] <- "sample"
      }
      array(x, c(d[-t_d], prod(d[t_d])), dimnames = dn)
    } else if (is.matrix(x)) {
      as.vector(x)
    }
  }
  
  if (combine_chains) {
    if (is.list(result)) {
      result <- lapply(result, combine_chains_f)
  } else {
      result <- combine_chains_f(result)
    }
  }
  
  result
}



fitted.mstan4bartFit <-
  function(object,
           value = c("mu", "ppd", "fixef", "mu.fixef", "ranef", "mu.ranef", "bart", "sigma"),
           sample = c("train", "test"),
           ...)
{
  samples <- extract(object, value, sample, combine_chains = TRUE)
  
  average_samples_f <- function(x) {
    if (is.array(x)) {
      d <- dim(x)
      l_d <- length(d)
      apply(x, seq(1L, l_d - 1L), mean)
    } else if (is.matrix(x)) {
      apply(x, 1L, mean)
    } else {
      mean(x)
    }
  }
  
  if (is.list(samples)) {
    result <- lapply(samples, average_samples_f)
  } else {
    result <- average_samples_f(samples)
  }
  result
}

