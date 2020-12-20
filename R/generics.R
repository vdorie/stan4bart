combine_chains_f <- function(x) {
  if (is.array(x) && length(dim(x)) > 2L) {
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
  } else if (is.matrix(x) || (!is.null(dim(x)) && length(dim(x)) == 2L)) {
    as.vector(x)
  }
}

extract.mstan4bartFit <-
  function(object,
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef", "indiv.bart", "sigma"),
           sample = c("train", "test"),
           combine_chains = TRUE,
           sample_new_levels = TRUE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  type   <- match.arg(type)
  sample <- match.arg(sample)
  
  n_samples <- dim(object$bart_train)[2L]
  n_obs     <- dim(object$bart_train)[1L]
  n_chains  <- dim(object$bart_train)[3L]
  n_obs_test <- if (!is.null(object$bart_test)) dim(object$bart_test)[1L] else 0L
  n_fixef <- if (!is.null(object$fixef)) dim(object$fixef)[1L] else 0L
  n_warmup <- if (!is.null(object$warmup$bart_train)) dim(object$warmup$bart_train)[2L] else 0L
  n_bart_vars <- dim(object$bart_varcount)[1L]
  
  if (!is.null(object$ranef)) {
    n_ranef_levels <- length(object$ranef)
  } else {
    n_ranef_levels <- 0L
  }
  
    
  if (sample == "train") {
    if (n_fixef > 0L)        X <- object$X
    if (n_ranef_levels > 0L) reTrms <- object$reTrms
    n_obs_inf <- n_obs
  } else {
    if (n_fixef > 0L)        X <- object$test$X
    if (n_ranef_levels > 0L) reTrms <- object$test$reTrms
    n_obs_inf <- n_obs_test
  }
  
  indiv.fixef <- indiv.ranef <- indiv.bart <- 0
  if (type %in% c("ev", "ppd", "indiv.fixef") && n_fixef > 0L) {
            
    indiv.fixef <- fitted_fixed(object, X)
    
  }
  if (type %in% c("ev", "ppd", "indiv.ranef") %% n_ranef_levels > 0L) {
    
    indiv.ranef <- fitted_random(object, reTrms, sample_new_levels)
  }
  indiv.bart <- if (sample == "train") object$bart_train else object$bart_test
  
  if (type %in% "ppd") {
    eps <- array(rnorm(n_obs_inf * n_samples * n_chains, 0, rep(as.vector(object$sigma), each = n_obs_inf)),
                 c(n_obs_inf, n_samples, n_chains))
  }
  
  offset <- 0
  if (type %in% c("ev", "ppd")) {
    if (sample == "train" && !is.null(object$offset) && length(object$offset) > 0L) {
      offset <- object$offset
    } else if (sample == "test" && !is.null(object$test$offset) && length(object$offset) > 0L) {
      offset <- object$test$offset
    }
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef + offset,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef + offset + eps,
                   indiv.fixef = indiv.fixef,
                   indiv.ranef = indiv.ranef,
                   indiv.bart  = indiv.bart,
                   varcount    = object$bart_varcount,
                   ranef       = object$ranef,
                   fixef       = object$fixef,
                   Sigma       = object$Sigma,
                   sigma       = object$sigma)
    
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
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef", "indiv.bart", "sigma"),
           sample = c("train", "test"),
           sample_new_levels = TRUE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  samples <- extract(object, type, sample, combine_chains = TRUE, sample_new_levels = sample_new_levels)
  
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

fitted_fixed <- function(object, x)
{
  fixef <- extract(object, "fixef", combine_chains = TRUE)
  
  a <- colnames(x)
  b <- dimnames(fixef)$predictor
  if (any(a %not_in% b))
    stop("newdata missing old fixef names: '", paste0(a[a %not_in% b], collapse = "', '"), "'")
  if (any(b %not_in% a))
    stop("newdata has unrecognized fixef names: '", paste0(b[b %not_in% a], collapse = "', '"), "'")
  if (any(a != b)) {
    warning("newdata columns not in same order as training data")
    x <- x[,b,drop = FALSE]
  }
  
  # b_0 + b_1 * (x_1 - mean(x_1)) + ... = 
  # b_0 - b_1 * mean(x_1) - b_2 * mean(x_2) + 
  
  x_means <- object$X_means
  
  keep_cols <- names(x_means) != "(Intercept)"
  intercept_delta <- apply(fixef[keep_cols,,drop = FALSE] * x_means[keep_cols], 2L, sum)
  
  n_obs     <- nrow(x)
  n_samples <- dim(object$fixef)[2L]
  n_chains  <- dim(object$fixef)[3L]
  
  array(t(t(x %*% fixef) - intercept_delta),
        c(n_obs, n_samples, n_chains),
          dimnames = list(observation = NULL, sample = NULL, chain = NULL))
}

fitted_random <- function(object, reTrms, sample_new_levels)
{
  n_obs     <- ncol(reTrms$Zt)
  n_samples <- dim(object$fixef)[2L]
  n_chains  <- dim(object$fixef)[3L] 
  
  ns.re <- names(re <- object$ranef)
  nRnms <- names(Rcnms <- reTrms$cnms)
  if (!all(nRnms %in% ns.re)) 
    stop("grouping factors specified in re.form that were not present in original model")
  
  new_levels <- lapply(reTrms$flist, function(x) levels(factor(x)))
  
  re_x <- Map(function(r, n, Sigma) levelfun(r, n, sample_new_levels, Sigma), 
              re[names(new_levels)], new_levels, object$Sigma[names(new_levels)])
  get_re <- function(rname, cnms) {
    nms <- dimnames(re[[rname]])$predictor
    if (identical(cnms, "(Intercept)") && length(nms) == 1 && grepl("^s(.*)$", nms)) {
      cnms <- nms
    }
    miss_names <- setdiff(cnms, nms)
    if (length(miss_names) > 0) {
      stop("random effects specified in re.form that were not present in original model ",
           paste(miss_names, collapse = ", "))
    }
    re_x[[rname]][cnms,,,,drop = FALSE]
  }
  re_new <- Map(get_re, nRnms, Rcnms)
  
  n_groups_at_level.test <- sapply(re_new, function(x) dim(x)[1L])
  n_ranef_at_level.test  <- sapply(re_new, function(x) dim(x)[2L])
  
  # put the dimensions of the random effects at the end so that they can be combined on
  # a per-sample, per-chain basis with the other random effects
  # then permute back to original
  b <- aperm(array(unlist(lapply(re_new, function(x) aperm(x, c(3L, 4L, 1L, 2L)))),
             c(n_samples, n_chains, sum(n_groups_at_level.test * n_ranef_at_level.test))), c(3L, 1L, 2L))
  
  # b_rows <- grep("^b\\.", rownames(object$chain_results[[1L]]$sample$stan$raw))
  # all(b[,,1] == object$chain_results[[1L]]$sample$stan$raw[b_rows,])
  # all(b[,,2] == object$chain_results[[2L]]$sample$stan$raw[b_rows,])
  
  b_mat <- matrix(b, dim(b)[1L], prod(dim(b)[-1L]))
  Zb <- Matrix::crossprod(reTrms$Zt, b_mat) # getting a memory error on this
  
  array(Zb, c(n_obs, n_samples, n_chains),
        dimnames = list(observation = NULL, sample = NULL, chain = NULL))
}

predict.mstan4bartFit <-
  function(object, newdata, offset,
           type = c("ev", "ppd", "indiv.fixef", "indiv.ranef", "indiv.bart"),
           combine_chains = TRUE,
           sample_new_levels = TRUE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  type <- match.arg(type)
  
  if (length(combine_chains) > 1 || !is.logical(combine_chains) || is.na(combine_chains))
    stop("'combine_chains' must be TRUE or FALSE")
  
  mc <- match.call()

  if (is.null(mc[["newdata"]]) || is.null(newdata))
    return(extract(object, type = type, combine_chains = combine_chains))
  
  testFrames <- switch(type,
                       indiv.fixef = "fixed",
                       indiv.ranef = "random",
                       indiv.bart  = "bart",
                       "all")
  testData <- getTestDataFrames(object, newdata, type = testFrames)
  
  n_samples <- dim(object$bart_train)[2L]
  n_chains  <- dim(object$bart_train)[3L]
  n_obs     <- dim(testData$X.bart)[1L]
  n_fixef   <- if (!is.null(object$fixef)) dim(object$fixef)[1L] else 0L
  n_warmup  <- if (!is.null(object$warmup$bart_train)) dim(object$warmup$bart_train)[2L] else 0L
  n_bart_vars <- dim(object$bart_varcount)[1L]
  
  if (!is.null(object$ranef)) {
    n_ranef_levels <- length(object$ranef)
    n_ranef_at_level  <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[1L])
    n_groups_at_level <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[2L])
  } else {
    n_ranef_levels <- n_ranef_at_level <- n_groups_at_level <- 0L
  }

  indiv.fixef <- indiv.ranef <- indiv.bart <- 0
  if (type %in% c("ev", "ppd", "indiv.fixef")) {
    if (!is.null(testData$X)) {
      indiv.fixef <- fitted_fixed(object, testData$X)
    } else {
      if (type == "indiv.fixef")
        stop("predict called with type 'indiv.fixef', but model does not include fixed effect terms")
      indiv.fixef <- 0
    }
  }
  if (type %in% c("ev", "ppd", "indiv.bart")) {
    if (is.null(object$sampler.bart))
      stop("predict for bart components requires 'bart_args' to contain 'keepTrees' as 'TRUE'")
    indiv.bart <- .Call(C_stan4bart_predictBART, object$sampler.bart, testData$X.bart, NULL)
    dimnames(indiv.bart) <-  list(observation = NULL, sample = NULL, chain = NULL)
    for (i_chain in seq_len(n_chains)) {
      indiv.bart[,,i_chain] <- object$range.bart["min",i_chain] +
        (0.5 + indiv.bart[,,i_chain]) * (object$range.bart["max",i_chain] - object$range.bart["min",i_chain])
    }
  }
  if (type %in% c("ev", "ppd", "indiv.ranef")) {
    
    if (!is.null(testData$reTrms)) {
      indiv.ranef <- fitted_random(object, testData$reTrms, sample_new_levels)
    } else {
      if (type == "indiv.ranef")
        stop("predict called with type 'indiv.ranef', but model does not include random effect terms")
      indiv.ranef <- 0
    }
  }
  
  if (type %in% "ppd") {
    eps <- array(rnorm(n_obs * n_samples * n_chains, 0, rep(as.vector(object$sigma), each = n_obs)),
                 c(n_obs, n_samples, n_chains))
  }
  
  if (type %in% c("ev", "ppd") && is.null(mc$offset)) {
    offset <- 0
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef + offset,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef + offset + eps,
                   indiv.fixef = indiv.fixef,
                   indiv.ranef = indiv.ranef,
                   indiv.bart  = indiv.bart)
  
  if (combine_chains) {
    if (is.list(result)) {
      result <- lapply(result, combine_chains_f)
  } else {
      result <- combine_chains_f(result)
    }
  }
  
  result
}
