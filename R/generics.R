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


extract.mstan4bartFit <-
  function(object,
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef", "indiv.bart", "sigma"),
           sample = c("train", "test"),
           combine_chains = TRUE,
           ...)
{
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
    n_ranef_at_level  <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[1L])
    n_groups_at_level <- sapply(object$ranef, function(ranef_i) dim(ranef_i)[2L])
  } else {
    n_ranef_levels <- n_ranef_at_level <- n_groups_at_level <- 0L
  }
  
    
  if (sample == "train") {
    if (n_fixef > 0L)        X <- object$X
    if (n_ranef_levels > 0L) Zt <- object$reTrms$Zt
    n_obs_inf <- n_obs
  } else {
    if (n_fixef > 0L)        X <- object$X.test
    if (n_ranef_levels > 0L) Zt <- object$reTrms$Zt.test
    n_obs_inf <- n_obs_test
  }
  
  indiv.fixef <- indiv.ranef <- indiv.bart <- 0
  if (type %in% c("ev", "ppd", "indiv.fixef") && n_fixef > 0L) {
    # b_0 + b_1 * (x_1 - mean(x_1)) + ... = 
    # b_0 - b_1 * mean(x_1) - b_2 * mean(x_2) + 
    fixef_flat <- matrix(object$fixef, dim(object$fixef)[1L], prod(dim(object$fixef)[-1L]),
                         dimnames = list(predictor = dimnames(object$fixef)[[1L]], sample = NULL))
    
    keep_cols <- names(object$X_means) != "(Intercept)"
    intercept_delta <- apply(fixef_flat[keep_cols,,drop = FALSE] * object$X_means[keep_cols], 2L, sum)
    
    indiv.fixef <- array(t(t(X %*% fixef_flat) - intercept_delta),
                         c(n_obs_inf, n_samples, n_chains),
                         dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  if (type %in% c("ev", "ppd", "indiv.ranef") %% n_ranef_levels > 0L) {
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
    indiv.ranef <- array(Zb,
                         c(n_obs_inf, n_samples, n_chains),
                         dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  indiv.bart <- if (sample == "train") object$bart_train else object$bart_test
  
  if (type %in% "ppd") {
    eps <- array(rnorm(n_obs_inf * n_samples * n_chains, 0, rep(as.vector(object$sigma), each = n_obs_inf)),
                 c(n_obs_inf, n_samples, n_chains))
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef + eps,
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
           ...)
{
  samples <- extract(object, type, sample, combine_chains = TRUE)
  
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

predict.mstan4bartFit <-
  function(object, newdata, offset,
           type = c("ev", "ppd", "indiv.fixef", "indiv.ranef", "indiv.bart"),
           combine_chains = TRUE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  type <- match.arg(type)
  
  if (length(combine_chains) > 1 || !is.logical(combine_chains) || is.na(combine_chains))
    stop("'combine_chains' must be TRUE or FALSE")
  
  mc <- match.call()

  if (is.null(mc[["newdata"]]) || is.null(newdata))
    return(extract(object, type = type, combine_chains = combine_chains))
  
  glcall <- mc
  glcall[[1L]] <- quoteInNamespace(glFormula)
  glcall$object <- NULL
  formula <- object$formula
  response_var <- formula[[2L]]
  formula[[2L]] <- NULL
  glcall[["formula"]] <- formula
  
  newdata[[response_var]] <- NULL
  
  names(glcall)[names(glcall) == "newdata"] <- "data"
  glcall[["data"]] <- newdata
  
  if ("type" %in% names(glcall)) glcall$type <- NULL
  if ("combine_chains" %in% names(glcall)) glcall$combine_chains <- NULL
  
  if (!is.null(object$treatment))
    glcall$treatment <- object$treatment
  
  glcall$na.action <- object$na.action
  glcall$control <- glmerControl(check.nobs.vs.nlev = "ignore",
                                 check.nlev.gtr.1   = "ignore",
                                 check.nobs.vs.nRE  = "ignore",
                                 check.rankX        = "ignore",
                                 check.scaleX       = "ignore",
                                 check.formula.LHS  = "ignore",
                                 check.response.not.const = "ignore")
  
  glmod <- eval(glcall, parent.frame())
  
  n_samples <- dim(object$bart_train)[2L]
  n_chains  <- dim(object$bart_train)[3L]
  n_obs     <- dim(glmod$bartData@x)[1L]
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
    fixef <- extract(object, "fixef", combine_chains = TRUE)
    
    a <- colnames(glmod$X)
    b <- dimnames(fixef)$predictor
    if (any(a %not_in% b))
      stop("newdata missing old fixef names: '", paste0(a[a %not_in% b], collapse = "', '"), "'")
    if (any(b %not_in% a))
      stop("newdata has unrecognized fixef names: '", paste0(b[b %not_in% a], collapse = "', '"), "'")
    glmod$X <- glmod$X[,b,drop = FALSE]
    rm(a, b)
    
    fixef <- extract(object, "fixef", combine_chains = TRUE)
    
    keep_cols <- names(object$X_means) != "(Intercept)"
    intercept_delta <- apply(fixef[keep_cols,,drop = FALSE] * object$X_means[keep_cols], 2L, sum)
    
    indiv.fixef <- array(t(t(glmod$X %*% fixef) - intercept_delta),
                         c(n_obs, n_samples, n_chains),
                         dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  if (type %in% c("ev", "ppd", "indiv.bart")) {
    if (is.null(object$sampler.bart))
      stop("predict for bart components requires 'bart_args' to contain 'keepTrees' as 'TRUE'")
    indiv.bart <- .Call(C_stan4bart_predictBART, object$sampler.bart, glmod$bartData@x, glmod$bartData@offset)
    dimnames(indiv.bart) <-  list(observation = NULL, sample = NULL, chain = NULL)
  }
  if (type %in% c("ev", "ppd", "indiv.ranef")) {
    stop("predict for random effects not yet implemented")
    re.form <- reOnly(formula(object))
    newRE <- get_test_reterms(object, newdata, re.form, allow.new.levels = TRUE)
    REvals <- base::drop(as(newRE$b %*% newRE$Zt, "matrix"))
    pred <- pred + REvals
     
    re_terms.test <- get_test_reterms(object, newdata)
  }
  
  if (type %in% "ppd") {
    eps <- array(rnorm(n_obs * n_samples * n_chains, 0, rep(as.vector(object$sigma), each = n_obs)),
                 c(n_obs, n_samples, n_chains))
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef + eps,
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

get_test_reterms <- function(object, newdata, re.form = NULL, na.action = na.pass, 
    allow.new.levels = FALSE, sparse = max(lengths(orig.random.levs)) > 100) 
{
  fixed.na.action <- NULL
  if (!identical(na.action, na.pass)) {
    fixed.terms <- terms(object, type = "fixed")
    mfnew <- model.frame(delete.response(terms(object, type = "fixed")), newdata, na.action = na.action)
    fixed.na.action <- attr(mfnew, "na.action")
  }
  newdata.NA <- newdata
  if (!is.null(fixed.na.action)) {
    newdata.NA <- newdata.NA[-fixed.na.action,]
  }
  tt <- delete.response(terms(object, type = "random"))
  orig.random.levs <- get.orig.levs(object, type = "random", newdata = newdata.NA)
  sparse_eval <- suppressWarnings(sparse)
  orig.random.cntr <- get.orig.levs(object, type = "random", FUN = contrasts, sparse = sparse_eval)
  if (inherits(re.form, "formula")) {
    pv <- attr(tt, "predvars")
    for (i in 2:(length(pv))) {
      missvars <- setdiff(all.vars(pv[[i]]), all.vars(re.form))
      for (mv in missvars) {
        newdata.NA[[mv]] <- NA
      }
    }
  }
  rfd <- suppressWarnings(model.frame(tt, newdata.NA, na.action = na.pass, 
                                      xlev = orig.random.levs))
  termvars <- unique(unlist(lapply(findbars(formula(object, type = "random")), function(x) all.vars(x[[2]]))))
  for (fn in Reduce(intersect, list(names(orig.random.cntr), termvars, names(rfd)))) {
    if (!is.factor(rfd[[fn]])) 
      rfd[[fn]] <- factor(rfd[[fn]])
      contrasts(rfd[[fn]]) <- orig.random.cntr[[fn]]
  }
  if (!is.null(fixed.na.action)) 
    attr(rfd, "na.action") <- fixed.na.action
  
  
  re_new <- NULL
  if (inherits(re.form, "formula")) {
    if (length(fixed.na.action) > 0) {
      newdata <- newdata[-fixed.na.action, ]
    }
    ReTrms <- mkReTrms(findbars(re.form[[2]]), rfd)
    ReTrms <- within(ReTrms, Lambdat@x <- unname(object$reTrms$theta[Lind]))
    if (!allow.new.levels && any(vapply(ReTrms$flist, anyNA, NA))) 
      stop("NAs are not allowed in prediction data", " for grouping variables unless allow.new.levels is TRUE")
    ns.re <- names(re <- object$ranef)
    nRnms <- names(Rcnms <- ReTrms$cnms)
    if (!all(nRnms %in% ns.re)) 
      stop("grouping factors specified in re.form that were not present in original model")
    new_levels <- lapply(ReTrms$flist, function(x) levels(factor(x)))
    # TODO: rewrite levelfun to use samples for levels and not set to 0
    re_x <- Map(function(r, n) levelfun(r, n, allow.new.levels = allow.new.levels), 
                re[names(new_levels)], new_levels)
    re_new <- lapply(seq_along(nRnms), function(i) {
      rname <- nRnms[i]
      if (!all(Rcnms[[i]] %in% dimnames(re[[rname]])[["predictor"]])) 
        stop("random effects specified in re.form that were not present in original model")
      re_x[[rname]][Rcnms[[i]],,,,drop = FALSE]
    })
  }
  Zt <- ReTrms$Zt
  attr(Zt, "na.action") <- attr(re_new, "na.action") <- fixed.na.action
  list(Zt = Zt, b = re_new, Lambdat = ReTrms$Lambdat)
}


