combine_chains_f <- function(x) {
  if (is.array(x) && length(dim(x)) > 2L) {
    d <- dim(x)
    l_d <- length(d)
    t_d <- seq.int(l_d - 1L, l_d)
    dn <- dimnames(x)
    if (!is.null(dn)) {
      dn <- dn[seq.int(1L, l_d - 2L)]
      dn[l_d - 1L] <- list(NULL)
      names(dn)[l_d - 1L] <- "iterations:chains"
    }
    array(x, c(d[-t_d], prod(d[t_d])), dimnames = dn)
  } else if (is.matrix(x) || (!is.null(dim(x)) && length(dim(x)) == 2L)) {
    as.vector(x)
  }
}

as.array.mstan4bartFit <- function (x, include_warmup = FALSE, ...) 
{
  include_warmup_orig <- include_warmup
  if (is.character(include_warmup)) {
    if (length(include_warmup) != 1L || include_warmup != "only")
      stop("'include_warmup' must be logical or \"only\"")
    include_warmup <- TRUE
    only_warmup <- TRUE
  } else if (!is.logical(include_warmup) || length(include_warmup) != 1L || is.na(include_warmup)) {
    stop("'include_warmup' must be logical or \"only\"")
  } else {
    only_warmup <- FALSE
  }

  result <- get_samples(x$stan[grep("^(?:gamma|beta|b|aux)\\.", dimnames(x$stan)[[1L]], perl = TRUE),,,drop = FALSE], include_warmup, only_warmup)
  
  par_names <- dimnames(result)[[1L]]
  if ("gamma.1" %in% par_names) {
    par_names[match("gamma.1", par_names)] <- "(Intercept)"
  }
  if (any(startsWith(par_names, "beta."))) {
    beta_names <- par_names[startsWith(par_names, "beta")]
    beta_col <- as.integer(sub("beta\\.", "", beta_names))
    has_intercept <- "(Intercept)" %in% colnames(x$X)
    
    par_names[startsWith(par_names, "beta")] <- colnames(x$X)[beta_col + as.integer(has_intercept)]
  }
  if (any(startsWith(par_names, "b."))) {
    numGroupingFactors <- length(x$reTrms$cnms)
    numRanefPerGroupingFactor <- unname(lengths(x$reTrms$cnms))
    group_names <- lapply(x$reTrms$flist, levels)
    
    ranef_names <- paste0(
      "b[",
      unlist(sapply(seq_along(x$reTrms$cnms), function(j) rep(x$reTrms$cnms[[j]], times = length(group_names[[j]])))),
      " ",
      rep(names(group_names), times = numRanefPerGroupingFactor * lengths(group_names)),
      ":",
      unlist(sapply(seq_along(group_names), function(j) rep(group_names[[j]], each = numRanefPerGroupingFactor[j]))),
      "]"
    )
    
    par_names[startsWith(par_names, "b.")] <- ranef_names
  }
  if ("aux.1" %in% par_names) {
    par_names[match("aux.1", par_names)] <- "sigma"
  }
  
  dimnames(result)[[1L]] <- par_names
  
  # TODO: modify this to avoid using extract and pull directly from theta_L
   # very low priority
  if (any(startsWith(par_names, "theta_L."))) {
    Sigmas <- extract(x, "Sigma", include_warmup = include_warmup_orig, combine_chains = FALSE)
    if (length(Sigmas) > 0L) {
      Sigmas.list <- lapply(Sigmas, function(Sigma) {
        r <- apply(Sigma, c(3L, 4L), function(x) x[upper.tri(x, diag = TRUE)])
        if (length(dim(r)) == 2L) r <- array(r, c(1L, dim(r)))
        
        dn <- outer(dimnames(Sigma)[[1]], dimnames(Sigma)[[2]], FUN = paste, sep = ",")
        dimnames(r)[[1L]] <- dn[upper.tri(dn, diag = TRUE)]
        aperm(r, c(2L, 3L, 1L))
      })
      n_samples <- dim(Sigmas.list[[1L]])[1L]
      n_chains  <- dim(Sigmas.list[[1L]])[2L]
      n_pars    <- sapply(Sigmas.list, function(x) dim(x)[3L])
      
      Sigma_names <- paste0("Sigma[",
                            rep(names(Sigmas), times = n_pars),
                          ":",
                            unlist(lapply(Sigmas.list, function(x) dimnames(x)[[3L]])),
                            "]")
      Sigmas.arr <- array(unlist(Sigmas.list),
                          dim = c(n_samples, n_chains, sum(n_pars)),
                          dimnames = list(iterations = NULL, chains = dimnames(Sigmas[[1L]])[[3L]],
                                          Sigma_names))
      result <- aperm(result, c(2L, 3L, 1L))
      result <- array(c(result, Sigmas.arr),
                      dim = c(n_samples, n_chains, dim(result)[3L] + dim(Sigmas.arr)[3L]),
                      dimnames = list(iterations = dimnames(result)[[1]],
                                      chains = dimnames(result)[[2]],
                                      parameters = c(dimnames(result)[[3L]], dimnames(Sigmas.arr)[[3L]])))
      
      result <- aperm(result, c(3L, 1L, 2L))
    }
  }
  aperm(result, c(2L, 3L, 1L))
}

as.matrix.mstan4bartFit <- function (x, ...) 
{
  result <- as.array(x, ...)
  
  as.matrix(result, ncol = dim(result)[3L])
}

array_bind_samples <- function(x, y) {
  if (is.null(dim(x)) || length(dim(x)) == 1L) {
    c(x, y)
  } else if (length(dim(x)) == 3L) {
    aperm(array(c(aperm(x, c(1L, 3L, 2L)),
                  aperm(y, c(1L, 3L, 2L))),
                dim = c(dim(x)[1L], dim(x)[3L], dim(x)[2L] + dim(y)[2L]),
                dimnames = list(parameters = dimnames(x)[[1L]],
                chain = dimnames(x)[[3L]],              
                iterations = NULL)),
          c(1L, 3L, 2L))
  } else if (length(dim(x)) == 2L) {
    t(matrix(c(t(x), t(y)),
             nrow = ncol(x), ncol = nrow(x) + nrow(y),
             dimnames = list(chain = colnames(x), iterations = NULL)))
  } else {
    stop("unrecognized dimensions of input x: ", paste0(dim(x), collapse = " x "))
  }
}

get_samples <- function(expr, include_warmup, only_warmup)
{
  get_warmup_expression <- function(expr) {
    if (length(expr) == 3L) {
      if ((expr[[1L]] == "$" || expr[[1L]] == "[[") && expr[[2L]] == "object") {
        lhs <- expr[[2L]] # it can be x$stan or object$stan, depending on context
        expr[[2L]] <- quote(object$warmup)
        expr[[2L]][[2L]] <- lhs
        expr[[3L]] <- get_warmup_expression(expr[[3L]])
      } else {
        expr[[2L]] <- get_warmup_expression(expr[[2L]])
        expr[[3L]] <- get_warmup_expression(expr[[3L]])
      }
    } else if (length(expr) > 1L) {
      for (i in seq.int(2L, length(expr)))
        expr[[i]] <- get_warmup_expression(expr[[i]])
    }
    expr
  }
  mc <- match.call()
  expr_sample <- mc$expr
  expr_warmup <- get_warmup_expression(expr_sample)
  
  if (!include_warmup) {
    result <- eval(expr_sample, envir = parent.frame())
  } else if (only_warmup) {
    result <- eval(expr_warmup, envir = parent.frame())
  } else {
    warmup <- eval(expr_warmup, envir = parent.frame())
    sample <- eval(expr_sample, envir = parent.frame())
    result <- array_bind_samples(warmup, sample)
  }
  result
}

extract.mstan4bartFit <-
  function(object,
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
                    "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan"),
           sample = c("train", "test"),
           combine_chains = TRUE,
           sample_new_levels = TRUE,
           include_warmup = FALSE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")
  
  type   <- match.arg(type)
  sample <- match.arg(sample)
  
  include_warmup_orig <- include_warmup
  if (is.character(include_warmup)) {
    if (length(include_warmup) != 1L || include_warmup != "only")
      stop("'include_warmup' must be logical or \"only\"")
    include_warmup <- TRUE
    only_warmup <- TRUE
  } else if (!is.logical(include_warmup) || length(include_warmup) != 1L || is.na(include_warmup)) {
    stop("'include_warmup' must be logical or \"only\"")
  } else {
    only_warmup <- FALSE
  }
  
  is_bernoulli <- object$family$family == "binomial"
  if (type == "sigma" && is_bernoulli)
    stop("cannot extract 'sigma': binary outcome model does not have a residual standard error parameter")
  if (type == "k" && is.null(object$k))
    stop("cannot extract 'k': model was not fit with end-node sensitivity as a modeled parameter")
  
  fixef_parameters <- grep("^beta|gamma", dimnames(object$stan)[[1L]])
  ranef_parameters <- startsWith(dimnames(object$stan)[[1L]], "b.")
  Sigma_parameters <- startsWith(dimnames(object$stan)[[1L]], "theta_L.")
  sigma_parameters <- dimnames(object$stan)[[1L]] %in% "aux.1"
  
  if (type == "fixef") {
    if (!any(fixef_parameters))
      stop("cannot extract fixef for model with no unmodeled parameters")
    result <- get_samples(object$stan[fixef_parameters,,,drop = FALSE],
                           include_warmup, only_warmup)
    dimnames(result)[[1L]] <- colnames(object$X)
    names(dimnames(result))[1L] <- "predictor"
  } else if (type == "ranef") {
    if (!any(ranef_parameters))
      stop("cannot extract ranef for model with no modeled parameters")
    
    ranef <- get_samples(object$stan[ranef_parameters,,,drop = FALSE],
                         include_warmup, only_warmup)
    
    numGroupingFactors <- length(object$reTrms$cnms)
    numRanefPerGroupingFactor <- unname(lengths(object$reTrms$cnms))
  
    result <- lapply(seq_len(numGroupingFactors), function(j) {
      ranef.group <- ranef[seq.int(object$reTrms$Gp[j] + 1L, object$reTrms$Gp[j + 1L]),,,drop = FALSE]
      result <- array(ranef.group, c(numRanefPerGroupingFactor[j], dim(ranef.group)[1L] %/% numRanefPerGroupingFactor[j], dim(ranef.group)[-1L]),
                      dimnames = list(
                        predictor = object$reTrms$cnms[[j]],
                        group = levels(object$reTrms$flist[[j]]),
                        iterations = NULL,
                        chains = dimnames(ranef.group)[[3L]]))
      result
    })
    names(result) <- names(object$reTrms$cnms)
  } else if (type == "Sigma") {
    if (!any(Sigma_parameters))
      stop("cannot extract Sigma for model with no modeled parameters")
    
    thetas <- get_samples(object$stan[Sigma_parameters,,,drop = FALSE],
                          include_warmup, only_warmup)
    numRanefPerGroupingFactor <- unname(lengths(object$reTrms$cnms))
    nms <- names(object$reTrms$cnms)
    Sigma <- apply(thetas, c(2L, 3L), function(theta) mkVarCorr(sc = 1, object$reTrms$cnms, numRanefPerGroupingFactor, theta, nms))
    result <- lapply(seq_along(Sigma[[1L]]), function(j) {
      raw <- sapply(Sigma, function(Sigma.i) Sigma.i[[j]])
      array(raw, c(NROW(Sigma[[1L]][[j]]), NCOL(Sigma[[1L]][[j]]), NROW(Sigma), NCOL(Sigma)),
            dimnames = list(rownames(Sigma[[1L]][[j]]), colnames(Sigma[[1L]][[j]]),
                            iterations = NULL, chains = dimnames(thetas)[[3L]]))
    })
    names(result) <- names(Sigma[[1L]])
  } else if (type == "sigma") {
    # test that this exists is above
    result <- get_samples(object$stan[sigma_parameters,,,drop = FALSE],
                          include_warmup, only_warmup)
    result <- matrix(result, dim(result)[2L], dim(result)[3L], dimnames = dimnames(result)[2L:3L])
  } else if (type %in% c("bart_varcount", "k", "stan")) {
    result <- get_samples(object[[type]],
                          include_warmup, only_warmup)
  }
  # return parametric components early, as they don't require building model matrices
  if (type %in% c("fixef", "ranef", "sigma", "Sigma", "k", "varcount", "stan")) {
    if (combine_chains) {
      if (is.list(result)) {
        result <- lapply(result, combine_chains_f)
      } else {
        result <- combine_chains_f(result)
      }
    }
    
    return(result)
  }
  
  n_samples <- dim(object$bart_train)[2L]
  n_obs     <- 0L
  n_chains  <- dim(object$bart_train)[3L]
  n_fixef <- sum(fixef_parameters)
  n_warmup <- if (!is.null(object$warmup$bart_train)) dim(object$warmup$bart_train)[2L] else 0L
  n_bart_vars <- dim(object$bart_varcount)[1L]
  
  if (!is.null(object$reTrms) && length(object$reTrms) > 0L) {
    n_ranef_levels <- length(object$reTrms$cnms)
  } else {
    n_ranef_levels <- 0L
  }
  
  offset <- NULL
  # The data for the fitted model will always be subsetted to only the
  # complet cases, across all three possible data frames/matrices.
  #
  # However, with na.action == na.exclude, it is reasonable to
  # return predictions for sub-components of the model when they
  # have additional data.
  
  if (sample == "train") {
    na.action.fixed  <- attr(object$frame, "na.action.fixed")
    na.action.bart   <- attr(object$frame, "na.action.bart")
    na.action.random <- attr(object$frame, "na.action.random")
    na.action.all    <- attr(object$frame, "na.action.all")
    n_all <- nrow(object$frame)
        
    if (!is.null(object$offset) && length(object$offset) > 0L)
      offset <- object$offset
    
    if (n_fixef > 0L) {
      if (!inherits(na.action.fixed, "exclude")) {
        X <- object$X
      } else {
        frame <- model.frame(object, type = "fixed")[-na.action.fixed,,drop = FALSE]
        attr(frame, "na.action") <- na.pass
        X <- model.matrix(formula(object, type = "fixed"), frame)
        intercept_col <- match("(Intercept)", colnames(X))
        if (!is.na(intercept_col))
          X <- X[, -intercept_col, drop = FALSE]
      }
    }
    if (n_ranef_levels > 0L) {
      if (!inherits(na.action.random, "exclude")) {
        reTrms <- object$reTrms
      } else {
        frame <- model.frame(object, type = "random")[-na.action.random,,drop = FALSE]
        attr(frame, "na.action") <- na.pass
        reTrms <- mkReTrms(findbars(RHSForm(formula(object))), frame)
      }
    }
    if (inherits(na.action.bart, "exclude")) {
      frame <- model.frame(object, type = "bart")[-na.action.bart,,drop = FALSE]
      if (dim(object$bart_train)[1L] != nrow(frame) && type == "indiv.bart")
        warning("cannot obtain training predictions for rows where BART component has no NAs but other components do; add to test data instead")
      # The bart training result will always be of size N - na.action.all, even if
      # the bart training frame can have more rows.
      na.action.bart <- na.action.all
    }
  } else {
    na.action.fixed  <- object$test$na.action.fixed
    na.action.bart   <- object$test$na.action.bart
    na.action.random <- object$test$na.action.random
    na.action.all <- object$test.na.action.all
    n_all <- nrow(object$test$frame)
    
    if (n_fixef > 0L)        X <- object$test$X
    if (n_ranef_levels > 0L) reTrms <- object$test$reTrms
    if (!is.null(object$test$offset) && length(object$test$offset) > 0L)
      offset <- object$test$offset
  }
  
  offset_type <- object$offset_type
  
  indiv.fixef <- indiv.ranef <- indiv.bart <- 0
  if (type %in% c("ev", "ppd", "indiv.fixef") && n_fixef > 0L) {
    if (!is.null(offset) && offset_type %in% c("fixed", "parametric") && type != "indiv.fixef")
      indiv.fixef <- offset
    else
      indiv.fixef <- fitted_fixed(object, X, include_warmup_orig)
    
    if (inherits(na.action.fixed, "exclude")) {
      indiv.fixef.all <- array(NA_real_, c(n_all, dim(indiv.fixef)[-1L]), dimnames(indiv.fixef))
      indiv.fixef.all[-na.action.fixed,,] <- indiv.fixef
      indiv.fixef <- indiv.fixef.all
    }
  }
  if (type %in% c("ev", "ppd", "indiv.ranef") && n_ranef_levels > 0L) {
    if (!is.null(offset) && offset_type %in% c("random", "parametric") && type != "indiv.ranef") {
      if (offset_type == "parametric")
        indiv.ranef <- 0
      else
        indiv.ranef <- offset
    }
    else
      indiv.ranef <- fitted_random(object, reTrms, include_warmup_orig, sample_new_levels)
    
    if (inherits(na.action.random, "exclude")) {
      indiv.ranef.all <- array(NA_real_, c(n_all, dim(indiv.ranef)[-1L]), dimnames(indiv.ranef))
      indiv.ranef.all[-na.action.random,,] <- indiv.ranef
      indiv.ranef <- indiv.ranef.all
    }
  }
  
  if (type %in% c("ev", "ppd", "indiv.bart")) {
    if (!is.null(offset) && offset_type %in% "bart" && type != "indiv.bart")
      indiv.bart <- offset
    else
      indiv.bart <- if (sample == "train")
        get_samples(object$bart_train, include_warmup, only_warmup)
      else
        get_samples(object$bart_test, include_warmup, only_warmup)
    
    if (inherits(na.action.bart, "exclude")) {
      indiv.bart.all <- array(NA_real_, c(n_all, dim(indiv.bart)[-1L]), dimnames(indiv.bart))
      indiv.bart.all[-na.action.bart,,] <- indiv.bart
      indiv.bart <- indiv.bart.all
    }
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef,
                   indiv.fixef = indiv.fixef,
                   indiv.ranef = indiv.ranef,
                   indiv.bart  = indiv.bart,
                   varcount    = object$bart_varcount,
                   ranef       = object$ranef,
                   fixef       = object$fixef,
                   Sigma       = object$Sigma,
                   sigma       = object$sigma,
                   k           = object$k)
  
  if (type %in% c("ev", "ppd") && !is.null(offset) && offset_type == "default")
    result <- result + offset
  
  if (type %in% c("ev", "ppd") && is_bernoulli)
    result <- pnorm(result)
  if (type %in% "ppd") {
    if (is_bernoulli) {
      result <- array(rbinom(length(result), 1L, result), dim(result), dimnames = dimnames(result))
    } else {
      sigma <- extract(object, "sigma", combine_chains = TRUE)
      result <- result + 
        array(rnorm(dim(result)[1L] * n_samples * n_chains, 0, rep(as.vector(sigma),
                                                                   each = dim(result)[1L])),
             c(dim(result)[1L], n_samples, n_chains))
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
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
                    "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan"),
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

fitted_fixed <- function(object, x, include_warmup)
{
  fixef <- extract(object, "fixef", include_warmup = include_warmup, combine_chains = TRUE)
  
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
  # If there is an intercept, center_x sweeps out the column means and bundles it into
  # that term. If there is not, x is used un-centered.
  #intercept_delta <- 
  #  if (!all(keep_cols))
  #    apply(fixef[keep_cols,,drop = FALSE] * x_means[keep_cols], 2L, sum)
  #  else
  #    0
  
  # stan4bart will always have an intercept, although it should never be present in the x matrix
  intercept_delta <- apply(fixef[keep_cols,,drop = FALSE] * x_means[keep_cols], 2L, sum)
  
  n_obs     <- nrow(x)
  # n_warmup  <- dim(object$warmup$bart_train)[2L]
  # n_samples <- dim(object$bart_train)[2L]
  n_chains  <- dim(object$bart_train)[3L]
  n_samples <- ncol(fixef) %/% n_chains
  
  array(t(t(x %*% fixef) - intercept_delta),
        c(n_obs, n_samples, n_chains),
          dimnames = list(observation = NULL, sample = NULL, chain = NULL))
}

fitted_random <- function(object, reTrms, include_warmup, sample_new_levels)
{
  n_obs <- ncol(reTrms$Zt)
  
  ns.re <- names(re <- extract(object, "ranef", include_warmup = include_warmup, combine_chains = FALSE))
  n_samples <- dim(re[[1L]])[3L]
  n_chains  <- dim(re[[1L]])[4L]
  
  nRnms <- names(Rcnms <- reTrms$cnms)
  if (!all(nRnms %in% ns.re)) 
    stop("grouping factors specified in re.form that were not present in original model")
  
  new_levels <- lapply(reTrms$flist, function(x) levels(factor(x)))
  
  Sigmas <- extract(object, "Sigma", combine_chains = FALSE)
  re_x <- Map(function(r, n, Sigma) levelfun(r, n, sample_new_levels, Sigma), 
              re[names(new_levels)], new_levels, Sigmas[names(new_levels)])
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
        dimnames = list(observation = NULL, iterations = NULL, chain = dimnames(re_new[[1L]])$chain))
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
  n_fixef <- sum(grepl("^beta|gamma", dimnames(object$stan)[[1L]]))
  n_warmup  <- if (!is.null(object$warmup$bart_train)) dim(object$warmup$bart_train)[2L] else 0L
  n_bart_vars <- dim(object$bart_varcount)[1L]
  is_bernoulli <- object$family$family == "binomial"
  
  if (!is.null(object$reTrms) && length(object$reTrms) > 0L) {
    n_ranef_levels <- length(object$reTrms$cnms)
  } else {
    n_ranef_levels <- 0L
  }
  
  indiv.fixef <- indiv.ranef <- indiv.bart <- 0
  if (type %in% c("ev", "ppd", "indiv.fixef")) {
    if (!is.null(testData$X)) {
      indiv.fixef <- fitted_fixed(object, testData$X, FALSE)
    } else {
      if (type == "indiv.fixef")
        stop("predict called with type 'indiv.fixef', but model does not include fixed effect terms")
      indiv.fixef <- 0
    }
  }
  if (type %in% c("ev", "ppd", "indiv.bart")) {
    if (is.null(object$sampler.bart))
      stop("predict for bart components requires 'bart_args' to contain 'keepTrees' as 'TRUE'")
    indiv.bart <- .Call("stan4bart_predictBART", object$sampler.bart, testData$X.bart, NULL, PACKAGE = "stan4bart")
    dimnames(indiv.bart) <-  list(observation = NULL, sample = NULL, chain = NULL)
    if (!is_bernoulli) for (i_chain in seq_len(n_chains)) {
      indiv.bart[,,i_chain] <- object$range.bart["min",i_chain] +
        (0.5 + indiv.bart[,,i_chain]) * (object$range.bart["max",i_chain] - object$range.bart["min",i_chain])
    }
  }
  if (type %in% c("ev", "ppd", "indiv.ranef")) {
    
    if (!is.null(testData$reTrms) && length(testData$reTrms) > 0L) {
      indiv.ranef <- fitted_random(object, testData$reTrms, FALSE, sample_new_levels)
    } else {
      if (type == "indiv.ranef")
        stop("predict called with type 'indiv.ranef', but model does not include random effect terms")
      indiv.ranef <- 0
    }
  }
  
  
  if (type %in% c("ev", "ppd") && is.null(mc$offset)) {
    offset <- 0
  }
  
  result <- switch(type,
                   ev          = indiv.bart + indiv.fixef + indiv.ranef + offset,
                   ppd         = indiv.bart + indiv.fixef + indiv.ranef + offset,
                   indiv.fixef = indiv.fixef,
                   indiv.ranef = indiv.ranef,
                   indiv.bart  = indiv.bart)
  
  if (type %in% c("ev", "ppd") && is_bernoulli)
    result <- pnorm(result)
  if (type %in% "ppd") {
    if (is_bernoulli) {
      result <- array(rbinom(length(result), 1L, result), dim(result), dimnames = dimnames(result))
    } else {
      sigma <- extract(object, "sigma", combine_chains = TRUE)
      result <- result + 
        array(rnorm(n_obs * n_samples * n_chains, 0, rep(as.vector(sigma), each = n_obs)),
             c(n_obs, n_samples, n_chains))
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

