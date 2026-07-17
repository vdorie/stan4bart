combine_chains_f <- function(x) {
  if (is.array(x) && length(dim(x)) > 2L) {
    d <- dim(x)
    l_d <- length(d)
    dn <- dimnames(x)
    if (!is.null(dn)) {
      dn <- dn[seq.int(1L, l_d - 2L)]
      dn[l_d - 1L] <- list(NULL)
      names(dn)[l_d - 1L] <- "iterations:chains"
    }
    # merging the last two dims (iterations, chains) into one never reorders
    # the underlying column-major data, so this is a pure reshape: dim<-/
    # dimnames<- attach the new shape in place instead of array()'s copy.
    dim(x) <- c(d[seq.int(1L, l_d - 2L)], prod(d[seq.int(l_d - 1L, l_d)]))
    dimnames(x) <- dn
    x
  } else if (is.matrix(x) || (!is.null(dim(x)) && length(dim(x)) == 2L)) {
    as.vector(x)
  }
}

# Guard the include_warmup accessors against a fit that did not store warmup
# (the save_warmup = FALSE default): error informatively rather than let the
# NULL $warmup surface as a cryptic subsetting failure.
check_warmup_stored <- function(object, include_warmup) {
  wants_warmup <- isTRUE(include_warmup) ||
    (is.character(include_warmup) && length(include_warmup) == 1L && include_warmup == "only")
  if (wants_warmup && is.null(object$warmup))
    stop("this fit was created with save_warmup = FALSE, so full warmup draws were not ",
         "stored; refit with save_warmup = TRUE to inspect them. A thinned warmup trace is ",
         "in 'fit$warmup_trace' and the frozen tuning summaries in 'fit$adaptation'.",
         call. = FALSE)
  invisible(NULL)
}

as.array.stan4bartFit <- function (x, include_warmup = FALSE, ...)
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
  check_warmup_stored(x, include_warmup_orig)

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

as.matrix.stan4bartFit <- function (x, ...) 
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

# Liveness probe for a StoredBARTSampler external pointer. predictBART null-
# checks the pointer (init.cpp:335) BEFORE its early return on a NULL test
# matrix (init.cpp:340), so predictBART(ptr, NULL, NULL) is a cheap, side-
# effect-free liveness test: it returns NULL for a live pointer and errors for
# a dead one - the state a saveRDS/readRDS round trip leaves behind (a non-NULL
# externalptr whose address is NULL).
bart_pointer_is_live <- function(ptr) {
  if (is.null(ptr)) return(FALSE)
  !inherits(tryCatch(.Call(C_stan4bart_predictBART, ptr, NULL, NULL),
                     error = function(e) e), "error")
}

# Return a live StoredBARTSampler pointer for `object`, rebuilding it from the
# retained serializable state (object$state.bart) when the live pointer has
# died in a fresh session after reload. The training design matrix, dropped
# from the retained bundle to keep it n-independent, is re-spliced from
# $bartData. The rebuilt pointer is cached in object$bart_env (a reference
# cell), so a reloaded fit rebuilds at most once and holds the pointer for the
# session; the in-session path returns the original live pointer untouched.
getBartSampler <- function(object) {
  if (bart_pointer_is_live(object$sampler.bart))
    return(object$sampler.bart)
  env <- object$bart_env
  if (!is.null(env) && bart_pointer_is_live(env$ptr))
    return(env$ptr)
  if (is.null(object$state.bart))
    stop("bart component requires stan4bart to be called with `bart_args = list('keepTrees' = TRUE)`")
  restore <- object$state.bart
  data.bart <- restore$data
  data.bart@x <- object$bartData@x
  data.bart@x.test <- object$bartData@x.test
  ptr <- .Call(C_stan4bart_createStoredBARTSampler,
               restore$control, data.bart, restore$model, restore$state)
  if (!is.null(env)) env$ptr <- ptr
  ptr
}

# The sampling run's draw x chain counts, read off whichever per-draw BART
# component the fit retained. store = "fits" keeps bart_train and it is the
# authoritative source; store = "trees" drops it but always keeps varcount (and
# the parametric "stan" block), which carry the same two counts. Used wherever
# the readers previously indexed dim(object$bart_train)[2:3].
bart_sample_chain_dims <- function(object) {
  for (comp in list(object$bart_train, object$bart_varcount, object$stan)) {
    d <- dim(comp)
    if (!is.null(d)) return(c(d[2L], d[3L]))
  }
  stop("cannot determine the BART sample and chain counts: no per-draw BART ",
       "component was stored", call. = FALSE)
}

# TRUE when the stored BART block for `sample` is absent but the kept trees can
# reproduce it (store = "trees"): the reader must route through recompute. The
# design-matrix check keeps a fit with no test data (bart_test legitimately
# NULL) from being treated as recompute-able for sample = "test".
bart_needs_recompute <- function(object, sample) {
  stored <- if (sample == "train") object$bart_train else object$bart_test
  if (!is.null(stored) || is.null(object$state.bart)) return(FALSE)
  X <- if (sample == "train") object$bartData@x else object$bartData@x.test
  !is.null(X) && nrow(X) > 0L
}

# Reproduce the n x draws x chains bart_train / bart_test block a store = "trees"
# fit did not retain, by replaying the kept trees over the training / test
# design matrix through dbarts's predict path (predictBART, init.cpp:332).
# Predictions arrive on the original response scale (the restored state carries
# each chain's fit transform), shaped and named like the stored block. `rows`
# restricts to a row subset (row-block streaming for fitted()); NULL is the full
# matrix - the cost extract() documents.
recompute_bart_block <- function(object, sample, rows = NULL) {
  X <- if (sample == "train") object$bartData@x else object$bartData@x.test
  if (!is.null(rows)) X <- X[rows, , drop = FALSE]
  result <- .Call(C_stan4bart_predictBART, getBartSampler(object), X, NULL)
  if (length(dim(result)) == 2L) dim(result) <- c(dim(result), 1L)
  dimnames(result) <- list(observation = NULL, iterations = NULL,
                           chain = paste0("chain:", seq_len(dim(result)[3L])))
  result
}

extract.stan4bartFit <-
  function(object,
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
                    "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan",
                    "trees", "callback"),
           sample = c("train", "test"),
           combine_chains = TRUE,
           sample_new_levels = TRUE,
           include_warmup = FALSE,
           ...)
{
  matchedCall <- match.call()
  type <- match.arg(type)
  
  if (type == "trees") {
    if (is.null(object$sampler.bart))
      stop("extracting trees requires stan4bart to be called with `bart.args = list('keepTrees' == TRUE)`")
    dotsList <- list(...)
    treeNums <- if ("treeNums" %in% names(dotsList)) as.integer(dotsList[["treeNums"]]) else NULL
    chainNums <- if ("chainNums" %in% names(dotsList)) as.integer(dotsList[["chainNums"]]) else NULL
    sampleNums <- if ("sampleNums" %in% names(dotsList)) as.integer(dotsList[["sampleNums"]]) else NULL
    return(.Call(C_stan4bart_getTrees, getBartSampler(object), chainNums, sampleNums, treeNums, FALSE))
  } else {
    if (length(list(...)) > 0) warning("unused arguments ignored")
  }

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
  check_warmup_stored(object, include_warmup_orig)

  if (type == "callback") {
    if ((include_warmup && is.null(object$warmup$callback)) ||
        (!only_warmup && is.null(object$sample$callback)))
      stop("cannot extract callback samples for model fit without callback function")
    if (any(sapply(matchedCall[c("sample", "sample_new_levels")], function(x) !is.null(x)))) {
      warning("'sample' and 'sample_new_levels' arguments ignored when extracting callback samples")
    }
    result <- get_samples(object$callback, include_warmup, only_warmup)
    return(if (combine_chains) combine_chains_f(result) else result)
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
  } else if (type %in% c("varcount", "k", "stan")) {
    if (type %in% "varcount")
      result <- get_samples(object$bart_varcount, include_warmup, only_warmup)
    else
      result <- get_samples(object[[type]], include_warmup, only_warmup)
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
  
  bart_dims <- bart_sample_chain_dims(object)
  n_samples <- bart_dims[1L]
  n_chains  <- bart_dims[2L]
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
      n_obs_bart <- if (!is.null(object$bart_train)) dim(object$bart_train)[1L] else nrow(object$bartData@x)
      if (n_obs_bart != nrow(frame) && type == "indiv.bart")
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
    if (!is.null(offset) && offset_type %in% "bart" && type != "indiv.bart") {
      indiv.bart <- offset
    } else if (bart_needs_recompute(object, sample)) {
      # store = "trees": replay the kept trees over the (training / test) design
      # to reproduce the block that was not stored. The trees carry only the
      # post-warmup draws, so warmup access is not recompute-able.
      if (isTRUE(include_warmup) || only_warmup)
        stop("cannot return warmup BART draws for a fit created with store = \"trees\"; ",
             "the kept trees reproduce only the post-warmup draws", call. = FALSE)
      indiv.bart <- recompute_bart_block(object, sample)
    } else {
      indiv.bart <- if (sample == "train")
        get_samples(object$bart_train, include_warmup, only_warmup)
      else
        get_samples(object$bart_test, include_warmup, only_warmup)
    }

    if (inherits(na.action.bart, "exclude")) {
      indiv.bart.all <- array(NA_real_, c(n_all, dim(indiv.bart)[-1L]), dimnames(indiv.bart))
      indiv.bart.all[-na.action.bart,,] <- indiv.bart
      indiv.bart <- indiv.bart.all
    }
  }
  
  if (type %in% c("ev", "ppd")) {
    # indiv.bart/indiv.fixef/indiv.ranef (and the additive offset, when
    # offset_type == "default") are full n x draws x chains arrays. Once dim/
    # dimnames are attached, R's arithmetic can no longer reuse an operand's
    # buffer for the answer, so a chain of k "+"s costs k full-size
    # allocations. Stripping dim/dimnames before summing - same terms, same
    # left-to-right order - lets the unbound intermediate sums be reused in
    # place; dim/dimnames are reattached once, on the finished sum, instead of
    # on every partial one.
    result_dim <- dim(indiv.bart)
    result_dimnames <- dimnames(indiv.bart)
    dim(indiv.bart) <- dim(indiv.fixef) <- dim(indiv.ranef) <- NULL

    result <- if (!is.null(offset) && offset_type == "default")
      indiv.bart + indiv.fixef + indiv.ranef + offset
    else
      indiv.bart + indiv.fixef + indiv.ranef

    dim(result) <- result_dim
    dimnames(result) <- result_dimnames
  } else {
    result <- switch(type,
                     indiv.fixef = indiv.fixef,
                     indiv.ranef = indiv.ranef,
                     indiv.bart  = indiv.bart,
                     varcount    = object$bart_varcount,
                     ranef       = object$ranef,
                     fixef       = object$fixef,
                     Sigma       = object$Sigma,
                     sigma       = object$sigma,
                     k           = object$k)
  }

  if (type %in% c("ev", "ppd") && is_bernoulli)
    result <- pnorm(result)
  if (type %in% "ppd") {
    weights <- if (sample == "test") object$test$frame[["(weights)"]] else object$frame[["(weights)"]]
    if (is.null(weights) || length(weights) == 0L) {
      if (is_bernoulli) {
        # array(x, dim, dimnames) always copies x even when its length already
        # matches; dim<-/dimnames<- attach the same shape in place.
        rb <- rbinom(length(result), 1L, result)
        dim(rb) <- dim(result)
        dimnames(rb) <- dimnames(result)
        result <- rb
      } else {
        sigma <- extract(object, "sigma", combine_chains = TRUE)
        noise <- rnorm(dim(result)[1L] * n_samples * n_chains, 0,
                       rep(as.vector(sigma), each = dim(result)[1L]))
        dim(noise) <- c(dim(result)[1L], n_samples, n_chains)
        result <- result + noise
      }
    } else {
      if (is_bernoulli) {
        rb <- rbinom(length(result), 1L, result)
        dim(rb) <- dim(result)
        dimnames(rb) <- dimnames(result)
        result <- rb
        result <- weights * result
      } else {
        n_obs <- dim(result)[1L]
        n_total_samples <- n_samples * n_chains
        sigma <- extract(object, "sigma", combine_chains = TRUE)
        sigma <- rep_len(sigma, n_obs * n_total_samples) * rep(sqrt(1 / weights), each = n_total_samples)
        # The rnorm draws come out in (sample:chain, obs) order to match sigma
        # above - that order is fixed by the RNG stream position, so it is the
        # (same-size) noise array, not `result`, that gets permuted into
        # observation-major order. One aperm of `noise` replaces the original
        # two aperms of `result` at the same values, measurably faster even
        # though total bytes allocated is about the same.
        noise <- rnorm(n_obs * n_total_samples, 0, sigma)
        dim(noise) <- c(n_samples, n_chains, n_obs)
        result <- result + aperm(noise, c(3L, 1L, 2L))
      }
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



.message_env <- new.env(parent = emptyenv())

fitted.stan4bartFit <-
  function(object,
           type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
                    "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan",
                    "callback"),
           sample = c("train", "test"),
           sample_new_levels = TRUE,
           ...)
{
  if (length(list(...)) > 0) warning("unused arguments ignored")

  type <- match.arg(type)
  sample <- match.arg(sample)

  # skipped for weighted binomial, where the ppd mean is the weights times the
  # ev mean rather than the ev mean itself
  if (type == "ppd" && !isTRUE(.message_env$ppd_mean)) {
    weights <- if (sample == "test") object$test$frame[["(weights)"]] else object$frame[["(weights)"]]
    if (!(object$family$family == "binomial" && !is.null(weights) && length(weights) > 0L)) {
      .message_env$ppd_mean <- TRUE
      message("the \"ppd\" mean equals the \"ev\" mean; fitted(type = \"ev\") computes it exactly and faster (once per session)")
    }
  }

  # store = "trees": for the BART-bearing linear surfaces, stream the posterior
  # mean in row-blocks so peak memory stays block_n x draws x chains rather than
  # materializing the full recomputed block that extract() would (M2). ppd is a
  # sampled quantity and the parametric-only types never read the bart block, so
  # both fall through to the extract()-then-average path below.
  if (type %in% c("ev", "indiv.bart") && bart_needs_recompute(object, sample))
    return(fitted_bart_streamed(object, type, sample, sample_new_levels))

  samples <- extract(object, type, sample, combine_chains = TRUE, sample_new_levels = sample_new_levels)

  average_samples_f <- function(x) {
    if (is.array(x)) {
      d <- dim(x)
      l_d <- length(d)
      if (l_d == 2L) {
        # matches apply(x, 1L, mean) - row means over the combined
        # iterations:chains column margin
        rowMeans(x)
      } else {
        # matches apply(x, seq(1L, l_d - 1L), mean) - average over the last
        # (iterations:chains) margin only; reshape to a 2D matrix so a single
        # rowMeans pass replaces the per-margin apply loop, then restore shape
        dn <- dimnames(x)
        dim(x) <- c(prod(d[-l_d]), d[l_d])
        result <- rowMeans(x)
        dim(result) <- d[-l_d]
        if (!is.null(dn)) dimnames(result) <- dn[-l_d]
        result
      }
    } else if (is.matrix(x)) {
      rowMeans(x)
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

# Posterior-mean fitted surface for a store = "trees" fit, computed by streaming
# row-blocks so no full n x draws x chains block is ever held (M2). For each
# block the BART draws are recomputed from the kept trees and the parametric
# fixef/ranef surfaces are evaluated on the same rows (by subsetting their design
# matrices, which keeps them block-sized too); the block is transformed by the
# response family and reduced to one posterior mean per observation. Matches the
# store = "fits" fitted() to numerical tolerance. The rarer training-side cases
# extract() handles specially - a component-replacing offset or na.exclude gaps -
# fall back to the materializing path (correct, just not memory-bounded).
fitted_bart_streamed <- function(object, type, sample, sample_new_levels) {
  is_bernoulli <- object$family$family == "binomial"
  offset_type  <- object$offset_type

  na_bart <- if (sample == "train") attr(object$frame, "na.action.bart") else object$test$na.action.bart
  offset  <- if (sample == "train") object$offset else object$test$offset

  has_fixef_par <- any(grepl("^beta|gamma", dimnames(object$stan)[[1L]]))
  has_ranef     <- !is.null(object$reTrms) && length(object$reTrms) > 0L
  X.fixef <- if (sample == "train") object$X else object$test$X
  reTrms  <- if (sample == "train") object$reTrms else object$test$reTrms
  # streaming needs the fixef/ranef designs for the requested sample; if a
  # parametric contribution ev requires is present but its design is not
  # available to evaluate per-block (or an offset replaces a component, or
  # na.exclude leaves gaps), fall back to the full-materialization path so no
  # contribution is silently dropped.
  cannot_stream <- (type == "ev" && ((has_fixef_par && is.null(X.fixef)) ||
                                     (has_ranef && is.null(reTrms))))
  if (cannot_stream ||
      (!is.null(offset) && length(offset) > 0L && offset_type != "default") ||
      inherits(na_bart, "exclude")) {
    samples <- extract(object, type, sample, combine_chains = TRUE, sample_new_levels = sample_new_levels)
    return(if (!is.null(dim(samples))) apply(samples, seq_len(length(dim(samples)) - 1L), mean) else mean(samples))
  }

  X.bart <- if (sample == "train") object$bartData@x else object$bartData@x.test
  n_obs  <- nrow(X.bart)
  dims   <- bart_sample_chain_dims(object)
  M      <- dims[1L] * dims[2L]

  add_fixef <- type == "ev" && has_fixef_par && !is.null(X.fixef)
  add_ranef <- type == "ev" && has_ranef
  ev_offset <- if (type == "ev" && !is.null(offset) && length(offset) > 0L && offset_type == "default") offset else NULL

  # bound peak at ~2e6 doubles (~16 MB) of recomputed draws per block
  block   <- max(1L, min(n_obs, as.integer(ceiling(2e6 / max(1L, M)))))
  result  <- numeric(n_obs)
  start   <- 1L
  while (start <= n_obs) {
    idx <- seq.int(start, min(start + block - 1L, n_obs))
    ev_b <- recompute_bart_block(object, sample, rows = idx)
    if (add_fixef)
      ev_b <- ev_b + fitted_fixed(object, X.fixef[idx, , drop = FALSE], FALSE)
    if (add_ranef) {
      reTrms_b <- reTrms
      reTrms_b$Zt <- reTrms$Zt[, idx, drop = FALSE]
      ev_b <- ev_b + fitted_random(object, reTrms_b, FALSE, sample_new_levels)
    }
    if (!is.null(ev_offset)) ev_b <- ev_b + ev_offset[idx]
    if (type == "ev" && is_bernoulli) ev_b <- pnorm(ev_b)
    result[idx] <- rowMeans(matrix(ev_b, length(idx), M))
    start <- start + block
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
  n_chains  <- bart_sample_chain_dims(object)[2L]
  n_samples <- ncol(fixef) %/% n_chains
  
  array(t(t(x %*% fixef) - intercept_delta),
        c(n_obs, n_samples, n_chains),
          dimnames = list(observation = NULL, sample = NULL, chain = NULL))
}

fitted_random <- function(object, reTrms, include_warmup, sample_new_levels)
{
  # quick version if the new data have the same form as that fitted in the model
  if (identical(reTrms$cnms, object$reTrms$cnms) && all(row.names(reTrms$Zt) == row.names(object$reTrms$Zt))) {
    b_rows <- grep("^b\\.", dimnames(object$stan)[[1L]])
    b_mat <- matrix(object$stan[b_rows,,], length(b_rows), prod(dim(object$stan)[-1L]))
    return(array(
      Matrix::crossprod(reTrms$Zt, b_mat), c(ncol(reTrms$Zt), dim(object$stan)[2L:3L]),
        dimnames = list(observation = NULL, iterations = NULL, chain = dimnames(object$stan)$chain)
    ))
  }

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

predict.stan4bartFit <-
  # 'offset' must stay behind the defaulted arguments: undefaulted and third,
  # it captures a positionally-supplied 'type', e.g. predict(fit, newdata, "ev").
  function(object, newdata,
           type = c("ev", "ppd", "indiv.fixef", "indiv.ranef", "indiv.bart"),
           combine_chains = TRUE,
           sample_new_levels = TRUE,
           offset,
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
  
  bart_dims <- bart_sample_chain_dims(object)
  n_samples <- bart_dims[1L]
  n_chains  <- bart_dims[2L]
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
    if (!is.null(testData$X) && n_fixef > 0L) {
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
    # predictions arrive on the original response scale: the restored
    # sampler's state carries each chain's fit transform
    indiv.bart <- .Call(C_stan4bart_predictBART, getBartSampler(object), testData$X.bart, NULL)
    if (length(dim(indiv.bart)) == 2L)
      dim(indiv.bart) <- c(dim(indiv.bart), 1L)
    dimnames(indiv.bart) <-  list(observation = NULL, sample = NULL, chain = NULL)
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
  
  dimnames(result)[[3L]] <- paste0("chain:", seq_len(dim(result)[3L]))

  if (combine_chains) {
    if (is.list(result)) {
      result <- lapply(result, combine_chains_f)
  } else {
      result <- combine_chains_f(result)
    }
  }

  result
}

# Sampling-phase mean-leapfrog-steps-per-transition value above which the
# frozen WALNUTS tuning is treated as poor. Derived from a warmup-grid
# measurement (six reference tiers, warmup in {20,50,100,200,400,1000}, two
# builds): every tier's sampling phase came in at <= 24 once warmup reached
# the documented floor of 100, while cells below that floor ranged from
# comparably tuned (tiers that never needed much warmup) up past 50. 30 sits
# with margin above the well-tuned ceiling and below the clearly mis-tuned
# band, so it does not fire on any reference tier run at or past the floor.
mean_leapfrog_warn_threshold <- 30

# Emit the one-line mis-tuning note if any chain's sampling-phase mean
# leapfrog count exceeds the threshold above. A no-op when the fit has no
# adaptation record (warmup = 0, or an object predating the diagnostic).
warn_if_leapfrog_high <- function(object) {
  mean_leapfrog <- object$adaptation$mean_leapfrog
  if (is.null(mean_leapfrog) || !any(is.finite(mean_leapfrog)))
    return(invisible(NULL))

  worst <- max(mean_leapfrog, na.rm = TRUE)
  if (worst > mean_leapfrog_warn_threshold) {
    message(sprintf(
      "parametric sampler tuning looks poor (mean leapfrog steps per transition = %.1f); consider increasing warmup.",
      worst))
  }
  invisible(NULL)
}

print.stan4bartFit <- function(x, ...)
{
  cat("stan4bart model fit\n")
  cat("formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("family:  ", x$family$family, "\n", sep = "")
  if (!is.null(x$stan)) {
    dims <- dim(x$stan)
    cat(dims[3L], " chains, ", dims[2L], " samples per chain\n", sep = "")
  }

  warn_if_leapfrog_high(x)

  invisible(x)
}

summary.stan4bartFit <- function(object, ...)
{
  warn_if_leapfrog_high(object)

  result <- list(call = object$call, formula = object$formula,
                 family = object$family$family,
                 dim = if (!is.null(object$stan)) dim(object$stan)[2L:3L] else NULL,
                 adaptation = object$adaptation[intersect(
                   c("step_size", "mean_leapfrog", "mean_leapfrog_warmup"),
                   names(object$adaptation))])
  class(result) <- "summary.stan4bartFit"
  result
}

print.summary.stan4bartFit <- function(x, ...)
{
  cat("stan4bart model fit\n")
  cat("formula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("family:  ", x$family, "\n", sep = "")
  if (!is.null(x$dim))
    cat(x$dim[2L], " chains, ", x$dim[1L], " samples per chain\n", sep = "")
  if (!is.null(x$adaptation$mean_leapfrog)) {
    cat("\nmean leapfrog steps per transition, by chain:\n")
    cat("  warmup:   ", paste(round(x$adaptation$mean_leapfrog_warmup, 1), collapse = ", "), "\n", sep = "")
    cat("  sampling: ", paste(round(x$adaptation$mean_leapfrog, 1), collapse = ", "), "\n", sep = "")
  }

  invisible(x)
}

