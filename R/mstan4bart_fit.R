getSigma <- function(cnms, samples) {
  thetas <- samples[grep("^theta_L", rownames(samples)),,drop = FALSE]
  nc <- sapply(cnms, FUN = length)
  nms <- names(cnms)
  Sigma <- apply(thetas, 2L, function(theta) mkVarCorr(sc = 1, cnms, nc, theta, nms))
  Sigmas <- lapply(seq_along(Sigma[[1L]]), function(j) {
    raw <- sapply(Sigma, function(Sigma.i) Sigma.i[[j]])
    array(raw, c(NROW(Sigma[[1L]][[j]]), NCOL(Sigma[[1L]][[j]]), length(Sigma)),
          dimnames = list(rownames(Sigma[[1L]][[j]]), colnames(Sigma[[1L]][[j]]), NULL))
  })
  names(Sigmas) <- names(Sigma[[1L]])
  Sigmas
}
getFixef <- function(samples)
  samples[grep("^(?:beta|gamma)\\.", rownames(samples), perl = TRUE),,drop = FALSE]

getRanef <- function(group, samples) {
  ranef <- samples[grep("^b\\.", rownames(samples)),,drop = FALSE]
  numGroupingFactors <- length(group$cnms)
  numRanefPerGroupingFactor <- unname(sapply(group$cnms, length))
  
  res <- lapply(seq_len(numGroupingFactors), function(j) {
    ranef.group <- ranef[seq.int(group$Gp[j] + 1L, group$Gp[j + 1L]),,drop = FALSE]
    array(ranef.group, c(numRanefPerGroupingFactor[j], nrow(ranef.group) %/% numRanefPerGroupingFactor[j], ncol(ranef.group)),
          dimnames = list(group$cnms[[j]], levels(group$flist[[j]]), NULL))
  })
  names(res) <- names(group$cnms)
  res
}

# putting this out here so we can export it when parallelizing
mstan4bart_fitforreal <- function(chain.num, control.bart, data.bart, model.bart, data.stan, control.stan, control.common, group)
{
  chain.num <- "ignored"
  # TODO: figure out why the C refs aren't reachable when dispatched to cluster
  ns <- asNamespace("stan4bart")
  
  sampler <- .Call(ns$C_stan4bart_create, control.bart, data.bart, model.bart, data.stan, control.stan, control.common)
  if (control.common$verbose > 0L)
    .Call(ns$C_stan4bart_printInitialSummary, sampler)
  results <- list()
  if (control.common$warmup > 0L)
    results$warmup  <- .Call(ns$C_stan4bart_run, sampler, control.common$warmup, TRUE, "both")
  .Call(ns$C_stan4bart_disengageAdaptation, sampler)
  results$sample <- .Call(ns$C_stan4bart_run, sampler, control.common$iter - control.common$warmup,
                          FALSE, "both")
  
  
  if (control.common$warmup > 0L) {
    stan_warmup <- list(raw = results$warmup$stan)
    stan_warmup$Sigma <- ns$getSigma(group$cnms, stan_warmup$raw)
    stan_warmup$fixef <- ns$getFixef(stan_warmup$raw)
    stan_warmup$ranef <- ns$getRanef(group, stan_warmup$raw)
    results$warmup$stan <- stan_warmup
  }
  stan_sample <- list(raw = results$sample$stan)
  stan_sample$Sigma <- ns$getSigma(group$cnms, stan_sample$raw)
  stan_sample$fixef <- ns$getFixef(stan_sample$raw)
  stan_sample$ranef <- ns$getRanef(group, stan_sample$raw)
  results$sample$stan <- stan_sample
  
  if (control.bart@keepTrees) {
    results$state.bart <- .Call(ns$C_stan4bart_exportBARTState, sampler)
    results$range.bart <- .Call(ns$C_stan4bart_getBARTDataRange, sampler)
  }
  
  results
}

mstan4bart_fit <- 
  function(object,
           family,
           bart_offset_init, sigma_init,
           verbose,
           iter,
           warmup,
           thin,
           chains,
           cores,
           refresh,
           stan_args,
           bart_args)
{  
  supported_families <- c("binomial", "gaussian")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  
  supported_links <- supported_glm_links(supported_families[fam])
  link <- which(supported_links == family$link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (is.null(stan_args))
    stan_args <- list()
  if (is.null(stan_args[["prior_covariance"]]))
    stan_args$prior_covariance <- decov()
  decov <- stan_args[["prior_covariance"]]
  if (is.null(stan_args[["prior"]]))
    stan_args$prior <- default_prior_coef_gaussian()
  if (is.null(stan_args[["prior_intercept"]]))
    stan_args$prior_intercept <- default_prior_intercept_gaussian()
  if (is.null(stan_args[["prior_aux"]]))
    stan_args$prior_aux <- exponential(autoscale = TRUE)
  
  
  x <- object$X
  y <- object$y
  data.bart   <- object$bartData
  weights     <- object$weights
  offset      <- object$offset
  offset_type <- object$offset_type
  group       <- object$reTrms
  
  
  x_stuff <- center_x(x, FALSE)
  xtemp <- has_intercept <- NULL # R CMD check
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars.fixef <- ncol(xtemp)
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso", "product_normal")
  ok_intercept_dists <- ok_dists[1:3]
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  
  prior_stuff <- handle_glm_prior(
    stan_args$prior,
    nvars = nvars.fixef,
    default_scale = 2.5,
    link = family$link, 
    ok_dists = ok_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}, 
  # global_prior_df, global_prior_scale, slab_df, slab_scale
  for (name in names(prior_stuff))
    stan_args[[name]] <- prior_stuff[[name]]
  
  if (isTRUE(is.list(stan_args$prior_intercept)) && 
      isTRUE(stan_args$prior_intercept$default)) {
    m_y <- 0
    if (family$family == "gaussian" && family$link == "identity") {
      if (!is.null(y)) m_y <- mean(y) # y can be NULL if prior_PD=TRUE
    }
    stan_args$prior_intercept$location <- m_y
  }
  prior_intercept_stuff <- handle_glm_prior(
    stan_args$prior_intercept,
    nvars = 1,
    default_scale = 2.5,
    link = family$link,
    ok_dists = ok_intercept_dists
  )
   # prior_{dist, mean, scale, df, dist_name, autoscale}_for_intercept
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (name in names(prior_intercept_stuff))
    stan_args[[name]] <- prior_intercept_stuff[[name]]
  
  prior_aux_stuff <-
    handle_glm_prior(
      stan_args$prior_aux,
      nvars = 1,
      default_scale = 1,
      link = NULL, # don't need to adjust scale based on logit vs probit
      ok_dists = ok_aux_dists
    )
  
  # prior_{dist, mean, scale, df, dist_name, autoscale}_for_aux
  names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
  if (is.null(stan_args$prior_aux))
    prior_aux_stuff$prior_scale_for_aux <- Inf
  
  for (name in names(prior_aux_stuff)) 
    stan_args[[name]] <-  prior_aux_stuff[[name]]
  
  famname <- supported_families[fam]
  is_bernoulli <- famname == "binomial"
  is_gaussian <- !is_bernoulli
  is_continuous <- is_gaussian
  
  # require intercept for certain family and link combinations
  if (!has_intercept) {
    linkname <- supported_links[link]
    needs_intercept <- !is_gaussian && linkname == "identity" ||
      # is_gamma && linkname == "inverse" ||
      is.binomial(famname) && linkname == "log"
    if (needs_intercept)
      stop("to use this combination of family and link ", 
           "the model must have an intercept")
  }
  
  if (is_gaussian) {
    ss <- sd(y)
    if (stan_args$prior_dist > 0L && stan_args$prior_autoscale) 
      stan_args$prior_scale <- ss * stan_args$prior_scale
    if (stan_args$prior_dist_for_intercept > 0L && stan_args$prior_autoscale_for_intercept) 
      stan_args$prior_scale_for_intercept <- ss * stan_args$prior_scale_for_intercept
    if (stan_args$prior_dist_for_aux > 0L && stan_args$prior_autoscale_for_aux)
      stan_args$prior_scale_for_aux <- ss * stan_args$prior_scale_for_aux
  }
  
  #if (!QR && prior_dist > 0L && prior_autoscale) {
  if (stan_args$prior_dist > 0L && stan_args$prior_autoscale) {
    min_prior_scale <- 1e-12
    stan_args$prior_scale <-
      pmax(min_prior_scale, stan_args$prior_scale / 
             apply(xtemp, 2L, FUN = function(x) {
               num.categories <- length(unique(x))
               x.scale <- 1
               if (num.categories == 1) {
                 x.scale <- 1
               } else {
                 x.scale <- sd(x)
               }
               return(x.scale)
            }))
  }
  
  stan_args$prior_scale <- 
    as.array(pmin(.Machine$double.xmax, stan_args$prior_scale))
  stan_args$prior_scale_for_intercept <- 
    min(.Machine$double.xmax, stan_args$prior_scale_for_intercept)
  
  #if (QR) {
  #  if (ncol(xtemp) <= 1)
  #    stop("'QR' can only be specified when there are multiple predictors.")
  #  if (sparse)
  #    stop("'QR' and 'sparse' cannot both be TRUE.")
  #  cn <- colnames(xtemp)
  #  decomposition <- qr(xtemp)
  #  Q <- qr.Q(decomposition)
  #  if (prior_autoscale) scale_factor <- sqrt(nrow(xtemp) - 1L)
  #  else scale_factor <- diag(qr.R(decomposition))[ncol(xtemp)]
  #  R_inv <- qr.solve(decomposition, Q) * scale_factor
  #  xtemp <- Q * scale_factor
  #  colnames(xtemp) <- cn
  #  xbar <- c(xbar %*% R_inv)
  #}
    
  if (length(weights) > 0L && all(weights == 1)) weights <- double()
  if (length(offset)  > 0L && all(offset  == 0)) offset  <- double()
  
  # create entries in the data block of the .stan file
  data.stan <- with(stan_args, nlist(
    N = length(y),
    has_weights = length(weights) > 0L,
    prior_dist_for_aux = prior_dist_for_aux,
    N = nrow(xtemp),
    K = ncol(xtemp),
    #family = stan_family_number(famname),  # hard-coded in model
    #link,                                  # hard-coded in model
    has_intercept,
    prior_dist,
    prior_mean,
    prior_scale,
    prior_df,
    prior_dist_for_intercept,
    prior_scale_for_intercept = c(prior_scale_for_intercept),
    prior_mean_for_intercept = c(prior_mean_for_intercept),
    prior_df_for_intercept = c(prior_df_for_intercept), 
    global_prior_df, global_prior_scale, slab_df, slab_scale, # for hs priors
    prior_df_for_intercept = c(prior_df_for_intercept),
    prior_dist_for_aux = prior_dist_for_aux,
    num_normals = if(prior_dist == 7) as.integer(prior_df) else integer(0)
    # mean,df,scale for aux added below depending on family
  ))

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- decov
  
  if (length(group) > 0 && length(group$flist) > 0) {
    check_reTrms(group)
    
    Z <- t(group$Zt)
    
    # don't need padding for now
    #group <-
    #  pad_reTrms(Ztlist = group$Ztlist,
    #             cnms = group$cnms,
    #             flist = group$flist)
    # Z <- group$Z
  
    p <- sapply(group$cnms, FUN = length)
    l <- sapply(attr(group$flist, "assign"), function(i) 
      nlevels(group$flist[[i]]))
    t <- length(l)
    b_nms <- make_b_nms(group)
    g_nms <- unlist(lapply(1:t, FUN = function(i) {
      paste(group$cnms[[i]], names(group$cnms)[i], sep = "|")
    }))
    data.stan$t <- t
    data.stan$p <- as.array(p)
    data.stan$l <- as.array(l)
    data.stan$q <- ncol(Z)
    data.stan$len_theta_L <- sum(choose(p, 2), p)
    if (is_bernoulli) {
      stop("implement me!")
      parts0 <- extract_sparse_parts(Z[y == 0, , drop = FALSE])
      parts1 <- extract_sparse_parts(Z[y == 1, , drop = FALSE])
      data.stan$num_non_zero <- c(length(parts0$w), length(parts1$w))
      data.stan$w0 <- as.array(parts0$w)
      data.stan$w1 <- as.array(parts1$w)
      data.stan$v0 <- as.array(parts0$v - 1L)
      data.stan$v1 <- as.array(parts1$v - 1L)
      data.stan$u0 <- as.array(parts0$u - 1L)
      data.stan$u1 <- as.array(parts1$u - 1L)
    } else {
      parts <- extract_sparse_parts(Z)
      data.stan$num_non_zero <- length(parts$w)
      data.stan$w <- parts$w
      data.stan$v <- parts$v - 1L
      data.stan$u <- parts$u - 1L
    }
    data.stan$shape <- as.array(maybe_broadcast(decov$shape, t))
    data.stan$scale <- as.array(maybe_broadcast(decov$scale, t))
    data.stan$len_concentration <- sum(p[p > 1])
    data.stan$concentration <- 
      as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
    data.stan$len_regularization <- sum(p > 1)
    data.stan$regularization <- 
      as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
  } else {
    data.stan$t <- 0L
    data.stan$p <- integer(0)
    data.stan$l <- integer(0)
    data.stan$q <- 0L
    data.stan$len_theta_L <- 0L
    if (is_bernoulli) {
      stop("implement me!")
      data.stan$num_non_zero <- rep(0L, 2)
      data.stan$w0 <- data.stan$w1 <- double(0)
      data.stan$v0 <- data.stan$v1 <- integer(0)
      data.stan$u0 <- data.stan$u1 <- integer(0)
    } else {
      data.stan$num_non_zero <- 0L
      data.stan$w <- double(0)
      data.stan$v <- integer(0)
      data.stan$u <- integer(0)
    }
    data.stan$shape <- data.stan$scale <- data.stan$concentration <-
      data.stan$regularization <- rep(0, 0)
    data.stan$len_concentration <- 0L
    data.stan$len_regularization <- 0L
  }
  
  if (!is_bernoulli) {
    data.stan$X <- xtemp
    
    data.stan$y <- y
    data.stan$weights <- weights
    data.stan$offset_ <- numeric(length(y))
  }

  if (is_continuous) {
    data.stan$ub_y <- Inf
    data.stan$lb_y <- if (is_gaussian) -Inf else 0
    data.stan$prior_scale_for_aux <- stan_args$prior_scale_for_aux %ORifINF% 0
    data.stan$prior_df_for_aux <- c(stan_args$prior_df_for_aux)
    data.stan$prior_mean_for_aux <- c(stan_args$prior_mean_for_aux)
    data.stan$len_y <- length(y)
  } else if (is_bernoulli) {
    data.stan$prior_scale_for_aux <- 
      if (!length(group) || stan_args$prior_scale_for_aux == Inf) 
        0 else stan_args$prior_scale_for_aux
    data.stan$prior_mean_for_aux <- 0
    data.stan$prior_df_for_aux <- 0
    if (is_bernoulli) {
      stop("implement me!")
      y0 <- y == 0
      y1 <- y == 1
      data.stan$N <- c(sum(y0), sum(y1))
      data.stan$X0 <- array(xtemp[y0, , drop = FALSE], dim = c(1, sum(y0), ncol(xtemp)))
      data.stan$X1 <- array(xtemp[y1, , drop = FALSE], dim = c(1, sum(y1), ncol(xtemp)))
      if (length(weights) > 0L) {
        # nocov start
        # this code is unused because weights are interpreted as number of 
        # trials for binomial glms
        data.stan$weights0 <- weights[y0]
        data.stan$weights1 <- weights[y1]
        # nocov end
      } else {
        data.stan$weights0 <- double(0)
        data.stan$weights1 <- double(0)
      }
      if (length(offset) > 0L) {
        # TODO: fix this!
        # nocov start
        data.stan$offset0 <- offset[y0]
        data.stan$offset1 <- offset[y1]
        # nocov end
      } else {
        data.stan$offset0 <- double(0)
        data.stan$offset1 <- double(0)
      }
    }
  }
    
  prior_info <- with(stan_args, summarize_glm_prior(
    user_prior = prior_stuff,
    user_prior_intercept = prior_intercept_stuff,
    user_prior_aux = prior_aux_stuff,
    user_prior_covariance = prior_covariance,
    has_intercept = has_intercept,
    has_predictors = nvars.fixef > 0,
    adjusted_prior_scale = prior_scale,
    adjusted_prior_intercept_scale = prior_scale_for_intercept,
    adjusted_prior_aux_scale = prior_scale_for_aux,
    family = family
  ))
 
  for (varName in names(data.stan))
    if (is.logical(data.stan[[varName]]))
      data.stan[[varName]] <- as.integer(data.stan[[varName]])
  for (varName in c("len_theta_L"))
    data.stan[[varName]] <- as.integer(data.stan[[varName]])
  
   
  if (!is.numeric(iter) || length(iter) != 1L || iter <= 0)
    stop("'iter' must be a positive integer")
  iter <- as.integer(iter) # could include a warning about info lost if coercing from double
  if (!is.numeric(warmup) || length(warmup) != 1L || warmup < 0)
    stop("'warmup' must be a non-negative integer")
  warmup <- as.integer(warmup)
  if (!is.numeric(thin)   || length(thin) < 1L || length(thin) > 2L || any(thin <= 0))
    stop("'thin' must be one or two positive integers")
  if (!is.null(names(thin)) && c("bart", "stan") %in% names(thin)) {
    thin.bart <- thin[["bart"]]
    thin.stan <- thin[["stan"]]
  } else {
    thin.bart <- thin[1L]
    thin.stan <- if (length(thin) > 1L) thin[2L] else thin[1L]
  }
  thin.bart <- as.integer(thin.bart)
  thin.stan <- as.integer(thin.stan)
  
  if (!is.numeric(chains) || length(chains) != 1L || chains <= 0)
    stop("'chains' must be a positive integer")
  chains <- as.integer(chains)
  if (!is.numeric(cores) || length(cores) != 1L || cores <= 0)
    stop("'cores' must be a positive integer")
  cores <- as.integer(cores)
  
  
  if (is.logical(verbose)) verbose <- as.integer(verbose)
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("'verbose' must an integer")
  verbose <- as.integer(verbose)
  
  if (is.na(refresh)) refresh <- max(iter %/% 10L, 1L)
  if (!is.numeric(refresh) || length(refresh) != 1L || (!is.na(refresh) && refresh < 0))
    stop("'refresh' must be a non-negative integer")
  refresh <- as.integer(refresh)
  
  offset_type <- which(offset_type == c("default", "fixef", "ranef", "bart")) - 1L
  
  if (is.null(bart_args)) bart_args <- list()
  data.bart@sigma <- sigma_init
  control_call <- quote(dbarts::dbartsControl(n.chains = 1L, n.samples = 1L, n.burn = 0L,
                                              n.thin = thin.bart, n.threads = 1L,
                                              updateState = FALSE))
  for (name in intersect(names(bart_args), setdiff(names(formals(dbarts::dbartsControl)),
                                                   names(control_call))))
  {
    control_call[[name]] <- bart_args[[name]]
  }
  control.bart <- eval(control_call)
  if (!is.null(bart_args[["n.cuts"]]))
    attr(control.bart, "n.cuts") <- bart_args[["n.cuts"]]
  
  data.bart@n.cuts <- rep_len(attr(control.bart, "n.cuts"), ncol(data.bart@x))
  control.bart@binary <- !is_continuous
  evalEnv <- sys.frame(sys.nframe())
  
  cgm <- normal <- fixed <- NULL # R CMD check
  bartPriors <- dbarts:::parsePriors(control.bart, data.bart, cgm, normal, fixed(1), evalEnv)
  model.bart <- new("dbartsModel", bartPriors$tree.prior, bartPriors$node.prior, bartPriors$resid.prior,
                    node.scale = if (!is_continuous) 3.0 else 0.5)
  
  control.stan <- list(
    seed = sample.int(.Machine$integer.max, 1L),
    init_r = stan_args[["init_r"]] %ORifNULL% 2.0,
    thin = thin.stan,
    adapt_gamma = stan_args[["adapt_gamma"]],
    adapt_delta = stan_args[["adapt_delta"]],
    adapt_kappa = stan_args[["adapt_kappa"]]
  )
  
  control.common <- nlist(iter, warmup, verbose, refresh,
                          offset = offset, offset_type = offset_type,
                          bart_offset_init = bart_offset_init, sigma_init = sigma_init)
  
  
  chainResults <- vector("list", chains)
  runSingleThreaded <- cores <= 1L || chains <= 1L
  if (!runSingleThreaded) {
    tryResult <- tryCatch(cluster <- makeCluster(min(cores, chains), "PSOCK"), error = function(e) e)
    if (is(tryResult, "error"))
      tryResult <- tryCatch(cluster <- makeCluster(min(cores, chains), "FORK"), error = function(e) e)
    
    if (is(tryResult, "error")) {
      warning("unable to multithread, defaulting to single: ", tryResult$message)
      runSingleThreaded <- TRUE
    } else {
      if (control.common$verbose > 0L)
        cat("starting multithreaded fit, futher output silenced\n")
      control.common$verbose <- 0L
      
      clusterExport(cluster, "mstan4bart_fitforreal", asNamespace("stan4bart"))
      clusterEvalQ(cluster, require(stan4bart))
      
      tryResult <- tryCatch(
        chainResults <- clusterMap(cluster, "mstan4bart_fitforreal", seq_len(chains), MoreArgs = nlist(control.bart, data.bart, model.bart, data.stan, control.stan, control.common, group)),
                            error = function(e) e)
    
      stopCluster(cluster)
      
      if (is(tryResult, "error")) {
        warning("error running multithreaded, defaulting to single: ", tryResult$message)
        runSingleThreaded <- TRUE
      }
    }
  }
  
  if (runSingleThreaded) {
    for (chainNum in seq_len(chains))
      chainResults[[chainNum]] <- mstan4bart_fitforreal(1L, control.bart, data.bart, model.bart, data.stan, control.stan, control.common, group)
  }
  
  if (!is.null(chainResults[[1L]]$state.bart)) {
    all_state <- chainResults[[1L]]$state.bart
    if (chains > 1L) for (i in seq.int(2L, chains))
      all_state[[i]] <- chainResults[[i]]$state.bart[[1L]]
    
    
    attr(chainResults, "sampler.bart") <- 
      .Call(C_stan4bart_createStoredBARTSampler, control.bart, data.bart, model.bart, all_state)
  }
  
  chainResults
}



# internal ----------------------------------------------------------------

# @param famname string naming the family
# @return character vector of supported link functions for the family
supported_glm_links <- function(famname) {
  switch(
    famname,
    binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
    gaussian = c("identity", "log", "inverse"),
    Gamma = c("identity", "log", "inverse"),
    inverse.gaussian = c("identity", "log", "inverse", "1/mu^2"),
    "neg_binomial_2" = , # intentional
    poisson = c("log", "identity", "sqrt"),
    "Beta regression" = c("logit", "probit", "cloglog", "cauchit"),
    stop("unsupported family")
  )
}

# Family number to pass to Stan
# @param famname string naming the family
# @return an integer family code
stan_family_number <- function(famname) {
  switch(
    famname,
    "gaussian" = 1L,
    "Gamma" = 2L,
    "inverse.gaussian" = 3L,
    "beta" = 4L,
    "Beta regression" = 4L,
    "binomial" = 5L,
    "poisson" = 6L,
    "neg_binomial_2" = 7L,
    stop("Family not valid.")
  )
}



# Verify that outcome values match support implied by family object
#
# @param y outcome variable
# @param family family object
# @return y (possibly slightly modified) unless an error is thrown
#
validate_glm_outcome_support <- function(y, family) {
  if (is.null(y)) {
    return(y)
  }
  
  .is_count <- function(x) {
    all(x >= 0) && all(abs(x - round(x)) < .Machine$double.eps^0.5)
  }
  
  fam <- family$family
  
  if (!is.binomial(fam)) {
    # make sure y has ok dimensions (matrix only allowed for binomial models)
    if (length(dim(y)) > 1) {
      if (NCOL(y) == 1) {
        y <- y[, 1]
      } else {
        stop("Except for binomial models the outcome variable ",
             "should not have multiple columns.", 
             call. = FALSE)
      }
    }
    
    # check that values match support for non-binomial models
    if (is.gaussian(fam)) {
      return(y)
    } else if (is.gamma(fam) && any(y <= 0)) {
      stop("All outcome values must be positive for gamma models.", 
           call. = FALSE)
    } else if (is.ig(fam) && any(y <= 0)) {
      stop("All outcome values must be positive for inverse-Gaussian models.", 
           call. = FALSE)
    } else if (is.poisson(fam) && !.is_count(y)) {
      stop("All outcome values must be counts for Poisson models",
           call. = FALSE)
    } else if (is.nb(fam) && !.is_count(y)) {
      stop("All outcome values must be counts for negative binomial models",
           call. = FALSE)
    }
  } else { # binomial models
    if (NCOL(y) == 1L) {
      if (is.numeric(y) || is.logical(y)) 
        y <- as.integer(y)
      if (is.factor(y)) 
        y <- fac2bin(y)
      if (!all(y %in% c(0L, 1L))) 
        stop("All outcome values must be 0 or 1 for Bernoulli models.", 
             call. = FALSE)
    } else if (isTRUE(NCOL(y) == 2L)) {
      if (!.is_count(y))
        stop("All outcome values must be counts for binomial models.",
             call. = FALSE)
    } else {
      stop("For binomial models the outcome should be a vector or ",
           "a matrix with 2 columns.", 
           call. = FALSE)
    }
  }
  
  return(y)
}

# Generate fake y variable to use if prior_PD and no y is specified
# @param N number of observations
# @param family family object
fake_y_for_prior_PD <- function(N, family) {
  fam <- family$family
  if (is.gaussian(fam)) {
    # if prior autoscaling is on then the value of sd(y) matters
    # generate a fake y so that sd(y) is 1
    fake_y <- as.vector(scale(rnorm(N)))
  } else if (is.binomial(fam) || is.poisson(fam) || is.nb(fam)) {
    # valid for all discrete cases
    fake_y <- rep_len(c(0, 1), N)
  } else {
    # valid for gamma, inverse gaussian, beta 
    fake_y <- runif(N)
  }
  return(fake_y)
}



# Add extra level _NEW_ to each group
# 
# @param Ztlist ranef indicator matrices
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], Matrix::Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
  }
  Z <- t(do.call(rbind, args = Ztlist))
  return(nlist(Z, cnms, flist))
}

# Drop the extra reTrms from a matrix x
#
# @param x A matrix or array (e.g. the posterior sample or matrix of summary
#   stats)
# @param columns Do the columns (TRUE) or rows (FALSE) correspond to the 
#   variables?
make_b_nms <- function(group, m = NULL, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  m_stub <- if (!is.null(m)) get_m_stub(m, stub = stub) else NULL
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(m_stub, nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(paste0(m_stub, nms_i), paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}


summarize_glm_prior <-
  function(user_prior,
           user_prior_intercept,
           user_prior_aux,
           user_prior_covariance,
           has_intercept, 
           has_predictors,
           adjusted_prior_scale,
           adjusted_prior_intercept_scale, 
           adjusted_prior_aux_scale,
           family) {
    rescaled_coef <-
      user_prior$prior_autoscale && 
      has_predictors &&
      !is.na(user_prior$prior_dist_name) &&
      !all(user_prior$prior_scale == adjusted_prior_scale)
    rescaled_int <-
      user_prior_intercept$prior_autoscale_for_intercept &&
      has_intercept &&
      !is.na(user_prior_intercept$prior_dist_name_for_intercept) &&
      (user_prior_intercept$prior_scale_for_intercept != adjusted_prior_intercept_scale)
    rescaled_aux <- user_prior_aux$prior_autoscale_for_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_aux) &&
      (user_prior_aux$prior_scale_for_aux != adjusted_prior_aux_scale)
    
    if (has_predictors && user_prior$prior_dist_name %in% "t") {
      if (all(user_prior$prior_df == 1)) {
        user_prior$prior_dist_name <- "cauchy"
      } else {
        user_prior$prior_dist_name <- "student_t"
      }
    }
    if (has_intercept &&
        user_prior_intercept$prior_dist_name_for_intercept %in% "t") {
      if (all(user_prior_intercept$prior_df_for_intercept == 1)) {
        user_prior_intercept$prior_dist_name_for_intercept <- "cauchy"
      } else {
        user_prior_intercept$prior_dist_name_for_intercept <- "student_t"
      }
    }
    if (user_prior_aux$prior_dist_name_for_aux %in% "t") {
      if (all(user_prior_aux$prior_df_for_aux == 1)) {
        user_prior_aux$prior_dist_name_for_aux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_aux <- "student_t"
      }
    }
    prior_list <- list(
      prior = 
        if (!has_predictors) NULL else with(user_prior, list(
          dist = prior_dist_name,
          location = prior_mean,
          scale = prior_scale,
          adjusted_scale = if (rescaled_coef)
            adjusted_prior_scale else NULL,
          df = if (prior_dist_name %in% c
                   ("student_t", "hs", "hs_plus", "lasso", "product_normal"))
            prior_df else NULL
        )),
      prior_intercept = 
        if (!has_intercept) NULL else with(user_prior_intercept, list(
          dist = prior_dist_name_for_intercept,
          location = prior_mean_for_intercept,
          scale = prior_scale_for_intercept,
          adjusted_scale = if (rescaled_int)
            adjusted_prior_intercept_scale else NULL,
          df = if (prior_dist_name_for_intercept %in% "student_t")
            prior_df_for_intercept else NULL
        ))
    )
    if (length(user_prior_covariance))
      prior_list$prior_covariance <- user_prior_covariance
    
    aux_name <- .rename_aux(family)
    prior_list$prior_aux <- if (is.na(aux_name)) 
      NULL else with(user_prior_aux, list(
        dist = prior_dist_name_for_aux,
        location = if (!is.na(prior_dist_name_for_aux) && 
                       prior_dist_name_for_aux != "exponential")
          prior_mean_for_aux else NULL,
        scale = if (!is.na(prior_dist_name_for_aux) && 
                    prior_dist_name_for_aux != "exponential")
          prior_scale_for_aux else NULL,
        adjusted_scale = if (rescaled_aux)
          adjusted_prior_aux_scale else NULL,
        df = if (!is.na(prior_dist_name_for_aux) && 
                 prior_dist_name_for_aux %in% "student_t")
          prior_df_for_aux else NULL, 
        rate = if (!is.na(prior_dist_name_for_aux) && 
                   prior_dist_name_for_aux %in% "exponential")
          1 / prior_scale_for_aux else NULL,
        aux_name = aux_name
      ))
      
    return(prior_list)
  }

# rename aux parameter based on family
.rename_aux <- function(family) {
  fam <- family$family
  if (fam == "gaussian") "sigma" else
    if (fam == "Gamma") "shape" else
      if (fam == "inverse.gaussian") "lambda" else 
        if (fam == "neg_binomial_2") "reciprocal_dispersion" else NA
}

.sample_indices <- function(wts, n_draws) {
  ## Stratified resampling
  ##   Kitagawa, G., Monte Carlo Filter and Smoother for Non-Gaussian
  ##   Nonlinear State Space Models, Journal of Computational and
  ##   Graphical Statistics, 5(1):1-25, 1996.
  K <- length(wts)
  w <- n_draws * wts # expected number of draws from each model
  idx <- rep(NA, n_draws)

  c <- 0
  j <- 0

  for (k in 1:K) {
    c <- c + w[k]
    if (c >= 1) {
      a <- floor(c)
      c <- c - a
      idx[j + 1:a] <- k
      j <- j + a
    }
    if (j < n_draws && c >= runif(1)) {
      c <- c - 1
      j <- j + 1
      idx[j] <- k
    }
  }
  return(idx)
}

extract_sparse_parts <- function (A) 
{
  if (!requireNamespace("Matrix")) 
    stop("You have to install the Matrix package to call 'extract_sparse_parts'")
  if (!is(A, "Matrix")) 
    A <- Matrix::Matrix(A, sparse = TRUE)
  A <- Matrix::t(A)
  A <- as(A, "dgCMatrix")
  
  return(list(w = A@x, v = A@i + 1L, u = A@p + 1L))
  # return(.Call(rstan:::extract_sparse_components, A))
}
