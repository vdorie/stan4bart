getSigma <- function(cnms, samples) {
  thetas <- samples[grep("^theta_L", rownames(samples)),]
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
getRanef <- function(group, samples) {
  ranef <- samples[grep("^b\\.", rownames(samples)),]
  numGroupingFactors <- length(group$cnms)
  numRanefPerGroupingFactor <- unname(sapply(group$cnms, length))
  
  res <- lapply(seq_len(numGroupingFactors), function(j) {
    ranef.group <- ranef[seq.int(group$Gp[j] + 1L, group$Gp[j + 1L]),]
    array(ranef.group, c(numRanefPerGroupingFactor[j], nrow(ranef.group) %/% numRanefPerGroupingFactor[j], ncol(ranef.group)),
          dimnames = list(group$cnms[[j]], levels(group$flist[[j]]), NULL))
  })
  names(res) <- names(group$cnms)
  res
}

# putting this out here so we can export it when parallelizing
mstan4bart_fitforreal <- function(chain.num, bartControl, bartData, bartModel, standata, stan_args, commonControl, group)
{
  chain.num <- "ignored"
  # TODO: figure out why the C refs aren't reachable when dispatched to cluster
  ns <- asNamespace("stan4bart")
  sampler <- .Call(ns$C_stan4bart_create, bartControl, bartData, bartModel, standata, stan_args, commonControl)
  if (commonControl$verbose > 0L)
    .Call(ns$C_stan4bart_printInitialSummary, sampler)
  results <- list()
  if (commonControl$warmup > 0L)
    results$warmup  <- .Call(ns$C_stan4bart_run, sampler, commonControl$warmup, TRUE, "both")
  .Call(ns$C_stan4bart_disengageAdaptation, sampler)
  results$sample <- .Call(ns$C_stan4bart_run, sampler, commonControl$iter, FALSE, "both")
  
  if (commonControl$warmup > 0L) {
    stan_warmup <- list(raw = results$warmup$stan)
    stan_warmup$Sigma <- ns$getSigma(group$cnms, stan_warmup$raw)
    stan_warmup$ranef <- ns$getRanef(group, stan_warmup$raw)
    results$warmup$stan <- stan_warmup
  }
  stan_sample <- list(raw = results$sample$stan)
  stan_sample$Sigma <- ns$getSigma(group$cnms, stan_sample$raw)
  stan_sample$ranef <- ns$getRanef(group, stan_sample$raw)
  results$sample$stan <- stan_sample
  
  results
}

mstan4bart_fit <- 
  function(bartData, y,
           ranef_init, sigma_init,
           weights = rep(1, NROW(y)), 
           offset = rep(0, NROW(y)), 
           family = family,
           ...,
           prior_aux = exponential(),
           group = list(),
           adapt_delta = NULL)
{  
  supported_families <- c("binomial", "gaussian")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  
  # useless assignments to pass R CMD check
  prior_df_for_aux <- prior_dist_for_aux <- prior_mean_for_aux <- prior_scale_for_aux <-
    prior_autoscale_for_aux <- NULL
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus", 
                    "laplace", "lasso", "product_normal")
  ok_aux_dists <- c(ok_dists[1:3], exponential = "exponential")
  
  prior_aux_stuff <-
    handle_glm_prior(
      prior_aux,
      nvars = 1,
      default_scale = 1,
      link = NULL, # don't need to adjust scale based on logit vs probit
      ok_dists = ok_aux_dists
    )
  
  names(prior_aux_stuff) <- paste0(names(prior_aux_stuff), "_for_aux")
  if (is.null(prior_aux))
    prior_aux_stuff$prior_scale_for_aux <- Inf
  
  for (i in names(prior_aux_stuff)) 
    assign(i, prior_aux_stuff[[i]])
  
  famname <- supported_families[fam]
  is_bernoulli <- famname == "binomial"
  is_gaussian <- !is_bernoulli
  is_continuous <- is_gaussian
  
  if (is_gaussian) {
    ss <- sd(y)
    if (prior_dist_for_aux > 0L && prior_autoscale_for_aux)
      prior_scale_for_aux <- ss * prior_scale_for_aux
  }
    
  if (length(weights) > 0L && all(weights == 1)) weights <- double()
  if (length(offset)  > 0L && all(offset  == 0)) offset  <- double()
  
  # create entries in the data block of the .stan file
  standata <- nlist(
    N = length(y),
    has_weights = length(weights) > 0L,
    prior_dist_for_aux = prior_dist_for_aux
  )

  # make a copy of user specification before modifying 'group' (used for keeping
  # track of priors)
  user_covariance <- if (length(group) == 0L) NULL else group[["decov"]]
  
  check_reTrms(group)
  decov <- group$decov
  
  Z <- t(group$Zt)
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
  standata$t <- t
  standata$p <- as.array(p)
  standata$l <- as.array(l)
  standata$q <- ncol(Z)
  standata$len_theta_L <- sum(choose(p, 2), p)
  if (is_bernoulli) {
    stop("implement me!")
    parts0 <- extract_sparse_parts(Z[y == 0, , drop = FALSE])
    parts1 <- extract_sparse_parts(Z[y == 1, , drop = FALSE])
    standata$num_non_zero <- c(length(parts0$w), length(parts1$w))
    standata$w0 <- as.array(parts0$w)
    standata$w1 <- as.array(parts1$w)
    standata$v0 <- as.array(parts0$v - 1L)
    standata$v1 <- as.array(parts1$v - 1L)
    standata$u0 <- as.array(parts0$u - 1L)
    standata$u1 <- as.array(parts1$u - 1L)
  } else {
    parts <- extract_sparse_parts(Z)
    standata$num_non_zero <- length(parts$w)
    standata$w <- parts$w
    standata$v <- parts$v - 1L
    standata$u <- parts$u - 1L
  }
  standata$shape <- as.array(maybe_broadcast(decov$shape, t))
  standata$scale <- as.array(maybe_broadcast(decov$scale, t))
  standata$len_concentration <- sum(p[p > 1])
  standata$concentration <- 
    as.array(maybe_broadcast(decov$concentration, sum(p[p > 1])))
  standata$len_regularization <- sum(p > 1)
  standata$regularization <- 
    as.array(maybe_broadcast(decov$regularization, sum(p > 1)))
    
  standata$y <- y
  standata$weights <- weights
  standata$offset_ <- numeric(length(y))
  

  # call stan() to draw from posterior distribution
  if (is_continuous) {
    standata$ub_y <- Inf
    standata$lb_y <- if (is_gaussian) -Inf else 0
    standata$prior_scale_for_aux <- prior_scale_for_aux %ORifINF% 0
    standata$prior_df_for_aux <- c(prior_df_for_aux)
    standata$prior_mean_for_aux <- c(prior_mean_for_aux)
    standata$len_y <- length(y)
  } else if (is_bernoulli) {
    standata$prior_scale_for_aux <- 
      if (!length(group) || prior_scale_for_aux == Inf) 
        0 else prior_scale_for_aux
    standata$prior_mean_for_aux <- 0
    standata$prior_df_for_aux <- 0
    if (is_bernoulli) {
      stop("implement me!")
      y0 <- y == 0
      y1 <- y == 1
      standata$N <- c(sum(y0), sum(y1))
      standata$X0 <- array(xtemp[y0, , drop = FALSE], dim = c(1, sum(y0), ncol(xtemp)))
      standata$X1 <- array(xtemp[y1, , drop = FALSE], dim = c(1, sum(y1), ncol(xtemp)))
      standata$nnz_X0 = 0L 
      standata$w_X0 = double(0)
      standata$v_X0 = integer(0)
      standata$u_X0 = integer(0)
      standata$nnz_X1 = 0L 
      standata$w_X1 = double(0)
      standata$v_X1 = integer(0)
      standata$u_X1 = integer(0)
      if (length(weights) > 0L) {
        # nocov start
        # this code is unused because weights are interpreted as number of 
        # trials for binomial glms
        standata$weights0 <- weights[y0]
        standata$weights1 <- weights[y1]
        # nocov end
      } else {
        standata$weights0 <- double(0)
        standata$weights1 <- double(0)
      }
      if (length(offset) > 0L) {
        # nocov start
        standata$offset0 <- offset[y0]
        standata$offset1 <- offset[y1]
        # nocov end
      } else {
        standata$offset0 <- double(0)
        standata$offset1 <- double(0)
      }
      standata$K_smooth <- ncol(S)
      standata$S0 <- S[y0, , drop = FALSE]
      standata$S1 <- S[y1, , drop = FALSE]
      standata$smooth_map <- smooth_map
    }
  }
    
  prior_info <- summarize_glm_prior(
    user_prior_aux = prior_aux_stuff,
    user_prior_covariance = user_covariance,
    adjusted_prior_aux_scale = prior_scale_for_aux,
    family = family
  )
 
  for (varName in names(standata))
    if (is.logical(standata[[varName]])) standata[[varName]] <- as.integer(standata[[varName]])
  for (varName in c("len_theta_L"))
    standata[[varName]] <- as.integer(standata[[varName]])
  
  dotsList <- list(...)
  sampling <- getMethod("sampling", "stanmodel", asNamespace("rstan"))
  
  # extract defaults from "sampling()" and use them if necessary
  tryResult <- tryCatch(samplingFormals <- formals(body(sampling)[[2L]][[3L]]), error = function(e) e)
  if (is(tryResult, "error")) {
    warning("sampling function in rstan has been updated - using hard coded defaults; please notify package author")
    samplingFormals <- list(iter = 2000L, warmup = quote(iter %/% 2L), thin = 1L, chains = 4L, cores = quote(getOption("mc.cores", 1L)))
  }
  
  iter <- warmup <- thin <- chains <- NULL
  varNames <- c("iter", "warmup", "thin", "chains", "cores")
  for (varName in varNames) {
    assign(varName, if (!is.null(dotsList[[varName]])) dotsList[[varName]] else eval(samplingFormals[[varName]]))
    dotsList[[varName]] <- NULL
  } 
  
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
  
  verbose <- if (!is.null(dotsList[["verbose"]])) dotsList[["verbose"]] else 1L
  if (is.logical(verbose)) verbose <- as.integer(verbose)
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("'verbose' must an integer")
  verbose <- as.integer(verbose)
  
  refresh <- if (!is.null(dotsList[["refresh"]])) dotsList[["refresh"]] else NA_integer_
  if (!is.numeric(refresh) || length(refresh) != 1L || (!is.na(refresh) && refresh < 0))
    stop("'refresh' must be a non-negative integer")
  refresh <- as.integer(refresh)
  
  offset_type <- if (!is.null(dotsList[["offset_type"]])) dotsList[["offset_type"]] else "default"
  offset_type <- match.arg(offset_type, c("default", "fixef", "ranef"))
  offset_type <- which(offset_type == c("default", "fixef", "ranef")) - 1L
  
  bartData@sigma <- sigma_init
  bartControl <- dbarts::dbartsControl(n.chains = 1L, n.samples = 1L, n.burn = 0L,
                                       n.thin = thin.bart, n.threads = 1L,
                                       updateState = FALSE)
  
  bartData@n.cuts <- rep_len(attr(bartControl, "n.cuts"), ncol(bartData@x))
  bartControl@binary <- !is_continuous
  evalEnv <- sys.frame(sys.nframe())
  bartPriors <- dbarts:::parsePriors(bartControl, bartData, cgm, normal, fixed(1), evalEnv)
  bartModel <- new("dbartsModel", bartPriors$tree.prior, bartPriors$node.prior, bartPriors$resid.prior,
                    node.scale = if (!is_continuous) 3.0 else 0.5)
  
  # TODO: make more stan args user setable
  stan_args <- list(
    seed = sample.int(.Machine$integer.max, 1L),
    init_r = 2.0,
    thin = thin.stan,
    adapt_delta = adapt_delta)
  
  commonControl <- nlist(iter, warmup, verbose, refresh,
                         offset = offset, offset_type = offset_type,
                         ranef_init = ranef_init, sigma_init = sigma_init)
  
  
  chainResults <- vector("list", chains)
  if (cores <= 1L || chains <= 1L) {
    for (chainNum in seq_len(chains))
      chainResults[[chainNum]] <- mstan4bart_fitforreal(1L, bartControl, bartData, bartModel, standata, stan_args, commonControl, group)
  } else {
    if (commonControl$verbose > 0L)
      cat("starting multithreaded fit, futher output silenced\n")
    commonControl$verbose <- 0L
    cluster <- makeCluster(min(cores, chains))
    
    clusterExport(cluster, "mstan4bart_fitforreal", asNamespace("stan4bart"))
    clusterEvalQ(cluster, require(stan4bart))
    
    tryResult <- tryCatch(
      chainResults <- clusterMap(cluster, "mstan4bart_fitforreal", seq_len(chains), MoreArgs = nlist(bartControl, bartData, bartModel, standata, stan_args, commonControl, group)))
    
    stopCluster(cluster)
  }
        
  return(chainResults)

  # below here is garbage
  check <- try(check_stanfit(stanfit))
  if (!isTRUE(check)) return(standata)
  
  if (standata$len_theta_L > 0L) {
    thetas <- extract(stanfit, pars = "theta_L", inc_warmup = TRUE, 
                      permuted = FALSE)
    cnms <- group$cnms
    nc <- sapply(cnms, FUN = length)
    nms <- names(cnms)
    Sigma <- apply(thetas, 1:2, FUN = function(theta) {
      Sigma <- mkVarCorr(sc = 1, cnms, nc, theta, nms)
      unlist(sapply(Sigma, simplify = FALSE, 
                    FUN = function(x) x[lower.tri(x, TRUE)]))
    })
    l <- length(dim(Sigma))
    end <- tail(dim(Sigma), 1L)
    shift <- grep("^theta_L", names(stanfit@sim$samples[[1]]))[1] - 1L
    if (l == 3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
      stanfit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain] 
    }
    else for (chain in 1:end) {
      stanfit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
    }
    Sigma_nms <- lapply(cnms, FUN = function(grp) {
      nm <- outer(grp, grp, FUN = paste, sep = ",")
      nm[lower.tri(nm, diag = TRUE)]
    })
    for (j in seq_along(Sigma_nms)) {
      Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
    }
    Sigma_nms <- unlist(Sigma_nms)
  }
  new_names <- c(if (has_intercept) "(Intercept)", 
                  colnames(xtemp),
                  if (ncol(S)) colnames(S),
                  if (length(group) && length(group$flist)) c(paste0("b[", b_nms, "]")),
                  if (is_gaussian) "sigma", 
                  if (is_gamma) "shape", 
                  if (is_ig) "lambda",
                  if (is_nb) "reciprocal_dispersion",
                  if (is_beta) "(phi)",
                  if (ncol(S)) paste0("smooth_sd[", names(x)[-1], "]"),
                  if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
                  if (mean_PPD && !standata$clogit) "mean_PPD", 
                  "log-posterior")
  stanfit@sim$fnames_oi <- new_names
  return(structure(stanfit, prior.info = prior_info))
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
unpad_reTrms <- function(x, ...) UseMethod("unpad_reTrms")
unpad_reTrms.default <- function(x, ...) {
  if (is.matrix(x) || is.array(x))
    return(unpad_reTrms.array(x, ...))
  keep <- !grepl("_NEW_", names(x), fixed = TRUE)
  x[keep]
}

unpad_reTrms.array <- function(x, columns = TRUE, ...) {
  ndim <- length(dim(x))
  if (ndim > 3)
    stop("'x' should be a matrix or 3-D array")
  
  nms <- if (columns) 
    last_dimnames(x) else rownames(x)
  keep <- !grepl("_NEW_", nms, fixed = TRUE)
  if (length(dim(x)) == 2) {
    x_keep <- if (columns) 
      x[, keep, drop = FALSE] else x[keep, , drop = FALSE]
  } else {
    x_keep <- if (columns) 
      x[, , keep, drop = FALSE] else x[keep, , , drop = FALSE]
  }
  return(x_keep)
}

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


# Create "prior.info" attribute needed for prior_summary()
#
# @param user_* The user's prior, prior_intercept, prior_covariance, and 
#   prior_aux specifications. For prior and prior_intercept these should be
#   passed in after broadcasting the df/location/scale arguments if necessary.
# @param has_intercept T/F, does model have an intercept?
# @param has_predictors T/F, does model have predictors?
# @param adjusted_prior_*_scale adjusted scales computed if using autoscaled priors
# @param family Family object.
# @return A named list with components 'prior', 'prior_intercept', and possibly 
#   'prior_covariance' and 'prior_aux' each of which itself is a list
#   containing the needed values for prior_summary.
summarize_glm_prior <-
  function(user_prior_aux,
           user_prior_covariance,
           adjusted_prior_aux_scale,
           family) {
    rescaled_aux <- user_prior_aux$prior_autoscale_for_aux &&
      !is.na(user_prior_aux$prior_dist_name_for_aux) &&
      (user_prior_aux$prior_scale_for_aux != adjusted_prior_aux_scale)
    
    if (user_prior_aux$prior_dist_name_for_aux %in% "t") {
      if (all(user_prior_aux$prior_df_for_aux == 1)) {
        user_prior_aux$prior_dist_name_for_aux <- "cauchy"
      } else {
        user_prior_aux$prior_dist_name_for_aux <- "student_t"
      }
    }
    prior_list <- list()
    
    if (length(user_prior_covariance) > 0L)
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
