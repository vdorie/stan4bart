if (FALSE) getSigma <- function(cnms, samples) {
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

if (FALSE) getFixef <- function(samples)
  samples[grep("^(?:beta|gamma)\\.", rownames(samples), perl = TRUE),,drop = FALSE]

if (FALSE) getRanef <- function(group, samples) {
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
mstan4bart_fit_worker <- function(chain.num, seed, control.bart, data.bart, model.bart, data.stan, control.stan, control.common, group)
{
  if (!is.na(seed))
    set.seed(seed)
  
  control.stan$seed <- sample.int(.Machine$integer.max, 1L)
  if (is.na(seed))
    orig_seed <- .Random.seed
  
  sampler <- .Call(C_stan4bart_create, control.bart, data.bart, model.bart, data.stan, control.stan, control.common)
  if (control.common$verbose > 0L) {
    cat("fitting chain ", chain.num, "\n", sep = "")
    .Call(C_stan4bart_printInitialSummary, sampler)
  }
  results <- list()
  if (control.common$warmup > 0L)
    results$warmup  <- .Call(C_stan4bart_run, sampler, control.common$warmup, TRUE, "both")
  .Call(C_stan4bart_disengageAdaptation, sampler)
  results$sample <- .Call(C_stan4bart_run, sampler, control.common$iter - control.common$warmup,
                          FALSE, "both")
  
  if (control.bart@keepTrees) {
    results$state.bart <- .Call(C_stan4bart_exportBARTState, sampler)
    results$range.bart <- .Call(C_stan4bart_getBARTDataRange, sampler)
  }
  
  results
}

mstan4bart_fit <- 
  function(object,
           family,
           bart_offset_init,
           sigma_init,
           verbose,
           iter,
           warmup,
           skip,
           chains,
           cores,
           refresh,
           seed,
           stan_args,
           bart_args)
{  
  matched_call <- match.call()
  
  supported_families <- c("binomial", "gaussian")
  fam <- which(pmatch(supported_families, family$family, nomatch = 0L) == 1L)
  
  # Use the rstanarm code to be consistent about ordering, but enforce
  # explicit constraint on allowable families.
  supported_links <- supported_glm_links(supported_families[fam])
  link <- which(supported_links == family$link)
  #if (!length(link)) 
  #  stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  if (family$family == "gaussian" && family$link != "identity")
    stop("'link' must be 'identity' if family is 'gaussian'") 
  if (family$family == "binomial" && family$link != "probit")
    stop("'link' must be 'probit' if family is 'binomial'") 
  
  
  if (is.null(stan_args))
    stan_args <- list()
  if (is.null(stan_args[["prior_covariance"]]))
    stan_args$prior_covariance <- decov()
  decov <- stan_args[["prior_covariance"]]
  if (is.null(stan_args[["prior"]]))
    stan_args$prior <- default_prior_coef_gaussian()
  #if (is.null(stan_args[["prior_intercept"]]))
  #  stan_args$prior_intercept <- default_prior_intercept_gaussian()
  if (!is.null(stan_args[["prior_intercept"]]))
    warning("intercepts for BART models are redundantly parameterized and not included")
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
  xtemp <- xbar <- has_intercept <- NULL # R CMD check
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
    # link = family$link,
    # mstan4bart: since we use latent variables, this is essentially a
    # gaussian/identity problem
    link = "identity", 
    ok_dists = ok_dists
  )
  # prior_{dist, mean, scale, df, dist_name, autoscale}, 
  # global_prior_df, global_prior_scale, slab_df, slab_scale
  for (name in names(prior_stuff))
    stan_args[[name]] <- prior_stuff[[name]]
  
  if (isTRUE(is.list(stan_args$prior_intercept)) && 
      isTRUE(stan_args$prior_intercept$default)) 
  {
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
    # link = family$link,
    # mstan4bart: since we use latent variables, this is essentially a
    # gaussian/identity problem
    link = "identity",
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
    needs_intercept <- (!is_gaussian && linkname == "identity") ||
      # is_gamma && linkname == "inverse" ||
      (is.binomial(famname) && linkname == "log")
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
  
  if (is.null(stan_args[["QR"]]))
    stan_args[["QR"]] <- FALSE
  if (!stan_args$QR && stan_args$prior_dist > 0L && stan_args$prior_autoscale) {
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
  
  if (stan_args$QR && ncol(xtemp) > 0L) {
    if (ncol(xtemp) <= 1)
      stop("'QR' can only be specified when there are multiple predictors.")
    #if (sparse)
    #  stop("'QR' and 'sparse' cannot both be TRUE.")
    cn <- colnames(xtemp)
    decomposition <- qr(xtemp)
    Q <- qr.Q(decomposition)
    if (stan_args$prior_autoscale) scale_factor <- sqrt(nrow(xtemp) - 1L)
    else scale_factor <- diag(qr.R(decomposition))[ncol(xtemp)]
    R_inv <- qr.solve(decomposition, Q) * scale_factor
    xtemp <- Q * scale_factor
    colnames(xtemp) <- cn
    xbar <- c(xbar %*% R_inv)
  }
    
  if (length(weights) > 0L && all(weights == 1)) weights <- double()
  if (length(offset)  > 0L && all(offset  == 0)) offset  <- double()
  
  # create entries in the data block of the .stan file
  data.stan <- with(stan_args, nlist(
    N = length(y),
    K = ncol(xtemp),
    has_weights = length(weights) > 0L,
    prior_dist_for_aux = prior_dist_for_aux,
    #family = stan_family_number(famname),  # hard-coded in model
    #link,                                  # hard-coded in model
    has_intercept,
    is_binary = is_bernoulli,
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
    
    parts <- extract_sparse_parts(Z)
    data.stan$num_non_zero <- length(parts$w)
    data.stan$w <- parts$w
    data.stan$v <- parts$v - 1L
    data.stan$u <- parts$u - 1L
    
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
    
    data.stan$num_non_zero <- 0L
    data.stan$w <- double(0)
    data.stan$v <- integer(0)
    data.stan$u <- rep.int(0L, length(y) + 1L)
    
    data.stan$shape <- data.stan$scale <- data.stan$concentration <-
      data.stan$regularization <- rep(0, 0)
    data.stan$len_concentration <- 0L
    data.stan$len_regularization <- 0L
  }
  
  data.stan$X <- xtemp
  data.stan$y <- y
  data.stan$weights <- weights
  data.stan$offset_ <- numeric(length(y))
  

  if (is_continuous) {
    data.stan$ub_y <- Inf
    data.stan$lb_y <- if (is_gaussian) -Inf else 0
    data.stan$prior_scale_for_aux <- stan_args$prior_scale_for_aux %ORifINF% 0
    data.stan$prior_df_for_aux <- c(stan_args$prior_df_for_aux)
    data.stan$prior_mean_for_aux <- c(stan_args$prior_mean_for_aux)
    data.stan$len_y <- length(y)
  } else if (is_bernoulli) {
    data.stan$ub_y <- Inf
    data.stan$lb_y <- -Inf
    data.stan$prior_scale_for_aux <- 0
    data.stan$prior_mean_for_aux <- 0
    data.stan$prior_df_for_aux <- 0
    data.stan$len_y <- length(y)
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
  iter <- as.integer(iter)
  if (!is.numeric(warmup) || length(warmup) != 1L || warmup < 0)
    stop("'warmup' must be a non-negative integer")
  warmup <- as.integer(warmup)
  if (!is.numeric(skip) || length(skip) < 1L || length(skip) > 2L || any(skip <= 0))
    stop("'skip' must be one or two positive integers")
  if (!is.null(names(skip)) && any(c("bart", "stan") %in% names(skip))) {
    skip.bart <- if ("bart" %in% names(skip)) skip[["bart"]] else 1L
    skip.stan <- if ("stan" %in% names(skip)) skip[["stan"]] else 1L
  } else {
    skip.bart <- skip[1L]
    skip.stan <- if (length(skip) > 1L) skip[2L] else skip[1L]
  }
  skip.bart <- as.integer(skip.bart)
  skip.stan <- as.integer(skip.stan)
  
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
  
  offset_type <- which(offset_type == c("default", "fixef", "ranef", "bart", "parametric")) - 1L
  
  # have to make sure that end node priors that are constructed as in chi(df, scale)
  # work correctly
  if (is.null(bart_args)) bart_args <- list()
  data.bart@sigma <- sigma_init
  control_call <- quote(dbarts::dbartsControl(n.chains = 1L, n.samples = 1L, n.burn = 0L,
                                              n.thin = skip.bart, n.threads = 1L,
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
  
  cgm <- normal <- fixed <- parsePriors <- NULL # R CMD check
  prior_call <- quote(parsePriors(control.bart, data.bart, cgm, normal, fixed(1), evalEnv))
  prior_call[[1L]] <- quoteInNamespace(parsePriors)
  prior_call[[1L]][[2L]] <- quote(dbarts)
  
  if (!is.null(bart_args[["k"]])) {
    end_node_pos <- which(as.character(prior_call) == "normal")
    end_node_prior <- quote(normal(k = k))
    end_node_prior[[2L]] <- bart_args[["k"]]
    prior_call[[end_node_pos]] <- end_node_prior
  }
  bart_priors <- eval(prior_call)
  model.bart <- new("dbartsModel", bart_priors$tree.prior, bart_priors$node.prior,
                    bart_priors$node.hyperprior, bart_priors$resid.prior,
                    node.scale = if (!is_continuous) 3.0 else 0.5)
  
  
  control.stan <- list(
    init_r = stan_args[["init_r"]] %ORifNULL% 2.0,
    skip = skip.stan,
    adapt_gamma = stan_args[["adapt_gamma"]],
    adapt_delta = stan_args[["adapt_delta"]],
    adapt_kappa = stan_args[["adapt_kappa"]]
  )
  
  control.common <- nlist(iter, warmup, verbose, refresh, is_binary = is_bernoulli,
                          offset = offset, offset_type = offset_type,
                          bart_offset_init = bart_offset_init, sigma_init = sigma_init)
  
  chainResults <- vector("list", chains)
  runSingleThreaded <- cores <= 1L || chains <= 1L
  if (!runSingleThreaded) {
    tryResult <- tryCatch(cluster <- makeCluster(min(cores, chains), "PSOCK"), error = function(e) e)
    if (inherits(tryResult, "error"))
      tryResult <- tryCatch(cluster <- makeCluster(min(cores, chains), "FORK"), error = function(e) e)
    
    if (inherits(tryResult, "error")) {
      warning("unable to multithread, defaulting to single: ", tryResult$message)
      runSingleThreaded <- TRUE
    } else {
      if (control.common$verbose > 0L)
        cat("starting multithreaded fit, futher output silenced\n")
      control.common$verbose <- -1L
      
      if (!is.na(seed)) {
        # We draw sequentially from the given seed, one for each thread. To be polite
        # (more to match bart), we set the seed back when we're are done.
        oldSeed <- .GlobalEnv[[".Random.seed"]]
        
        set.seed(seed)
        randomSeeds <- sample.int(.Machine$integer.max, chains)
        
        if (!is.null(oldSeed))
          .Random.seed <- oldSeed
      } else {
        randomSeeds <- rep.int(NA_integer_, chains)
      }
      
      clusterExport(cluster, "mstan4bart_fit_worker", asNamespace("stan4bart"))
      clusterEvalQ(cluster, require(stan4bart))
      
      tryResult <- tryCatch(
        chainResults <- clusterMap(
          cluster, "mstan4bart_fit_worker",
          seq_len(chains), randomSeeds,
          MoreArgs = nlist(control.bart, data.bart, model.bart, data.stan,
                           control.stan, control.common, group)),
        error = function(e) e)
    
      stopCluster(cluster)
      
      if (inherits(tryResult, "error")) {
        warning("error running multithreaded, defaulting to single: ", tryResult$message)
        runSingleThreaded <- TRUE
      }
    }
  }
  
  if (runSingleThreaded) {
    if (!is.na(seed)) {
      # If the seed was passed in, since we're running single threaded everything will draw
      # from the built-in generator. In that case, we just have to set.seed and set it
      # back when done.
      oldSeed <- .GlobalEnv[[".Random.seed"]]
      set.seed(seed)
    }
    
    for (chainNum in seq_len(chains))
      chainResults[[chainNum]] <- mstan4bart_fit_worker(chainNum, NA_integer_, control.bart, data.bart, model.bart, data.stan, control.stan, control.common, group)
    
    if (exists("oldSeed"))
      .Random.seed <- oldSeed
  }
  
  if (stan_args$QR) {
    fixef_rows <- dimnames(chainResults[[1L]]$sample$stan)[[1L]]
    fixef_rows <- startsWith(fixef_rows, "beta.")
    if (any(fixef_rows) && exists("R_inv")) for (chainNum in seq_len(chains)) {
      if (!is.null(chainResults[[chainNum]]$warmup))
        chainResults[[chainNum]]$warmup$stan[fixef_rows,] <-
          R_inv %*% chainResults[[chainNum]]$warmup$stan[fixef_rows,]
      chainResults[[chainNum]]$sample$stan[fixef_rows,] <-
        R_inv %*% chainResults[[chainNum]]$sample$stan[fixef_rows,]
    }
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



