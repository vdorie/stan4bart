mstan4bart <- 
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           test = NULL,
           treatment = NULL,
           offset_test = NULL,
           verbose = FALSE,
           iter = 2000L,
           warmup = iter %/% 2L,
           thin = 1L,
           chains = 4L,
           cores = getOption("mc.cores", 1L),
           refresh = max(iter %/% 10L, 1L),
           offset_type = c("default", "fixef", "ranef", "bart", "parametric"),
           stan_args = NULL,
           bart_args = NULL)
{
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  if (!is.null(data)) {
    data <- as.data.frame(data)
    data <- drop_redundant_dims(data)
  }
  gl_call <- mc
  gl_call[[1L]] <- quoteInNamespace(glFormula)
  gl_call$control <- make_glmerControl()
  gl_call$data <- data
  gl_call$na.action <- na.action
  
  if (identical(na.action, stats::na.pass) || identical(na.action, "na.pass"))
    stop("na.action of 'na.pass' not allowed")
  
  for (name in setdiff(names(formals(mstan4bart)), names(formals(glFormula)))) {
    if (name %in% names(gl_call)) gl_call[[name]] <- NULL
  }
  
  glmod <- eval(gl_call, parent.frame())
  
  bartData <- glmod$bartData
  if ("b" %in% colnames(bartData@x)) {
    stop("mstan4bart does not allow the name 'b' for predictor variables", 
         call. = FALSE)
  }
  
  if (!has_outcome_variable(formula))
    stop("bart model requires a response variable")
  
  y <- model.response(glmod$fr)
  if (is.matrix(y)) {
    if (ncol(y) != 1L) stop("response variable must be a vector")
    y <- as.vector(y)
  }
  y <- as.double(y)
  if (!is.null(attr(glmod$fr, "na.action.all")))
    y <- y[setdiff(seq_len(nrow(glmod$fr)), attr(glmod$fr, "na.action.all"))]
  if (length(y) > 0L) {
    u.y <- unique(y)
    family <- if (length(u.y) == 2L && all(sort(u.y) == c(0, 1))) binomial(link = "probit") else gaussian()
  } else {
    stop("response required to fit bart model")
  }
  
  is_bernoulli <- family$family == "binomial"
  

  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(as.vector(model.weights(glmod$fr)))
  offset_type <- match.arg(offset_type)
  
  result <- nlist(y, weights, offset, offset_type = offset_type,
                  frame = glmod$fr, formula, na.action,
                  call = mc, family, bartData)
  if (!is.null(glmod$X)) {
    result$X <- glmod$X
    result$X_means <- apply(glmod$X, 2L, mean)
  }
  if (!is.null(glmod$reTrms)) {
    result$reTrms <- glmod$reTrms
  }
  class(result) <- "mstan4bartFit"
  
  if (!is.null(mc$test) && !is.null(mc$treatment))
    stop("only one of 'treatment' or 'test' can be specified")
  
  if (!is.null(mc$treatment)) {
    trtCall <- quote(getTreatmentData(result, treatment))
    trtCall[[1L]] <- quoteInNamespace(getTreatmentData)
    trtCall[[2L]] <- result
    trtCall[[3L]] <- mc$treatment
    
    trtResult <- eval(trtCall, parent.frame())
    result$treatment <- trtResult$treatment
    test <- trtResult$fr.test
  }
  
  if (!is.null(test)) {
    testDataFrames <- getTestDataFrames(result, test, na.action)
    result$test <- testDataFrames
    
    if (!is.null(offset_test)) result$test$offset <- offset_test
    
    result$bartData@x.test <- testDataFrames$X.bart
    result$bartData@testUsesRegularOffset <- FALSE
    
    result$test$X.bart <- NULL
  }
  
  if (iter == 0L)
    return(result)
  
  init_call <- mc
  
  bart_offset_init <- NULL
  sigma_init <- 1.0
  
  if (nzchar(system.file(package = "lme4"))) {
    if (family$family == "gaussian") {
      init_call[[1L]] <- quote(lme4::lmer)
    } else {
      init_call[[1L]] <- quote(lme4::glmer)
      init_call[["family"]] <- stats::binomial("probit")
    }
    
    formals_diff <-
      setdiff(names(formals(mstan4bart)), names(formals(eval(init_call[[1L]]))))
    for (name in formals_diff)
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    
    init_call$control <- lme4::lmerControl(check.conv.grad     = "ignore",
                                           check.conv.singular = "ignore",
                                           check.conv.hess     = "ignore")
    init_call$formula <- nobart(mc$formula)
    init_call$verbose <- FALSE
    init_call$na.action <- quote(stats::na.omit)
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!is(try_result, "error")) {
      bart_offset_init <- fitted(init_fit)
      if (!is_bernoulli)
        sigma_init <- sigma(init_fit)
    }
  }
  if (is.null(bart_offset_init)) {
    if (family$family == "gaussian") {
      init_call[[1L]] <- quote(stats::lm)
    } else {
      init_call[[1L]] <- quote(stats::glm)
      init_call[["family"]] <- stats::binomial("probit")
    }

    formals_diff <-
      setdiff(names(formals(mstan4bart)), names(formals(eval(init_call[[1L]]))))
    for (name in formals_diff)
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    
    init_call$control <- NULL
    init_call$formula <- subbars(nobart(mc$formula))
    init_call$verbose <- NULL
    init_call$na.action <- quote(stats::na.omit)
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!is(try_result, "error")) {
      bart_offset_init <- fitted(init_fit, type = "link")
      if (!is_bernoulli)
        sigma_init <- sigma(init_fit)
    }
  }
  if (is.null(bart_offset_init)) {
    init_call$formula <- nobars(nobart(mc$formula))
    init_fit <- eval(init_call, parent.frame())
    bart_offset_init <- fitted(init_fit, type = "link")
    if (!is_bernoulli)
      sigma_init <- sigma(init_fit)
  }
  
  # This can happen if parts of the formula got dropped during the init because we
  # drop down to lm/glm. In that case, the na patterns may no longer be the same
  # and the init vector will be too long.
  if (!is.null(bart_offset_init) && length(bart_offset_init) != length(y)) {
    if (inherits(init_fit, "merMod")) {
      na.action.init <- attr(init_fit@frame, "na.action")
    } else {
      na.action.init <- init_fit$na.action
    }
    init_rows <- seq_len(nrow(result$frame)) %not_in% na.action.init
    fit_rows <- seq_len(nrow(result$frame)) %not_in% attr(result$frame, "na.action.all")
    bart_offset_init <- bart_offset_init[fit_rows[init_rows]]
  }
   
  # Some trickery to allow calls like bart_args = list(k = chi(2, Inf))
  # to work. Since 'chi' is never defined, use a function that returns
  # a quoted version of the call.
  defn_env <- new.env(parent = parent.frame())
  defn_env$chi <- function(degreesOfFreedom = 1.25, scale = Inf)
    match.call()
  
  bart_args <- eval(mc[["bart_args"]], envir = defn_env)
  
  chain_results <- mstan4bart_fit(result,
                                  family,
                                  bart_offset_init,
                                  sigma_init,
                                  verbose,
                                  iter,
                                  warmup,
                                  thin,
                                  chains,
                                  cores,
                                  refresh,
                                  stan_args,
                                  bart_args)
  
  samples <- package_samples(chain_results, colnames(glmod$X), colnames(bartData@x))
  for (name in names(samples)) 
    result[[name]] <- samples[[name]]
  
  parameter_names <- rownames(chain_results[[1L]]$sample$stan$raw)
  # count the number of unconstrained parameters
  first_upar <- which.max(!endsWith(parameter_names, "__"))
  if ((aux_pos <- match("aux.1", parameter_names, 0L)) != 0L) {
    last_upar <- aux_pos - 1L
  } else if ((beta_pos <- match("beta.1", parameter_names, 0L)) != 0L) {
    last_upar <- beta_pos - 1L
  } else {
    last_upar <- match("b.1", parameter_names) - 1L
  }
  n_upars <- last_upar - first_upar + 1L

  if (as.integer(verbose) >= 0L)
    check_sampler_diagnostics(result, stan_args, n_upars)
  
  result
}

check_sampler_diagnostics <- function(object, stan_args, n_upars)
{
  n_d <- if ("divergent__" %in% dimnames(object$diagnostics)$diagnostic)
    sum(object$diagnostics["divergent__",,])
  else
    0L
  if (n_d > 0) {
    ad <- stan_args$adapt_delta %ORifNULL% 0.8
    
    warning("There were ", n_d, " divergent transitions after warmup. See\n",
            "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup\n", 
            "to find out why this is a problem and how to eliminate them.", call. = FALSE)
  }
  
  max_td <- stan_args$max_treedepth %ORifNULL% 10
  n_m <- if ("treedepth__" %in% dimnames(object$diagnostics)$diagnostic)
    sum(object$diagnostics["treedepth__",,] >= max_td)
  else
    0
  
  if (n_m > 0)
    warning("There were ", n_m,
            " transitions after warmup that exceeded the maximum treedepth.",
            " Increase max_treedepth above ", max_td, ". See\n",
            "http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded", call. = FALSE)
  
  
  n_e <- 0L
  if ("energy__" %in% dimnames(object$diagnostics)$diagnostic) {
    E <- as.matrix(object$diagnostics["energy__",,,drop = TRUE])
    threshold <- 0.2
    if (nrow(E) > 1) {
      EBFMI <- n_upars / apply(E, 2, var)
      n_e <- sum(EBFMI < threshold, na.rm = TRUE)
    }
    else n_e <- 0L
    if (n_e > 0)
      warning("There were ", n_e, 
              " chains where the estimated Bayesian Fraction of Missing Information",
              " was low. See\n", 
              "http://mc-stan.org/misc/warnings.html#bfmi-low", call. = FALSE)
  }
}

package_samples <- function(chain_results, fixef_names, bart_var_names) {
  result <- list()
  n_chains  <- length(chain_results)
  # must have a bart component, so we use that to determine the number of samples
  
  n_warmup <- 0L
  if (!is.null(chain_results[[1L]]$warmup) && !is.null(dim(chain_results[[1L]]$warmup$bart$train)))
    n_warmup <- dim(chain_results[[1L]]$warmup$bart$train)[2L]
  n_obs      <- dim(chain_results[[1L]]$sample$bart$train)[1L]
  n_obs_test <- 0L
  if (!is.null(chain_results[[1L]]$sample$bart$test))
    n_obs_test <- dim(chain_results[[1L]]$sample$bart$test)[1L]
  n_samples <- dim(chain_results[[1L]]$sample$bart$train)[2L]
  n_fixef   <- dim(chain_results[[1L]]$sample$stan$fixef)[1L]
  n_ranef_levels <- length(chain_results[[1L]]$sample$stan$ranef)
  if (n_ranef_levels > 0L) {
    n_ranef_at_level  <- sapply(chain_results[[1L]]$sample$stan$ranef, function(ranef_i) dim(ranef_i)[1L])
    n_groups_at_level <- sapply(chain_results[[1L]]$sample$stan$ranef, function(ranef_i) dim(ranef_i)[2L])
  }
  n_bart_vars <- dim(chain_results[[1L]]$sample$bart$varcount)[1L]
  aux_row <- which(rownames(chain_results[[1]]$sample$stan$raw) == "aux.1")
  
  result$bart_train <- array(sapply(seq_len(n_chains), function(i_chains)
                               chain_results[[i_chains]]$sample$bart$train),
                             dim = c(n_obs, n_samples, n_chains),
                             dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  if (n_obs_test > 0L) {
      result$bart_test <- array(sapply(seq_len(n_chains), function(i_chains)
                                  chain_results[[i_chains]]$sample$bart$test),
                                dim = c(n_obs_test, n_samples, n_chains),
                                dimnames = list(observation = NULL, sample = NULL, chain = NULL))
  }
  result$bart_varcount <- array(sapply(seq_len(n_chains), function(i_chains)
                               chain_results[[i_chains]]$sample$bart$varcount),
                             dim = c(n_bart_vars, n_samples, n_chains),
                             dimnames = list(predictor = bart_var_names, sample = NULL, chain = NULL))
  
  
  if (length(aux_row) > 0L) {
    result$sigma <- sapply(seq_len(n_chains), function(i_chain) {
      chain_results[[i_chain]]$sample$stan$raw[aux_row,]
    })
    dimnames(result$sigma) <- list(NULL)
    names(dimnames(result$sigma)) <- c("sample", "chain")
  }
  if (n_ranef_levels > 0L) {
    result$ranef <- lapply(seq_len(n_ranef_levels), function(i_level) {
      names <- dimnames(chain_results[[1L]]$sample$stan$ranef[[i_level]])
      names[4L] <- list(NULL)
      names(names) <- c("predictor", "group", "sample", "chain")
      array(sapply(seq_len(n_chains), function(i_chain)
              chain_results[[i_chain]]$sample$stan$ranef[[i_level]]),
            c(n_ranef_at_level[i_level], n_groups_at_level[i_level], n_samples, n_chains),
            dimnames = names)
    })
    names(result$ranef) <- names(chain_results[[1L]]$sample$stan$ranef)
    result$Sigma <- lapply(seq_len(n_ranef_levels), function(i_level) {
      names <- dimnames(chain_results[[1L]]$sample$stan$Sigma[[i_level]])
      names[4L] <- list(NULL)
      names(names)[3L:4L] <- c("sample", "chain")
      array(sapply(seq_len(n_chains), function(i_chain)
              chain_results[[i_chain]]$sample$stan$Sigma[[i_level]]),
            c(n_ranef_at_level[i_level], n_ranef_at_level[i_level], n_samples, n_chains),
            dimnames = names)
    })
    names(result$Sigma) <- names(chain_results[[1L]]$sample$stan$Sigma)
    #b_rows <- grep("^b\\.", rownames(chain_results[[1L]]$sample$stan$raw))
    #result$b <- array(sapply(seq_len(n_chains), function(i_chain)
    #                     chain_results[[i_chain]]$sample$stan$raw[b_rows,]),
    #                  c(sum(n_groups_at_level * n_ranef_at_level), n_samples, n_chains))
  }

  if (n_fixef > 0L) {
    names <- dimnames(chain_results[[1L]]$sample$stan$fixef)
    names[[1L]] <- fixef_names
    names[3L] <- list(NULL)
    names(names) <- c("predictor", "sample", "chain")
    
    result$fixef <- array(sapply(seq_len(n_chains), function(i_chain)
                            chain_results[[i_chain]]$sample$stan$fixef),
                          c(n_fixef, n_samples, n_chains),
                          dimnames = names)
    
  }
  
  if (!is.null(chain_results[[1L]]$sample$bart$k)) {
    result$k <- matrix(sapply(seq_len(n_chains), function(i_chains)
                              chain_results[[i_chains]]$sample$bart$k),
                       n_samples, n_chains,
                       dimnames = list(sample = NULL, chain = NULL))
  }
  
  diagnostic_names <- rownames(chain_results[[1]]$sample$stan$raw)
  diagnostic_names <- diagnostic_names[endsWith(diagnostic_names, "__")]
  result$diagnostics <- array(sapply(seq_len(n_chains), function(i_chain)
                                chain_results[[i_chain]]$sample$stan$raw[diagnostic_names,,drop = FALSE]),
                              c(length(diagnostic_names), n_samples, n_chains),
                              dimnames = list(diagnostic = diagnostic_names, sample = NULL, chain = NULL))  
  if (n_warmup > 0L) {
    result$warmup <- list()
    result$warmup$bart_train <- array(sapply(seq_len(n_chains), function(i_chains)
                                        chain_results[[i_chains]]$warmup$bart$train),
                                      dim = c(n_obs, n_warmup, n_chains),
                                      dimnames = list(observation = NULL, sample = NULL, chain = NULL))
    if (n_obs_test > 0L) {
      result$warmup$bart_test <- array(sapply(seq_len(n_chains), function(i_chains)
                                         chain_results[[i_chains]]$warmup$bart$test),
                                       dim = c(n_obs_test, n_warmup, n_chains),
                                       dimnames = list(observation = NULL, sample = NULL, chain = NULL))
    }
    result$warmup$bart_varcount <- array(sapply(seq_len(n_chains), function(i_chains)
                                           chain_results[[i_chains]]$warmup$bart$varcount),
                                         dim = c(n_bart_vars, n_warmup, n_chains),
                                         dimnames = list(predictor = bart_var_names, sample = NULL, chain = NULL))
    if (length(aux_row) > 0L) {
      result$warmup$sigma <- sapply(seq_len(n_chains), function(i_chain) {
        chain_results[[i_chain]]$warmup$stan$raw[aux_row,]
      })
      dimnames(result$warmup$sigma) <- list(NULL)
      names(dimnames(result$warmup$sigma)) <- c("sample", "chain")
    }
    if (n_ranef_levels > 0L) {
      result$warmup$ranef <- lapply(seq_len(n_ranef_levels), function(i_level) {
        names <- dimnames(chain_results[[1L]]$warmup$stan$ranef[[i_level]])
        names[4L] <- list(NULL)
        names(names) <- c("predictor", "group", "sample", "chain")
        array(sapply(seq_len(n_chains), function(j_chain)
                chain_results[[j_chain]]$warmup$stan$ranef[[i_level]]),
              c(n_ranef_at_level[i_level], n_groups_at_level[i_level], n_warmup, n_chains),
              dimnames = names)
      })
      names(result$warmup$ranef) <- names(chain_results[[1L]]$warmup$stan$ranef)
      result$warmup$Sigma <- lapply(seq_len(n_ranef_levels), function(i_level) {
        names <- dimnames(chain_results[[1L]]$warmup$stan$Sigma[[i_level]])
        names[4L] <- list(NULL)
        names(names)[3L:4L] <- c("sample", "chain")
        array(sapply(seq_len(n_chains), function(i_chain)
                chain_results[[i_chain]]$warmup$stan$Sigma[[i_level]]),
              c(n_ranef_at_level[i_level], n_ranef_at_level[i_level], n_warmup, n_chains),
              dimnames = names)
      })
      names(result$warmup$Sigma) <- names(chain_results[[1L]]$warmup$stan$Sigma)
    }
    if (n_fixef > 0L) {
      names <- dimnames(chain_results[[1L]]$warmup$stan$fixef)
      names[[1L]] <- fixef_names
      names[3L] <- list(NULL)
      names(names) <- c("predictor", "sample", "chain")
      
      result$warmup$fixef <- array(sapply(seq_len(n_chains), function(i_chain)
                                     chain_results[[i_chain]]$warmup$stan$fixef),
                                   c(n_fixef, n_warmup, n_chains),
                                   dimnames = names)
    }
    if (!is.null(chain_results[[1L]]$warmup$bart$k)) {
      result$warmup$k <- matrix(sapply(seq_len(n_chains), function(i_chains)
                                  chain_results[[i_chains]]$warmup$bart$k),
                                n_warmup, n_chains,
                                dimnames = list(sample = NULL, chain = NULL))
    }
  }
  
  if (!is.null(attr(chain_results, "sampler.bart"))) {
    result$sampler.bart <- attr(chain_results, "sampler.bart")
    result$range.bart <- matrix(sapply(seq_len(n_chains), function(i_chains)
                                       chain_results[[i_chains]]$range.bart),
                                nrow = 2L, ncol = n_chains,
                                dimnames = list(c("min", "max"), chain = NULL))
  }
  
  result
}

getTreatmentData <- function(object, treatment) {
  mc <- match.call()
  if (!is.null(mc$treatment)) {
    if (is.symbol(mc$treatment)) treatment <- as.character(mc$treatment)
    if (!is.character(treatment))
      stop("'treament' argument must be a character or symbol")
  }
  
  if (treatment %not_in% colnames(object$frame))
      stop("treament must be the name of a column in data")
  uq <- sort(unique(object$frame[[treatment]]))
  if (length(uq) != 2L || !all(uq == c(0, 1)))
    stop("treatment must in { 0, 1 }^n")
  fr.test <- object$frame
  fr.test[[treatment]] <- 1 - fr.test[[treatment]]
  
  nlist(fr.test, treatment)
}

