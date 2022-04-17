stan4bart <- 
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
           skip = 1L,
           chains = 4L,
           cores = getOption("mc.cores", 1L),
           refresh = max(iter %/% 10L, 1L),
           offset_type = c("default", "fixef", "ranef", "bart", "parametric"),
           seed = NA_integer_,
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
  
  for (name in setdiff(names(formals(stan4bart)), names(formals(glFormula)))) {
    if (name %in% names(gl_call)) gl_call[[name]] <- NULL
  }
  
  glmod <- eval(gl_call, parent.frame())
  
  bartData <- glmod$bartData
  if ("b" %in% colnames(bartData@x)) {
    stop("stan4bart does not allow the name 'b' for predictor variables", 
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
  class(result) <- "stan4bartFit"
  
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
      setdiff(names(formals(stan4bart)), names(formals(eval(init_call[[1L]]))))
    for (name in formals_diff)
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    
    init_call$control <- lme4::lmerControl(check.conv.grad     = "ignore",
                                           check.conv.singular = "ignore",
                                           check.conv.hess     = "ignore")
    init_call$formula <- nobart(mc$formula)
    init_call$verbose <- FALSE
    init_call$na.action <- quote(stats::na.omit)
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!inherits(try_result, "error")) {
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
      setdiff(names(formals(stan4bart)), names(formals(eval(init_call[[1L]]))))
    for (name in formals_diff)
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    
    init_call$control <- NULL
    init_call$formula <- subbars(nobart(mc$formula))
    init_call$verbose <- NULL
    init_call$na.action <- quote(stats::na.omit)
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!inherits(try_result, "error")) {
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
  
  chain_results <- stan4bart_fit(result,
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
  
  samples <- package_samples(chain_results, colnames(bartData@x))
  for (name in names(samples)) 
    result[[name]] <- samples[[name]]
  
  stan_par_names <- dimnames(result$stan)[[1L]]
  diagnostic_names <- stan_par_names[endsWith(stan_par_names, "__")]
  
  upar_names <- stan_par_names[grepl("^(?:gamma|z_beta|global|local|caux|mix|one_over_lambda|z_b|z_T|rho|zeta|tau|aux_unscaled)\\.", stan_par_names, perl = TRUE)]
  tpar_names <- setdiff(stan_par_names, c(diagnostic_names, upar_names))
  
  result$par_names <- list(diagnostic = diagnostic_names, upar = upar_names, tpar = tpar_names)
  
  if (as.integer(verbose) >= 0L)
    check_sampler_diagnostics(result, stan_args, length(upar_names))
  
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

package_samples <- function(chain_results, bart_var_names) {
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
  
  n_bart_vars <- dim(chain_results[[1L]]$sample$bart$varcount)[1L]
  
  chain_names <- paste0("chain:", seq_len(n_chains))
  
  # grab the bart bits
  result$bart_train <- array(sapply(seq_len(n_chains), function(i_chains)
                               chain_results[[i_chains]]$sample$bart$train),
                             dim = c(n_obs, n_samples, n_chains),
                             dimnames = list(observation = NULL, iterations = NULL, chain = chain_names))
  
  if (n_obs_test > 0L) {
    result$bart_test <- array(sapply(seq_len(n_chains), function(i_chains)
                                chain_results[[i_chains]]$sample$bart$test),
                              dim = c(n_obs_test, n_samples, n_chains),
                              dimnames = list(observation = NULL, iterations = NULL, chain = chain_names))
  }
  result$bart_varcount <- array(sapply(seq_len(n_chains), function(i_chains)
                               chain_results[[i_chains]]$sample$bart$varcount),
                             dim = c(n_bart_vars, n_samples, n_chains),
                             dimnames = list(predictor = bart_var_names, iterations = NULL, chain = chain_names))
  
  if (!is.null(chain_results[[1L]]$sample$bart$k)) {
    result$k <- matrix(sapply(seq_len(n_chains), function(i_chains)
                              chain_results[[i_chains]]$sample$bart$k),
                       n_samples, n_chains,
                       dimnames = list(iterations = NULL, chain = chain_names))
  }

  # treat stan as a big array with no processing, until someone asks
  stan_par_names <- dimnames(chain_results[[1L]]$sample$stan)[[1L]]
  n_stan_pars <- length(stan_par_names)
  
  result$stan <- array(sapply(seq_len(n_chains), function(i_chains)
                              chain_results[[i_chains]]$sample$stan),
                       dim = c(n_stan_pars, n_samples, n_chains),
                       dimnames = list(parameters = stan_par_names, iterations = NULL, chain = chain_names))
  
  # Do it all again if there are warmup samples
  if (n_warmup > 0L) {
    warmup <- list()
    
    warmup$bart_train <- array(sapply(seq_len(n_chains), function(i_chains)
                                 chain_results[[i_chains]]$warmup$bart$train),
                               dim = c(n_obs, n_warmup, n_chains),
                               dimnames = list(observation = NULL, iterations = NULL, chain = chain_names))
    if (n_obs_test > 0L) {
      warmup$bart_test <- array(sapply(seq_len(n_chains), function(i_chains)
                                  chain_results[[i_chains]]$warmup$bart$test),
                                dim = c(n_obs_test, n_warmup, n_chains),
                                dimnames = list(observation = NULL, iterations = NULL, chain = chain_names))
    }
    warmup$bart_varcount <- array(sapply(seq_len(n_chains), function(i_chains)
                                    chain_results[[i_chains]]$warmup$bart$varcount),
                                  dim = c(n_bart_vars, n_warmup, n_chains),
                                  dimnames = list(predictor = bart_var_names, iterations = NULL, chain = chain_names))
    
    if (!is.null(chain_results[[1L]]$warmup$bart$k)) {
      warmup$k <- matrix(sapply(seq_len(n_chains), function(i_chains)
                           chain_results[[i_chains]]$warmup$bart$k),
                         n_warmup, n_chains,
                         dimnames = list(iterations = NULL, chain = chain_names))
    }
    
    warmup$stan <- array(sapply(seq_len(n_chains), function(i_chains)
                          chain_results[[i_chains]]$warmup$stan),
                        dim = c(n_stan_pars, n_warmup, n_chains),
                        dimnames = list(parameters = stan_par_names, iterations = NULL, chain = chain_names))
    
    result$warmup <- warmup
  }
  
  if (!is.null(attr(chain_results, "sampler.bart"))) {
    result$sampler.bart <- attr(chain_results, "sampler.bart")
    result$range.bart <- matrix(sapply(seq_len(n_chains), function(i_chains)
                                       chain_results[[i_chains]]$range.bart),
                                nrow = 2L, ncol = n_chains,
                                dimnames = list(c("min", "max"), chain = chain_names))
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
  if (is.logical(fr.test[[treatment]])) {
    fr.test[[treatment]] <- !fr.test[[treatment]]
  } else {
    fr.test[[treatment]] <- 1 - fr.test[[treatment]]
  }
  
  nlist(fr.test, treatment)
}

