mstan4bart <- 
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           test = NULL,
           offset_test = NULL,
           verbose = FALSE,
           iter = 2000L,
           warmup = iter %/% 2L,
           thin = 1L,
           chains = 4L,
           cores = getOption("mc.cores", 1L),
           refresh = max(iter %/% 10L, 1L),
           treatment = NULL,
           prior_covariance = decov(),
           offset_type = c("default", "fixef", "ranef"),
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
  
  for (name in setdiff(names(formals(mstan4bart)), names(formals(glFormula)))) {
    if (name %in% names(gl_call)) gl_call[[name]] <- NULL
  }
  
  glmod <- eval(gl_call, parent.frame())
  family <- glmod$family
  
  bartData <- glmod$bartData
  if ("b" %in% colnames(bartData@x)) {
    stop("mstan4bart does not allow the name 'b' for predictor variables", 
         call. = FALSE)
  }
  
  if (!has_outcome_variable(formula))
    stop("bart model requires a response variable")
  y <- glmod$fr[, as.character(glmod$formula[2L])]  
  if (is.matrix(y)) {
    if (ncol(y) != 1L) stop("response variable must be a vector")
    y <- as.vector(y)
  }

  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(as.vector(model.weights(glmod$fr)))
  offset_type <- match.arg(offset_type)
  
  if (is.null(prior_covariance))
    stop("'prior_covariance' cannot be NULL", call. = FALSE)
  group <- glmod$reTrms
  group$decov <- prior_covariance
    
  init_call <- mc
    
  bart_offset_init <- NULL
  sigma_init <- NULL
  
  if (nzchar(system.file(package = "lme4"))) {
    init_call[[1L]] <- quote(lme4::lmer)
    for (name in setdiff(names(formals(mstan4bart)), names(formals(lme4::lmer)))) {
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    }

    init_call$control <- lme4::lmerControl()
    init_call$formula <- nobart(mc$formula)
    init_call$verbose <- FALSE
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!is(try_result, "error")) {
      bart_offset_init <- fitted(init_fit)
      sigma_init <- sigma(init_fit)
    }
  }
  if (is.null(bart_offset_init)) {
    init_call[[1L]] <- quote(stats::lm)
    for (name in setdiff(names(formals(mstan4bart)), names(formals(stats::lm)))) {
      if (name %in% names(init_call)) init_call[[name]] <- NULL
    }
    init_call$control <- NULL
    init_call$formula <- subbars(nobart(mc$formula))
    init_call$verbose <- NULL
    try_result <- tryCatch(init_fit <- suppressWarnings(eval(init_call, parent.frame())), error = function(e) e)
    if (!is(try_result, "error")) {
      bart_offset_init <- fitted(init_fit)
      sigma_init <- sigma(init_fit)
    }
  }
  if (is.null(bart_offset_init)) {
    init_call$formula <- nobars(nobart(mc$formula))
    init_fit <- eval(init_call, parent.frame())
    bart_offset_init <- fitted(init_fit)
    sigma_init <- sigma(init_fit)
  }
  
  chain_results <- mstan4bart_fit(bartData, glmod$X, y,
                                  weights, offset,
                                  family,
                                  group,
                                  bart_offset_init,
                                  sigma_init,
                                  verbose,
                                  iter,
                                  warmup,
                                  thin,
                                  chains,
                                  cores,
                                  refresh,
                                  offset_type,
                                  stan_args, bart_args)
  
  result <- package_samples(chain_results, colnames(glmod$X), colnames(bartData@x))
  
  if (!is.null(glmod$X)) {
    result$X <- glmod$X
    result$X_means <- apply(glmod$X, 2L, mean)
  }
  if (!is.null(glmod$X.test))      result$X.test <- glmod$X.test
  if (!is.null(glmod$reTrms)) {
    result$reTrms <- glmod$reTrms
    if (!is.null(glmod$reTrms.test)) result$reTrms$Zt.test <- glmod$reTrms.test$Zt
  }
  if (!is.null(offset_test))       result$offset_test <- offset_test
  if (!is.null(glmod$treatment))   result$treatment <- glmod$treatment
  
  result$formula   <- formula
  result$terms     <- glmod$terms
  result$na.action <- na.action
  result$call      <- mc
  result$frame     <- glmod$fr
  
  class(result) <- "mstan4bartFit"
  result
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
  aux_row <- which(rownames(chain_results[[1]]$sample$stan$raw) == "aux")
  
  result$bart_train <- array(sapply(seq_along(n_chains), function(i_chains)
                               chain_results[[i_chains]]$sample$bart$train),
                             dim = c(n_obs, n_samples, n_chains),
                             dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
  if (n_obs_test > 0L) {
      result$bart_test <- array(sapply(seq_along(n_chains), function(i_chains)
                                  chain_results[[i_chains]]$sample$bart$test),
                                dim = c(n_obs_test, n_samples, n_chains),
                                dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
  }
  result$bart_varcount <- array(sapply(seq_along(n_chains), function(i_chains)
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
  if (n_warmup > 0L) {
    result$warmup <- list()
    result$warmup$bart_train <- array(sapply(seq_along(n_chains), function(i_chains)
                                        chain_results[[i_chains]]$warmup$bart$train),
                                      dim = c(n_obs, n_warmup, n_chains),
                                      dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
    if (n_obs_test > 0L) {
      result$warmup$bart_test <- array(sapply(seq_along(n_chains), function(i_chains)
                                         chain_results[[i_chains]]$warmup$bart$test),
                                       dim = c(n_obs_test, n_warmup, n_chains),
                                       dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
    }
    result$warmup$bart_varcount <- array(sapply(seq_along(n_chains), function(i_chains)
                                           chain_results[[i_chains]]$sample$bart$varcount),
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
  }
  
  if (!is.null(attr(chain_results, "sampler.bart")))
    result$sampler.bart <- attr(chain_results, "sampler.bart")
  
  
  # TODO: turn into test
  # all(as.vector(result$ranef[[1]][,,,1]) == as.vector(chain_results[[1L]]$sample$stan$ranef[[1]]))
  # all(as.vector(result$ranef[[1]][,,,2]) == as.vector(chain_results[[2L]]$sample$stan$ranef[[1]]))
  # all(as.vector(result$Sigma[[2]][,,,1]) == as.vector(chain_results[[1]]$sample$stan$Sigma[[2]]))
  # all(as.vector(result$Sigma[[2]][,,,2]) == as.vector(chain_results[[2]]$sample$stan$Sigma[[2]]))
  
  result
}


default_prior_intercept_gaussian = function() {
  out <- normal(0, 2.5, autoscale = TRUE)
  out$location <- NULL # not determined yet
  out$default <- TRUE
  out$version <- utils::packageVersion("stan4bart")
  out
}

default_prior_coef_gaussian = function() {
  out <- normal(0, 2.5, autoscale = TRUE)
  out$default <- TRUE
  out$version <- utils::packageVersion("stan4bart")
  out
}

normal <- function(location = 0, scale = NULL, autoscale = FALSE) {
  validate_parameter_value(scale)
  nlist(dist = "normal", df = NA, location, scale, autoscale)
}

student_t <- function(df = 1, location = 0, scale = NULL, autoscale = FALSE) {
  validate_parameter_value(scale)
  validate_parameter_value(df)
  nlist(dist = "t", df, location, scale, autoscale)
}

cauchy <- function(location = 0, scale = NULL, autoscale = FALSE) {
  student_t(df = 1, location = location, scale = scale, autoscale)
}

hs <- function(df = 1, global_df = 1, global_scale = 0.01,
               slab_df = 4, slab_scale = 2.5) {
  validate_parameter_value(df)
  validate_parameter_value(global_df)
  validate_parameter_value(global_scale)
  validate_parameter_value(slab_df)
  validate_parameter_value(slab_scale)
  nlist(dist = "hs", df, location = 0, scale = 1, 
        global_df, global_scale, slab_df, slab_scale)
}

hs_plus <- function(df1 = 1, df2 = 1, global_df = 1, global_scale = 0.01,
                    slab_df = 4, slab_scale = 2.5) {
  validate_parameter_value(df1)
  validate_parameter_value(df2)
  validate_parameter_value(global_df)
  validate_parameter_value(global_scale)
  validate_parameter_value(slab_df)
  validate_parameter_value(slab_scale)
  # scale gets used as a second df hyperparameter
  nlist(dist = "hs_plus", df = df1, location = 0, scale = df2, global_df, 
        global_scale, slab_df, slab_scale)
}


laplace <- function(location = 0, scale = NULL, autoscale = FALSE) {
  nlist(dist = "laplace", df = NA, location, scale, autoscale)
}

lasso <- function(df = 1, location = 0, scale = NULL, autoscale = FALSE) {
  nlist(dist = "lasso", df, location, scale, autoscale)
}

exponential <- function(rate = 1, autoscale = TRUE) {
  stopifnot(length(rate) == 1)
  validate_parameter_value(rate)
  nlist(dist = "exponential", df = NA, location = NA, scale = 1 / rate, 
        autoscale)
}

validate_parameter_value <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x)) 
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0)) 
        stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}

