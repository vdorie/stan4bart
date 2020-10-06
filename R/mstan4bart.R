 # results are:
  # named list of length = num chains
  #   warmup - has same elements as sample
  #   sample
  #     stan
  #       raw   - matrix of num pars x num samples
  #       Sigma - list of length num grouping factors
  #         g.XX - array of cov rows x cov cols x num samples
  #       ranef - list of length num groupin factors
  #         g.XX - array of num ranef @ factor x num groups x num samples
  #     bart
  #       sigma - length of num samples (use aux in stan$raw instead)
  #       train - num obs x num samples
  #       test  - num obs x num samples, contains counterfactuals
  #       varcount - num predictors x num samples
  # also present are two attributes:
  #   Zt.obs - random effects design matrix for observed
  #   Zt.cf  - random effects design matrix for counterfactuals

mstan4bart <- 
  function(formula,
           data = NULL,
           subset,
           weights,
           na.action = getOption("na.action", "na.omit"),
           offset,
           contrasts = NULL,
           test = NULL,
           test.offset = NULL,
           treatment = NULL,
           ...,
           prior = default_prior_coef_gaussian(),
           prior_intercept = default_prior_intercept_gaussian(),
           prior_aux = exponential(),
           prior_covariance = decov(),
           adapt_delta = NULL)
{
  
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  if (!is.null(data)) {
    data <- as.data.frame(data)
    data <- drop_redundant_dims(data)
  }
  mc[[1]] <- quoteInNamespace(glFormula)
  mc$control <- make_glmerControl()
  mc$data <- data
  mc$prior_covariance <- mc$prior_aux <-
    mc$adapt_delta <- mc$... <- NULL
  
  glmod <- eval(mc, parent.frame())
  family <- glmod$family
  
  bartData <- glmod$bartData
  if ("b" %in% colnames(bartData@x)) {
    stop("mstan4bart does not allow the name 'b' for predictor variables.", 
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
  
  if (is.null(prior_aux)) 
    prior_aux <- list()
  if (is.null(prior_covariance))
    stop("'prior_covariance' can't be NULL.", call. = FALSE)
  group <- glmod$reTrms
  group$decov <- prior_covariance
    
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- quote(lme4::lmer)
  mc$control <- lme4::lmerControl()
  mc$formula <- nobart(mc$formula)
  mc$data <- data
  mc$verbose <- FALSE
  mc$prior_covariance <- mc$prior_aux <-
    mc$adapt_delta <- mc$treatment <- mc$... <- NULL
  lmerFit <- suppressWarnings(eval(mc, parent.frame()))
  chain_results <- mstan4bart_fit(bartData = bartData, x = glmod$X,
                                  y = y, weights = weights, offset = offset,
                                  family = family,
                                  prior,
                                  prior_intercept,
                                  prior_aux,
                                  prior_covariance,
                                  adapt_delta = adapt_delta,
                                  group = group,
                                  bart_offset_init = fitted(lmerFit),
                                  sigma_init = sigma(lmerFit),
                                  ...)
  # TODO: de-list results into md-array
  result <- package_samples(chain_results, colnames(glmod$X), colnames(bartData@x))
  
  if (!is.null(glmod$X)) {
    result$X <- glmod$X
    result$X_means <- apply(glmod$X, 2L, mean)
  }
  if (!is.null(glmod$X.test))      result$X.test <- glmod$X.test
  if (!is.null(glmod$reTrms))      result$Zt <- glmod$reTrms$Zt
  if (!is.null(glmod$reTrms.test)) result$Zt.test <- glmod$reTrms.test$Zt
  
  # result$chain_results <- chain_results
  
  class(result) <- "mstan4bartFit"
  return(result)
  
  sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  X <- X[, !sel, drop = FALSE]
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  # colnames(Z) <- b_names(names(fit), value = TRUE)
  colnames(Z) <- colnames(fit[[1L]]$ranef)
  
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = cbind(X, Z), y = y, data, call, terms = NULL, model = NULL,
               na.action = attr(glmod$fr, "na.action"), contrasts, glmod)
  
  return(fit)
  
  #out <- stanreg(fit)
  #class(out) <- c(class(out), add_classes)
  #class(out) <- "rstanbartarm"
  
  #return(out)
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
    result$Sigma <- lapply(seq_len(n_ranef_levels), function(i_level) {
      names <- dimnames(chain_results[[1L]]$sample$stan$Sigma[[i_level]])
      names[4L] <- list(NULL)
      names(names)[3L:4L] <- c("sample", "chain")
      array(sapply(seq_len(n_chains), function(i_chain)
              chain_results[[i_chain]]$sample$stan$Sigma[[i_level]]),
            c(n_ranef_at_level[i_level], n_ranef_at_level[i_level], n_samples, n_chains),
            dimnames = names)
    })
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
                                      dim = c(n_obs, n_samples, n_chains),
                                      dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
    if (n_obs_test > 0L) {
      result$warmup$bart_test <- array(sapply(seq_along(n_chains), function(i_chains)
                                         chain_results[[i_chains]]$warmup$bart$test),
                                       dim = c(n_obs_test, n_samples, n_chains),
                                       dimnames = list(obseration = NULL, sample = NULL, chain = NULL))
    }
    result$warmup$bart_varcount <- array(sapply(seq_along(n_chains), function(i_chains)
                                           chain_results[[i_chains]]$sample$bart$varcount),
                                         dim = c(n_bart_vars, n_samples, n_chains),
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
              c(n_ranef_at_level[i_level], n_groups_at_level[i_level], n_samples, n_chains),
              dimnames = names)
      })
      result$warmup$Sigma <- lapply(seq_len(n_ranef_levels), function(i_level) {
        names <- dimnames(chain_results[[1L]]$warmup$stan$Sigma[[i_level]])
        names[4L] <- list(NULL)
        names(names)[3L:4L] <- c("sample", "chain")
        array(sapply(seq_len(n_chains), function(i_chain)
                chain_results[[i_chain]]$warmup$stan$Sigma[[i_level]]),
              c(n_ranef_at_level[i_level], n_ranef_at_level[i_level], n_samples, n_chains),
              dimnames = names)
      })
    }
    if (n_fixef > 0L) {
      names <- dimnames(chain_results[[1L]]$warmup$stan$fixef)
      names[[1L]] <- fixef_names
      names[3L] <- list(NULL)
      names(names) <- c("predictor", "sample", "chain")
      
      result$warmup$fixef <- array(sapply(seq_len(n_chains), function(i_chain)
                                     chain_results[[i_chain]]$warmup$stan$fixef),
                                   c(n_fixef, n_samples, n_chains),
                                   dimnames = names)
    }
  }
  
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

