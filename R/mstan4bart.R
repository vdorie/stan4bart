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
           treatment = NULL,
           ...,
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
  mc$data <- data
  mc$verbose <- FALSE
  mc$prior_covariance <- mc$prior_aux <-
    mc$adapt_delta <- mc$treatment <- mc$... <- NULL
  lmerFit <- suppressWarnings(eval(mc, parent.frame()))
  fit <- mstan4bart_fit(bartData = bartData, y = y, weights = weights,
                        offset = offset, family = family,
                        prior_aux = prior_aux, 
                        adapt_delta = adapt_delta,
                        group = group,
                        ranef_init = fitted(lmerFit, ranef.only = TRUE),
                        sigma_init = sigma(lmerFit),
                        ...)
  attr(fit, "Zt.obs") <- glmod$reTrms$Zt
  if (!is.null(glmod$reTrms.cf))
    attr(fit, "Zt.cf") <- glmod$reTrms.cf$Zt
  return(fit)
  
  sel <- apply(X, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  X <- X[ , !sel, drop = FALSE]
  Z <- pad_reTrms(Ztlist = group$Ztlist, cnms = group$cnms, 
                  flist = group$flist)$Z
  # colnames(Z) <- b_names(names(fit), value = TRUE)
  colnames(Z) <- colnames(fit[[1L]]$ranef)
  
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = cbind(X, Z), y = y, data, call, terms = NULL, model = NULL,
               na.action = attr(glmod$fr, "na.action"), contrasts, glmod)
  
  # TODO: get this to work with stanreg
  return(fit)
  
  #out <- stanreg(fit)
  #class(out) <- c(class(out), add_classes)
  #class(out) <- "rstanbartarm"
  
  #return(out)
}

