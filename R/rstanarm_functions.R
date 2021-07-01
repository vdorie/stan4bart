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
  if (!inherits(A, "Matrix")) 
    A <- Matrix::Matrix(A, sparse = TRUE)
  A <- Matrix::t(A)
  A <- as(A, "dgCMatrix")
  
  return(list(w = A@x, v = A@i + 1L, u = A@p + 1L))
  # return(.Call(rstan:::extract_sparse_components, A))
}



center_x <- function(x, sparse) {
  x <- as.matrix(x)
  has_intercept <- if (ncol(x) == 0) 
    FALSE else grepl("(Intercept", colnames(x)[1L], fixed = TRUE)
  
  xtemp <- if (has_intercept) x[, -1L, drop=FALSE] else x
  
  # stan4bart modification: always has an implicit intercept thanks to having
  # a bart component
  #if (has_intercept && !sparse) {
  if (!sparse) {
    xbar <- colMeans(xtemp)
    xtemp <- sweep(xtemp, 2, xbar, FUN = "-")
  }
  else xbar <- rep(0, ncol(xtemp))
  
  sel <- apply(xtemp, 2L, function(x) !all(x == 1) && length(unique(x)) < 2)
  if (any(sel)) {
    # drop any column of x with < 2 unique values (empty interaction levels)
    # exception is column of 1s isn't dropped 
    warning("Dropped empty interaction levels: ",
            paste(colnames(xtemp)[sel], collapse = ", "))
    xtemp <- xtemp[, !sel, drop = FALSE]
    xbar <- xbar[!sel]
  }
  
  return(nlist(xtemp, xbar, has_intercept))
}

handle_glm_prior <- function(prior, nvars, default_scale, link,
                             ok_dists = nlist("normal", student_t = "t", 
                                              "cauchy", "hs", "hs_plus", 
                                              "laplace", "lasso", "product_normal")) {
  if (!length(prior))
    return(list(prior_dist = 0L, prior_mean = as.array(rep(0, nvars)),
                prior_scale = as.array(rep(1, nvars)),
                prior_df = as.array(rep(1, nvars)), prior_dist_name = NA,
                global_prior_scale = 0, global_prior_df = 0,
                slab_df = 0, slab_scale = 0,
                prior_autoscale = FALSE))

  if (!is.list(prior)) 
    stop(sQuote(deparse(substitute(prior))), " should be a named list")
  
  prior_dist_name <- prior$dist
  prior_scale <- prior$scale
  prior_mean <- prior$location
  prior_df <- prior$df
  prior_mean[is.na(prior_mean)] <- 0.0
  prior_df[is.na(prior_df)] <- 1
  global_prior_scale <- 0
  global_prior_df <- 0
  slab_df <- 0
  slab_scale <- 0
  if (!prior_dist_name %in% unlist(ok_dists)) {
    stop("The prior distribution should be one of ",
         paste(names(ok_dists), collapse = ", "))
  } else if (prior_dist_name %in% 
             c("normal", "t", "cauchy", "laplace", "lasso", "product_normal")) {
    if (prior_dist_name == "normal") prior_dist <- 1L
    else if (prior_dist_name == "t") prior_dist <- 2L
    else if (prior_dist_name == "laplace") prior_dist <- 5L
    else if (prior_dist_name == "lasso") prior_dist <- 6L
    else if (prior_dist_name == "product_normal") prior_dist <- 7L
    prior_scale <- set_prior_scale(prior_scale, default = default_scale, 
                                   link = link)
  } else if (prior_dist_name %in% c("hs", "hs_plus")) {
    prior_dist <- ifelse(prior_dist_name == "hs", 3L, 4L)
    global_prior_scale <- prior$global_scale
    global_prior_df <- prior$global_df
    slab_df <- prior$slab_df
    slab_scale <- prior$slab_scale
  } else if (prior_dist_name %in% "exponential") {
    prior_dist <- 3L # only used for scale parameters so 3 not a conflict with 3 for hs
  }
  
  prior_df <- maybe_broadcast(prior_df, nvars)
  prior_df <- as.array(pmin(.Machine$double.xmax, prior_df))
  prior_mean <- maybe_broadcast(prior_mean, nvars)
  prior_mean <- as.array(prior_mean)
  prior_scale <- maybe_broadcast(prior_scale, nvars)

  nlist(prior_dist, 
        prior_mean, 
        prior_scale, 
        prior_df, 
        prior_dist_name, 
        global_prior_scale,
        global_prior_df,
        slab_df,
        slab_scale,
        prior_autoscale = isTRUE(prior$autoscale))
}

set_prior_scale <- function(scale, default, link) {
  stopifnot(is.numeric(default), is.character(link) || is.null(link))
  if (is.null(scale)) 
    scale <- default
  if (isTRUE(link == "probit"))
    scale <- scale * dnorm(0) / dlogis(0)
  
  return(scale)
}

drop_redundant_dims <- function(data) {
  drop_dim <- sapply(data, function(v) is.matrix(v) && NCOL(v) == 1)
  data[, drop_dim] <- lapply(data[, drop_dim, drop=FALSE], drop)
  return(data)
}

`%ORifNULL%` <- function(a, b) {
  if (is.null(a)) b else a
}
`%ORifINF%` <- function(a, b) {
  if (a == Inf) b else a
}

has_outcome_variable <- function(f) {
  return(f[[1]] == "~" && length(f) > 2L)
  
  tt <- terms(as.formula(f))
  if (attr(tt, "response") == 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

validate_weights <- function(w) {
  if (missing(w) || is.null(w)) {
    w <- double(0)
  } else {
    if (!is.numeric(w)) 
      stop("'weights' must be a numeric vector.", 
           call. = FALSE)
    if (any(w < 0)) 
      stop("Negative weights are not allowed.", 
           call. = FALSE)
  }
  
  return(w)
}

nlist <- function(...) {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) FALSE else nzchar(names(out))
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  } else {
    names(out)[!has_name] <- nms[!has_name]
  } 
  
  return(out)
}


maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

check_reTrms <- function(reTrms) {
  stopifnot(is.list(reTrms))
  nms <- names(reTrms$cnms)
  dupes <- duplicated(nms)
  for (i in which(dupes)) {
    original <- reTrms$cnms[[nms[i]]]
    dupe <- reTrms$cnms[[i]]
    overlap <- dupe %in% original
    if (any(overlap))
      stop("stan4bart does not permit formulas with duplicate group-specific terms.\n", 
           "In this case ", nms[i], " is used as a grouping factor multiple times and\n",
           dupe[overlap], " is included multiple times.\n", 
           "Consider using || or -1 in your formulas to prevent this from happening.")
  }
  return(invisible(NULL))
}


set_sampling_args <- function(object, prior, user_dots = list(), 
                              user_adapt_delta = NULL, ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  defaults <- default_stan_control(prior = prior, 
                                   adapt_delta = user_adapt_delta)
  
  if (!"control" %in% unms) { 
    # no user-specified 'control' argument
    args$control <- defaults
  } else { 
    # user specifies a 'control' argument
    if (!is.null(user_adapt_delta)) { 
      # if user specified adapt_delta argument to stan_* then 
      # set control$adapt_delta to user-specified value
      args$control$adapt_delta <- user_adapt_delta
    } else {
      # use default adapt_delta for the user's chosen prior
      args$control$adapt_delta <- defaults$adapt_delta
    }
    if (is.null(args$control$max_treedepth)) {
      # if user's 'control' has no max_treedepth set it to rstanarm default
      args$control$max_treedepth <- defaults$max_treedepth
    }
  }
  args$save_warmup <- FALSE
  
  return(args)
}

default_stan_control <- function(prior, adapt_delta = NULL, 
                                 max_treedepth = 15L) {
  if (!length(prior)) {
    if (is.null(adapt_delta)) adapt_delta <- 0.95
  } else if (is.null(adapt_delta)) {
    adapt_delta <- switch(prior$dist, 
                          "R2" = 0.99,
                          "hs" = 0.99,
                          "hs_plus" = 0.99,
                          "lasso" = 0.99,
                          "product_normal" = 0.99,
                          0.95) # default
  }
  nlist(adapt_delta, max_treedepth)
}

is.gaussian <- function(x) x == "gaussian"
is.binomial <- function(x) x == "binomial"
is.poisson  <- function(x) x == "poisson"
is.gamma    <- function(x) x == "Gamma"
is.ig       <- function(x) x == "inverse.gaussian"
is.nb       <- function(x) x == "neg_binomial_2"


get_m_stub <- function(m, stub = "Long") 
{
  if (is.null(m))           return(NULL)
  else if (is.numeric(m))   return(paste0(stub, m, "|"))
  else if (is.character(m)) return(paste0(m, "|"))
}

fac2bin <- function(y) 
{
  if (!is.factor(y)) 
    stop("Bug found: non-factor as input to fac2bin.", call. = FALSE)
  if (!identical(nlevels(y), 2L)) 
    stop("Bug found: factor with nlevels != 2 as input to fac2bin.", 
         call. = FALSE)
  as.integer(y != levels(y)[1L])
}

