# The parametric part's arguments: the priors and QR, as rstanarm spells them,
# and the sampler's adapt_delta and init_r. Each is read off the call, so a
# prior constructor written in one resolves with nothing attached.

# The constructors a prior argument may call. Only the names a call holds in
# call position are bound, in an environment whose parent is the caller's
# frame, so a caller's own value called 'normal' is untouched.
parametric_constructors <- function()
  list(normal = normal, student_t = student_t, cauchy = cauchy, hs = hs,
       hs_plus = hs_plus, laplace = laplace, lasso = lasso,
       product_normal = product_normal, exponential = exponential,
       decov = decov)

eval_parametric_arg <- function(expr, env) {
  constructors <- parametric_constructors()
  bound <- new.env(parent = env)
  for (name in intersect(names(constructors), called_names(expr)))
    assign(name, constructors[[name]], envir = bound)
  eval(expr, envir = bound)
}

# 0.0-13's stan_args entries and where each lives now
stan_args_homes <- c(prior = "prior", prior_aux = "prior_aux",
                     prior_covariance = "prior_covariance", QR = "QR",
                     adapt_delta = "adapt_delta", init_r = "init_r")
# the entries that parse and do nothing
stan_args_ignored <- c("prior_intercept", "adapt_gamma", "adapt_kappa",
                       "adapt_t0", "adapt_init_buffer", "adapt_term_buffer",
                       "adapt_window", "stepsize", "stepsize_jitter",
                       "max_treedepth")

# Reads the parametric arguments from 'exprs', the expressions of those the
# call supplied (substitute() of each, so an argument a wrapper forwarded in
# its dots is still the expression written), the deprecated stan_args among
# them, evaluated in 'env', the caller's frame. Returns what the fit
# reads: prior, prior_aux, prior_covariance, QR, init_r, adapt_delta, and the
# two internal switches. An argument the call does not name takes the default
# rstanarm has for it; NULL for prior or prior_aux is rstanarm's flat prior.
resolve_parametric_args <- function(exprs, env) {
  given <- function(name) name %in% names(exprs)
  read <- function(name) eval_parametric_arg(exprs[[name]], env)

  # a NULL prior or prior_aux stays an entry: it is the flat prior, not a
  # missing name
  result <- list(
    prior = if (given("prior")) read("prior") else default_prior_coef_gaussian(),
    prior_aux = if (given("prior_aux")) read("prior_aux") else exponential(autoscale = TRUE),
    prior_covariance = if (given("prior_covariance")) read("prior_covariance") else decov(),
    QR = if (given("QR")) read("QR") else FALSE,
    init_r = if (given("init_r")) read("init_r") else 2.0,
    adapt_delta = if (given("adapt_delta")) read("adapt_delta") else NULL)
  explicit <- vapply(names(stan_args_homes), given, logical(1L))

  if (given("stan_args")) {
    legacy <- read("stan_args")
    if (!is.null(legacy)) {
      if (!is.list(legacy) || (length(legacy) > 0L && (is.null(names(legacy)) || any(!nzchar(names(legacy))))))
        stop("'stan_args' must be a named list; its entries are now arguments of stan4bart()",
             call. = FALSE)
      names_legacy <- names(legacy)
      if (anyDuplicated(names_legacy))
        stop("stan_args names '", names_legacy[anyDuplicated(names_legacy)], "' twice", call. = FALSE)
      accepted <- c(names(stan_args_homes), stan_args_ignored)
      unknown <- setdiff(names_legacy, accepted)
      if (length(unknown) > 0L)
        stop("stan_args has unrecognized name", if (length(unknown) > 1L) "s" else "", ": ",
             paste0("'", unknown, "'", collapse = ", "), "; it accepts ",
             paste0("'", accepted, "'", collapse = ", "),
             ", which are now arguments of stan4bart()", call. = FALSE)
      both <- intersect(names_legacy, names(stan_args_homes)[explicit])
      if (length(both) > 0L)
        stop(paste0("'", both, "' is given both as an argument of stan4bart() and in stan_args; ",
                    "give it as the argument", collapse = "\n"), call. = FALSE)
      routed <- intersect(names_legacy, names(stan_args_homes))
      for (name in routed) {
        # a NULL entry is the old spelling of "not given"
        if (!is.null(legacy[[name]])) result[name] <- list(legacy[[name]])
      }
      if (length(routed) > 0L && !isTRUE(.message_env$stan_args)) {
        .message_env$stan_args <- TRUE
        warning("'stan_args' is deprecated and will be removed in a future release; its entries are now arguments of stan4bart(): ",
                paste0("'", routed, "'", collapse = ", "),
                " (once per session)", call. = FALSE)
      }
      warn_ignored_nuts_args(legacy)
      if (!is.null(legacy[["prior_intercept"]]))
        warning("intercepts for BART models are redundantly parameterized and not included")
    }
  }
  result
}
