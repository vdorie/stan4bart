# The parametric part's arguments (prior, prior_aux, prior_covariance, QR,
# adapt_delta, init_r) and the deprecated stan_args that routes to them.

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(60, TRUE, FALSE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y))
rm(testData)

fit_call <- function(...) {
  withCallingHandlers(
    stan4bart(y ~ bart(. - g.1 - g.2 - X4 - X5) + X4 + X5 + (1 + X4 | g.1) + (1 | g.2), df,
              cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6, seed = 7L,
              bart_args = list(n.trees = 3), ...),
    warning = function(w)
      if (startsWith(conditionMessage(w), "lme4 is not installed")) invokeRestart("muffleWarning"))
}
# the fit and the warnings it raised
catch <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}
reset_session <- function() assign("stan_args", NULL, envir = stan4bart:::.message_env)

# the new formals follow the old ones, so no positional call changes
local({
  nms <- names(formals(stan4bart))
  expect_identical(tail(nms, 7L),
                   c("bart_range", "prior", "prior_aux", "prior_covariance", "QR",
                     "adapt_delta", "init_r"))
  expect_identical(nms[match("bart_range", nms) - 1L], "bart_args")
  expect_true("stan_args" %in% nms)
})

# a prior constructor resolves inside its argument with nothing attached; a
# caller's own value of another constructor's name is untouched
local({
  expect_false("normal" %in% getNamespaceExports("stan4bart"))
  expect_inherits(fit_call(prior = normal(0, 1)), "stan4bartFit")
  expect_inherits(fit_call(prior = student_t(3, 0, 1), prior_aux = cauchy(0, 2),
                           prior_covariance = decov(2, 1, 1, 1)), "stan4bartFit")
  # a caller's value named for a constructor that is not called is its own
  wrapper <- function(p, ...)
    withCallingHandlers(
      stan4bart(y ~ bart(. - g.1 - g.2 - X4 - X5) + X4 + X5 + (1 | g.2), df,
                cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                bart_args = list(n.trees = 3), prior = p, ...),
      warning = function(w)
        if (startsWith(conditionMessage(w), "lme4 is not installed")) invokeRestart("muffleWarning"))
  # a prior held in a variable, as a wrapper passes it on
  expect_inherits(wrapper(stan4bart:::normal(0, 1)), "stan4bartFit")
  own_value <- function() {
    hs <- 5
    withCallingHandlers(
      stan4bart(y ~ bart(. - g.1 - g.2 - X4 - X5) + X4 + X5 + (1 | g.2), df,
                cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                bart_args = list(n.trees = 3), prior = normal(0, scale = hs)),
      warning = function(w)
        if (startsWith(conditionMessage(w), "lme4 is not installed")) invokeRestart("muffleWarning"))
  }
  expect_inherits(own_value(), "stan4bartFit")
  # an argument resolves its own set: exponential is no coefficient prior
  expect_error(fit_call(prior = exponential()), "prior distribution should be one of")
})

# a wrapper forwarding its dots: the prior arguments are read where the
# wrapper's caller wrote them, with the constructors bound, so a local of that
# caller is found, and the plain arguments are values, so are its locals, one
# named for a function included; the fit is the one the literals fit
local({
  wrap1 <- function(...) fit_call(...)
  wrap2 <- function(...) wrap1(...)
  literal <- fit_call(prior = normal(0, 1.5), prior_aux = exponential(rate = 2),
                      prior_covariance = decov(2, 1, 1, 1), QR = TRUE,
                      adapt_delta = 0.9, init_r = 1)
  expect_false(identical(literal$stan, fit_call()$stan))
  with_locals <- function(call_it) {
    s <- 1.5; rate <- 2; reg <- 2; ad <- 0.9; ir <- 1
    qr <- TRUE
    call_it(prior = normal(0, s), prior_aux = exponential(rate = rate),
            prior_covariance = decov(reg, 1, 1, 1), QR = qr,
            adapt_delta = ad, init_r = ir)
  }
  for (layers in list(fit_call, wrap1, wrap2))
    expect_identical(with_locals(layers)$stan, literal$stan)
  # a prior held in a local of the caller, by name
  literal_prior <- fit_call(prior = normal(0, 1.5))
  with_named <- function(call_it) {
    p <- stan4bart:::normal(0, 1.5)
    call_it(prior = p)
  }
  for (layers in list(fit_call, wrap2))
    expect_identical(with_named(layers)$stan, literal_prior$stan)
  # stan_args forwarded the same way
  with_stan_args <- function(call_it) {
    s <- 1.5
    call_it(stan_args = list(prior = normal(0, s)))
  }
  for (layers in list(fit_call, wrap2)) {
    reset_session()
    seen <- catch(with_stan_args(layers))
    expect_identical(seen$value$stan, literal_prior$stan)
  }
})

# the defaults are the ones stan_args took: naming them is the same fit
local({
  fit_unset <- fit_call()
  fit_named <- fit_call(prior = normal(0, 2.5, autoscale = TRUE),
                        prior_aux = exponential(autoscale = TRUE),
                        prior_covariance = decov(1, 1, 1, 1), QR = FALSE,
                        adapt_delta = 0.8, init_r = 2)
  expect_identical(fit_named$stan, fit_unset$stan)
  expect_identical(fit_named$bart_train, fit_unset$bart_train)
})

# QR, as rstanarm has it, and the flat prior a NULL names
local({
  expect_inherits(fit_call(QR = TRUE), "stan4bartFit")
  expect_inherits(fit_call(prior = NULL), "stan4bartFit")
  expect_inherits(fit_call(prior_aux = NULL), "stan4bartFit")
  # the covariance prior has no flat spelling
  expect_error(fit_call(prior_covariance = NULL), "'prior_covariance' can't be NULL")
})

# each stan_args entry routes to its argument, and the fit is the same one
local({
  routes <- list(
    prior = list(prior = stan4bart:::student_t(4, 0, 1)),
    prior_aux = list(prior_aux = stan4bart:::normal(0, 3)),
    prior_covariance = list(prior_covariance = stan4bart:::decov(2, 1, 1, 1)),
    QR = list(QR = TRUE),
    adapt_delta = list(adapt_delta = 0.9),
    init_r = list(init_r = 1))
  for (name in names(routes)) {
    reset_session()
    seen <- catch(do.call(fit_call, list(stan_args = routes[[name]])))
    old <- seen$value
    expect_identical(length(seen$warnings), 1L, info = name)
    expect_true(grepl(paste0("'stan_args' is deprecated.*'", name, "'"), seen$warnings), info = name)
    new <- do.call(fit_call, routes[[name]])
    expect_identical(old$stan, new$stan, info = name)
    expect_identical(old$bart_train, new$bart_train, info = name)
  }

  # one warning a session, naming every home the call used
  reset_session()
  expect_warning(
    fit_call(stan_args = list(prior = stan4bart:::normal(0, 1), QR = TRUE)),
    "'prior', 'QR'")
  expect_silent(fit_call(stan_args = list(QR = TRUE)))

  # a legacy constructor needs nothing attached, as 0.0-13 did not manage
  reset_session()
  expect_warning(fit_call(stan_args = list(prior = normal(0, 1))), "deprecated")

  # NULL stan_args and an empty list are silent
  expect_silent(fit_call(stan_args = NULL))
  expect_silent(fit_call(stan_args = list()))
})

# a thing named in both places is an error
local({
  expect_error(fit_call(QR = TRUE, stan_args = list(QR = TRUE)),
               "'QR' is given both as an argument of stan4bart\\(\\) and in stan_args")
  expect_error(fit_call(adapt_delta = 0.9, stan_args = list(adapt_delta = 0.9, init_r = 1)),
               "'adapt_delta' is given both")
  expect_error(fit_call(prior = normal(), stan_args = list(prior = stan4bart:::normal())),
               "'prior' is given both")
})

# unknown names are refused, naming the accepted ones; the two that never
# shipped are not accepted
local({
  expect_error(fit_call(stan_args = list(adapt_deta = 0.9)),
               "unrecognized name: 'adapt_deta'.*'prior', 'prior_aux'.*'adapt_delta'")
  expect_error(fit_call(stan_args = list(ridge_move = FALSE)), "unrecognized name: 'ridge_move'")
  expect_error(fit_call(stan_args = list(save_raw_parameters = TRUE)),
               "unrecognized name: 'save_raw_parameters'")
  expect_error(fit_call(stan_args = list(a = 1, b = 2)), "unrecognized names: 'a', 'b'")
  # the message says of each group what is true of it: six are arguments now,
  # ten are ignored
  msg <- tryCatch(fit_call(stan_args = list(adapt_deta = 0.9)), error = conditionMessage)
  expect_identical(msg, paste0(
    "stan_args has unrecognized name: 'adapt_deta'; it accepts 'prior', 'prior_aux', ",
    "'prior_covariance', 'QR', 'adapt_delta', 'init_r', which are now arguments of stan4bart(), and ",
    "'prior_intercept', 'adapt_gamma', 'adapt_kappa', 'adapt_t0', 'adapt_init_buffer', ",
    "'adapt_term_buffer', 'adapt_window', 'stepsize', 'stepsize_jitter', 'max_treedepth', ",
    "which are ignored"))
  expect_error(fit_call(stan_args = 3), "must be a named list")
  expect_error(fit_call(stan_args = list(1)), "must be a named list")
})

# the dead names keep their warning; prior_intercept keeps its own
local({
  expect_warning(fit_call(stan_args = list(adapt_kappa = 0.75)),
                 "'adapt_kappa' is ignored by the gradient-based sampler")
  expect_warning(fit_call(stan_args = list(prior_intercept = stan4bart:::normal())),
                 "redundantly parameterized")
})

# the internal switches are unexported options, not arguments
local({
  expect_false(any(c("ridge_move", "save_raw_parameters") %in% names(formals(stan4bart))))
  op <- options(stan4bart.save_raw_parameters = TRUE)
  on.exit(options(op))
  fit_raw <- fit_call()
  expect_true(any(startsWith(dimnames(fit_raw$stan)$parameters, "z_beta.")))
})
