# the response range of the BART component: one pair, fixed before sampling
# from an initial fit that holds the forest's columns, shared by every chain

has_lme4 <- requireNamespace("lme4", quietly = TRUE)

# the value of an expression, the warnings it raised and the lines it printed
observe <- function(expr) {
  warned <- character()
  printed <- capture.output(value <- withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  }))
  list(value = value, warned = warned, printed = printed)
}
# the paths that run without lme4, on a machine that has it
without_lme4 <- function(expr) {
  old <- options(stan4bart.lme4 = FALSE)
  on.exit(options(old))
  expr
}
route_of <- function(...) stan4bart:::bart_range_route(...)
# a short fit that keeps its trees, the call evaluated where it is written so
# that 'weights' and 'offset' name columns of the data
fit_short <- function(...) {
  call <- match.call()
  call[[1L]] <- quote(stan4bart)
  defaults <- list(cores = 1L, verbose = -1L, chains = 2L, warmup = 10L, iter = 20L, seed = 5L,
                   bart_args = quote(list(n.trees = 5L, keepTrees = TRUE)))
  for (name in setdiff(names(defaults), names(call))) call[[name]] <- defaults[[name]]
  eval(call, parent.frame())
}
state_pairs <- function(fit)
  vapply(fit$state.bart$state, function(state) state[[1L]]$fit.scale, double(2L))

set.seed(21L)
n <- 200L
d <- data.frame(g = factor(sample(16L, n, TRUE)), h = factor(sample(5L, n, TRUE)),
                x = rnorm(n), s = runif(n), u = runif(n), v = runif(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)), o = rnorm(n))
d$w <- runif(n, 0.5, 2)
d$w[seq(3L, n, 7L)] <- 0
# a grouping level with no row of positive weight
d$w[d$g == "4"] <- 0
effects <- matrix(rnorm(32L), 16L) %*% chol(matrix(c(4, 1, 1, 1), 2L))
d$y <- 50 + 4 * d$x + 8 * sin(pi * d$u * d$v) + c(0, 1, -1)[d$f] + effects[d$g, 1L] +
  d$s * effects[d$g, 2L] + rnorm(5L)[d$h] + d$o + rnorm(n)
# rows of weight zero hold the extremes of the response less the offset
d$y[3L] <- d$y[3L] + 40
d$y[10L] <- d$y[10L] - 40
d$xcopy <- d$x
rows_in <- d$w > 0
own_range <- range((d$y - d$o)[rows_in])

mixed <- y ~ x + (1 + s | g) + (1 | h) + bart(u + v + f + xcopy)
plain <- y ~ x + bart(u + v + f + xcopy)
mixed0 <- stan4bart(mixed, d, weights = w, offset = o, iter = 0L)
plain0 <- stan4bart(plain, d, weights = w, offset = o, iter = 0L)

# what "lm" is defined as, computed here: the weighted fit of the response
# less the offset on the model's own column and the forest's, over the rows in
lm_block <- cbind(1, d$x, d$u, d$v, d$f == "b", d$f == "c")
lm_fit <- lm.wfit(lm_block[rows_in, ], (d$y - d$o)[rows_in], d$w[rows_in])
lm_own <- (d$x - mean(d$x)) * lm_fit$coefficients[[2L]]
lm_range <- range((d$y - d$o - lm_own)[rows_in])

# the definition ----------------------------------------------------------

# "mixed" is lmer on the written-out formula over the rows of positive weight
if (has_lme4) local({
  control <- lme4::lmerControl(check.conv.grad = "ignore", check.conv.singular = "ignore",
                               check.conv.hess = "ignore")
  oracle <- lme4::lmer(y ~ x + u + v + f + (1 + s | g) + (1 | h), d[rows_in, ], weights = w,
                       offset = o, control = control)
  effects <- lme4::ranef(oracle)
  # the level with no row of positive weight has no effect in the oracle, and
  # starts at zero
  expect_true(any(d$g == "4") && !("4" %in% rownames(effects$g)))
  effects_g <- as.matrix(effects$g)[match(as.character(d$g), rownames(effects$g)), ]
  effects_g[is.na(effects_g)] <- 0
  own <- (d$x - mean(d$x)) * lme4::fixef(oracle)[["x"]] + effects_g[, 1L] + d$s * effects_g[, 2L] +
    effects$h[as.character(d$h), 1L]
  pair <- range((d$y - d$o - own)[rows_in])
  width <- diff(pair)

  route <- route_of(mixed0, "mixed")
  expect_true(all(abs(route$range - pair) < 1e-6 * width))
  expect_true(max(abs(route$own - own)) < 1e-6 * width)
  expect_true(max(abs(route$coefficients - lme4::fixef(oracle))) < 1e-6 * width)
  expect_identical(route$info[c("n", "p", "fit", "dropped", "left.out")],
                   list(n = sum(rows_in), p = 6L, fit = "lme4", dropped = "xcopy",
                        left.out = character()))
  expect_equal(route$sigma, sd((d$y - d$o - own)[rows_in]), tolerance = 1e-6)
  # lme4 is given copies of what it writes into: the parsed structure is as
  # it was, and a second estimate is the first
  expect_identical(mixed0$reTrms$theta, c(1, 0, 1, 1))
  expect_identical(route_of(mixed0, "mixed")[c("range", "own")], route[c("range", "own")])

  # the fit records what the route gave
  fit <- fit_short(mixed, d, weights = w, offset = o)
  expect_identical(fit$range.bart, route$range)
  expect_identical(fit$range.route, "mixed")
  expect_null(fit$range.fallback)
  expect_identical(fit$range.info, route$info)
  expect_identical(fit$state.bart$data@sigma, route$sigma)
  expect_identical(suppressMessages(capture.output(print(summary(fit))))[5L],
                   paste0("BART response range: ", format(route$range[1L], digits = 6L), " to ",
                          format(route$range[2L], digits = 6L), " (mixed)"))
})

# what stan4bart() hands the chains is the route's range and starting values
local({
  route <- route_of(plain0, "lm")
  direct <- stan4bart:::stan4bart_fit(plain0, plain0$family, route$own, route$sigma, route$range,
                                      -1L, 20L, 10L, 1L, 2L, 1L, 2L, 5L, TRUE, FALSE, NULL, NULL,
                                      NULL, list(n.trees = 5L))
  fit <- stan4bart(plain, d, weights = w, offset = o, bart_range = "lm", cores = 1L,
                   verbose = -1L, chains = 2L, warmup = 10L, iter = 20L, seed = 5L,
                   bart_args = list(n.trees = 5L))
  for (i in 1:2) {
    expect_identical(unname(fit$stan[, , i]), unname(direct[[i]]$sample$stan))
    expect_identical(unname(fit$bart_train[, , i]), unname(direct[[i]]$sample$bart$train))
  }
})

# the columns -------------------------------------------------------------

local({
  set.seed(22L)
  n <- 120L
  d <- data.frame(g = factor(sample(12L, n, TRUE)), x = rnorm(n), u = runif(n), v = runif(n))
  d$y <- d$x + 3 * d$u + rnorm(12L)[d$g] + rnorm(n)
  d$lg <- d$v > 0.4
  d$od <- factor(sample(c("lo", "mid", "hi"), n, TRUE), levels = c("lo", "mid", "hi"),
                 ordered = TRUE)
  d$one <- 1
  d$su <- d$u + d$v
  d$xcopy <- d$x
  d$f <- factor(sample(c("a", "b", "c", "d"), n, TRUE))
  # finer than the grouping factor, and coarser than it
  d$fine <- factor(paste(d$g, sample(2L, n, TRUE)))
  d$coarse <- factor(as.integer(d$g) %% 3L)
  fit0 <- suppressWarnings(stan4bart(
    y ~ x + (1 | g) + bart(u + v + lg + od + log(u + 1) + one + su + xcopy + f + g + fine + coarse),
    d, iter = 0L))
  columns <- stan4bart:::range_fixed_block(fit0, seq_len(n), rep_len(1, n))

  # a constant, a sum of two columns and a copy of the model's own column are
  # dropped, the model's own kept
  expect_identical(columns$dropped, c("one", "su", "xcopy"))
  expect_identical(colnames(columns$X)[columns$own], "x")
  expect_identical(columns$own_cols, 1L)
  expect_identical(unname(columns$X[, 2L]), d$x)
  # the grouping factor and a factor finer than it are left out and named; an
  # unordered factor is indicators of its levels less the first, a coarser one
  # among them
  expect_identical(columns$left.out, c("g", "fine"))
  expect_identical(columns$added, c("u", "v", "lg", "od", "log(u + 1)", "fb", "fc", "fd",
                                    "coarse1", "coarse2"))
  expect_identical(colnames(columns$X), c("(Intercept)", "x", columns$added))
  # a logical, an ordered factor and a transformed term enter as the value the
  # forest splits on
  expect_identical(unname(columns$X[, "lg"]), as.double(d$lg))
  expect_identical(unname(columns$X[, "od"]), as.double(as.integer(d$od) - 1L))
  expect_identical(unname(columns$X[, "log(u + 1)"]), log(d$u + 1))
  for (level in c("b", "c", "d"))
    expect_identical(unname(columns$X[, paste0("f", level)]), as.double(d$f == level))
  expect_identical(unname(columns$X[, "coarse2"]), as.double(d$coarse == "2"))

  # a level present on no row in has no indicator
  rows <- which(d$f != "d")
  columns <- stan4bart:::range_fixed_block(fit0, rows, rep_len(1, length(rows)))
  expect_identical(grep("^f", columns$added, value = TRUE), c("fb", "fc"))
  expect_identical(nrow(columns$X), length(rows))
})

# a forest holding only the model's own column: the fit is the model's own,
# and a verbose line says so
local({
  alone0 <- stan4bart(y ~ x + bart(xcopy), d, weights = w, offset = o, iter = 0L)
  seen <- observe(route_of(alone0, "mixed", verbose = 1L))
  expect_identical(seen$printed, c(
    paste("estimating the BART response range by a linear model on", sum(rows_in),
          "rows and 2 fixed columns; see 'bart_range'"),
    "no BART predictor enters the initial fit, so the range is that of the model's own terms; see 'bart_range'"))
  expect_identical(seen$warned, character())
  fit <- lm.wfit(cbind(1, d$x)[rows_in, ], (d$y - d$o)[rows_in], d$w[rows_in])
  pair <- range((d$y - d$o - (d$x - mean(d$x)) * fit$coefficients[[2L]])[rows_in])
  expect_true(all(abs(seen$value$range - pair) < 1e-10 * diff(pair)))
  expect_identical(seen$value$info[c("added", "dropped")],
                   list(added = character(), dropped = "xcopy"))
  expect_identical(observe(route_of(alone0, "mixed", verbose = 0L))$printed, character())
})

# the routes --------------------------------------------------------------

local({
  # "lm" is the linear fit's range, also on a model with random terms
  for (fit0 in list(plain0, mixed0)) {
    route <- route_of(fit0, "lm")
    expect_true(all(abs(route$range - lm_range) < 1e-10 * diff(lm_range)))
    expect_true(max(abs(route$own - lm_own)) < 1e-10 * diff(lm_range))
    expect_identical(route[c("route", "fallback")], list(route = "lm", fallback = NULL))
    expect_identical(route$info$fit, "linear")
    expect_equal(route$sigma, sd((d$y - d$o - lm_own)[rows_in]), tolerance = 1e-10)
  }
  # the range read no row of weight zero, and is narrower than the response's
  expect_true(diff(lm_range) < diff(own_range) &&
              diff(own_range) + 50 < diff(range(d$y - d$o)))

  # a model with no random term asked for "mixed" is the linear fit
  route <- route_of(plain0, "mixed")
  expect_identical(route$range, route_of(plain0, "lm")$range)
  expect_identical(route[c("route", "fallback")], list(route = "mixed", fallback = NULL))
  expect_identical(route$info$fit, "linear")

  # "response" fits nothing
  route <- route_of(mixed0, "response", fitter = function(...) stop("fitted"))
  expect_identical(route$range, own_range)
  expect_identical(route[c("route", "fallback")], list(route = "response", fallback = NULL))
  expect_identical(route$own, numeric(n))
  expect_identical(route$sigma, sd((d$y - d$o)[rows_in]))
  expect_identical(route$info$fit, "none")

  # a pair is used as it is, far wider than the response's own or inside it,
  # and starts from the linear fit
  for (pair in list(c(-150.5, 300.25), c(49, 51))) {
    route <- route_of(mixed0, pair)
    expect_identical(route$range, pair)
    expect_identical(route[c("route", "fallback")], list(route = "supplied", fallback = NULL))
    expect_true(max(abs(route$own - lm_own)) < 1e-10 * diff(lm_range))
    expect_equal(route$sigma, sd((d$y - d$o - lm_own)[rows_in]), tolerance = 1e-10)
  }
})

# lme4's builder of the criterion is called once under "mixed" on a model with
# random terms, and by nothing else
if (has_lme4) local({
  calls <- 0L
  invisible(suppressMessages(trace("mkLmerDevfun", where = asNamespace("lme4"), print = FALSE,
                                   tracer = function() calls <<- calls + 1L)))
  on.exit(invisible(suppressMessages(untrace("mkLmerDevfun", where = asNamespace("lme4")))))
  route_of(mixed0, "mixed")
  expect_identical(calls, 1L)
  for (bart_range in list("lm", "response", c(0, 100))) route_of(mixed0, bart_range)
  route_of(plain0, "mixed")
  expect_identical(calls, 1L)
})

# when the range is the response's own though another was asked -----------------

# lme4 not installed: a warning, for a model with random terms under "mixed"
# and for nothing else
local({
  seen <- without_lme4(observe(route_of(mixed0, "mixed", verbose = 1L)))
  expect_identical(seen$value$range, own_range)
  expect_identical(seen$value[c("route", "fallback")], list(route = "mixed", fallback = "lme4"))
  expect_identical(seen$warned, paste0(
    "lme4 is not installed, so the BART response range is the response's own, not the mixed ",
    "model's; install lme4, or set 'bart_range'"))
  expect_identical(seen$printed, character())
  expect_identical(seen$value$own, numeric(n))
  expect_identical(seen$value$sigma, sd((d$y - d$o)[rows_in]))

  for (bart_range in list("lm", "response", c(0, 100))) {
    seen <- without_lme4(observe(route_of(mixed0, bart_range)))
    expect_identical(seen$warned, character())
    expect_null(seen$value$fallback)
  }
  seen <- without_lme4(observe(route_of(plain0, "mixed")))
  expect_identical(seen$warned, character())
  expect_null(seen$value$fallback)
  expect_true(all(abs(seen$value$range - lm_range) < 1e-10 * diff(lm_range)))

  # the fit records it, and its chains hold the response's own range
  seen <- without_lme4(observe(fit_short(mixed, d, weights = w, offset = o)))
  expect_identical(length(seen$warned), 1L)
  expect_identical(seen$value$range.bart, own_range)
  expect_identical(seen$value$range.fallback, "lme4")
  expect_identical(state_pairs(seen$value), matrix(own_range, 2L, 2L))
  expect_identical(suppressMessages(capture.output(print(summary(seen$value))))[5L],
                   paste0("BART response range: ", format(own_range[1L], digits = 6L), " to ",
                          format(own_range[2L], digits = 6L), " (mixed; the response's own: lme4)"))
})

# the fit is not defined: a verbose line and no warning
local({
  undefined_line <- function(pair, n)
    paste0("the BART response range is the response's own, ", format(pair[1L], digits = 6L),
           " to ", format(pair[2L], digits = 6L), ": the initial fit is not defined on ", n,
           " rows; see 'bart_range'")

  # 30 rows and 40 forest columns
  set.seed(23L)
  wide <- data.frame(y = rnorm(30L), x = rnorm(30L), matrix(rnorm(30L * 40L), 30L))
  wide0 <- stan4bart(y ~ x + bart(. - x), wide, iter = 0L)
  for (bart_range in c("mixed", "lm")) {
    seen <- observe(route_of(wide0, bart_range, verbose = 1L))
    expect_identical(seen$value$range, range(wide$y))
    expect_identical(seen$value[c("route", "fallback")],
                     list(route = bart_range, fallback = "undefined"))
    expect_identical(seen$value$info[c("n", "p", "fit")], list(n = 30L, p = 30L, fit = "none"))
    expect_identical(seen$warned, character())
    expect_identical(seen$printed, undefined_line(range(wide$y), 30L))
    expect_identical(observe(route_of(wide0, bart_range, verbose = 0L))$printed, character())
  }

  # one random effect per row, which lmer refuses and lme4's builder does not
  if (has_lme4) {
    d$id <- factor(seq_len(n))
    each0 <- stan4bart(y ~ x + (1 | id) + bart(u + v), d, weights = w, offset = o, iter = 0L)
    seen <- observe(route_of(each0, "mixed", verbose = 1L,
                             fitter = function(...) stop("fitted")))
    expect_identical(seen$value$range, own_range)
    expect_identical(seen$value[c("route", "fallback")],
                     list(route = "mixed", fallback = "undefined"))
    expect_identical(seen$warned, character())
    expect_identical(seen$printed, undefined_line(own_range, sum(rows_in)))
    # one effect fewer than rows in is defined
    d$id[which(rows_in)[1:2]] <- d$id[which(rows_in)[1L]]
    each0 <- stan4bart(y ~ x + (1 | id) + bart(u + v), d, weights = w, offset = o, iter = 0L)
    expect_null(route_of(each0, "mixed", fitter = function(z, w, X, reTrms)
      list(beta = numeric(ncol(X)), b = numeric(nrow(reTrms$Zt))))$fallback)
  }

  # a fitted range of no width
  flat <- data.frame(x = rep(c(-1, 1), 10L), u = runif(20L))
  flat$y <- 3 + 2 * flat$x
  flat0 <- stan4bart(y ~ x + bart(u), flat, iter = 0L)
  seen <- observe(route_of(flat0, "lm", verbose = 1L, fitter = function(z, w, X)
    list(beta = c(3, 2, 0), sigma = 0)))
  expect_identical(seen$value$range, c(1, 5))
  expect_identical(seen$value$fallback, "undefined")
  expect_identical(seen$value$own, numeric(20L))
  expect_identical(seen$warned, character())
  expect_identical(seen$printed[2L], undefined_line(c(1, 5), 20L))
})

# the fit fails: a warning
local({
  for (fitter in list(function(...) stop("no optimum"),
                      function(z, w, X, reTrms = NULL) list(beta = rep_len(NA_real_, ncol(X))))) {
    seen <- observe(route_of(plain0, "mixed", verbose = 1L, fitter = fitter))
    expect_identical(seen$value$range, own_range)
    expect_identical(seen$value[c("route", "fallback")], list(route = "mixed", fallback = "failed"))
    expect_identical(length(seen$warned), 1L)
    expect_identical(seen$value$own, numeric(n))
    expect_identical(length(seen$printed), 1L)
  }
  seen <- observe(route_of(plain0, "lm", fitter = function(...) stop("no optimum")))
  expect_identical(seen$warned, paste0("the initial fit for the BART response range failed ",
                                       "(no optimum); the response's own range is used; ",
                                       "see 'bart_range'"))
  # a pair whose start fails is still the pair, and says nothing
  seen <- observe(route_of(plain0, c(0, 100), fitter = function(...) stop("no optimum")))
  expect_identical(seen$value$range, c(0, 100))
  expect_identical(seen$value$own, numeric(n))
  expect_identical(seen$warned, character())
  expect_null(seen$value$fallback)
})

# the fitted range is wider than the response's own: a forest column 0.99995
# correlated with the model's own
local({
  set.seed(24L)
  n <- 300L
  near <- data.frame(x = rnorm(n), u = runif(n), v = runif(n))
  near$xnear <- near$x + 0.01 * rnorm(n)
  near$y <- 2 * near$x + 10 * sin(pi * near$u * near$v) + rnorm(n)
  near0 <- stan4bart(y ~ x + bart(xnear + u + v), near, iter = 0L)
  fit <- lm.wfit(cbind(1, near$x, near$xnear, near$u, near$v), near$y, rep_len(1, n))
  fitted_range <- range(near$y - (near$x - mean(near$x)) * fit$coefficients[[2L]])
  # PRECONDITION: the fit's own range is the wider of the two
  expect_true(diff(fitted_range) > diff(range(near$y)))

  for (bart_range in c("mixed", "lm")) {
    seen <- observe(route_of(near0, bart_range, verbose = 1L))
    expect_identical(seen$value$range, range(near$y))
    expect_identical(seen$value[c("route", "fallback")],
                     list(route = bart_range, fallback = "wider"))
    expect_identical(seen$value$info$fit, "linear")
    expect_identical(seen$value$own, numeric(n))
    expect_identical(seen$value$sigma, sd(near$y))
    expect_identical(seen$warned, character())
    expect_identical(seen$printed[2L], paste0(
      "the BART response range is the response's own, ", format(min(near$y), digits = 6L),
      " to ", format(max(near$y), digits = 6L), ": the fitted range was wider; see 'bart_range'"))
    expect_identical(observe(route_of(near0, bart_range, verbose = 0L))$printed, character())
  }
  # a pair is not bounded, and where its start's range is wider starts at zero
  route <- route_of(near0, c(-100, 100))
  expect_identical(route$range, c(-100, 100))
  expect_identical(route$own, numeric(n))
  expect_null(route$fallback)
})

# a singular fit is used as it comes
if (has_lme4) local({
  set.seed(25L)
  n <- 80L
  flat <- data.frame(g = factor(rep(1:20, 4L)), x = rnorm(n), u = runif(n))
  flat$y <- flat$x + rep(c(-1, 1, 1, -1), each = 20L) * 0.5 + rnorm(n, 0, 0.1)
  flat0 <- stan4bart(y ~ x + (1 | g) + bart(u), flat, iter = 0L)
  seen <- observe(route_of(flat0, "mixed", verbose = 1L))
  # PRECONDITION: the group variance is estimated at zero
  expect_identical(seen$value$info$theta, 0)
  expect_identical(seen$value[c("route", "fallback")], list(route = "mixed", fallback = NULL))
  expect_identical(seen$value$info$fit, "lme4")
  expect_identical(seen$warned, character())
  expect_identical(length(seen$printed), 1L)
  expect_true(diff(seen$value$range) < diff(range(flat$y)))
})

# the argument ------------------------------------------------------------

local({
  words <- "'bart_range' must be c(lower, upper), lower below upper, with a finite width"
  for (iter in c(0L, 20L)) {
    fit_with <- function(bart_range)
      stan4bart(plain, d, weights = w, offset = o, bart_range = bart_range, iter = iter,
                cores = 1L, verbose = -1L, chains = 1L, bart_args = list(n.trees = 5L))
    expect_error(fit_with("splines"), "should be one of", fixed = TRUE)
    for (bart_range in list(function(x) x, c(1, 2, 3), NA, c(NA, 1), c(2, 1), c(1, 1),
                            c(-1e308, 1e308), c(0, Inf), list(0, 1)))
      expect_error(fit_with(bart_range), words, fixed = TRUE)
  }
  # the names abbreviate, as those of 'offset_type' and 'store' do
  expect_identical(stan4bart:::resolve_bart_range("resp"), "response")
  expect_identical(stan4bart:::resolve_bart_range(NULL), "mixed")
  expect_identical(stan4bart:::resolve_bart_range(c(2L, 5L)), c(2, 5))
  seen <- observe(fit_short(plain, d, weights = w, offset = o, bart_range = "resp"))
  expect_identical(seen$value$range.route, "response")
  expect_identical(seen$value$range.bart, own_range)
  expect_identical(seen$warned, character())

  # 'bart' no longer abbreviates 'bart_args'
  expect_error(stan4bart(plain, d, iter = 0L, bart = list(n.trees = 5L)),
               "matches multiple formal arguments", fixed = TRUE)

  # a pair that misses the data wholly is used, with one warning
  seen <- observe(fit_short(plain, d, weights = w, offset = o, bart_range = c(1000, 1001)))
  expect_identical(length(seen$warned), 1L)
  expect_true(grepl("'bart_range' (1000, 1001) lies wholly outside the values BART is fitted to (",
                    seen$warned, fixed = TRUE))
  expect_identical(seen$value$range.bart, c(1000, 1001))
  expect_identical(seen$value$range.route, "supplied")
  expect_identical(state_pairs(seen$value), matrix(c(1000, 1001), 2L, 2L))
  # one that touches them is not warned of
  touching <- c(lm_range[2L], lm_range[2L] + 1)
  expect_identical(observe(route_of(plain0, touching))$warned, character())

  # a binary response has no range
  d$yb <- as.double(d$y > median(d$y))
  seen <- observe(stan4bart(yb ~ x + s + bart(u + v), d, bart_range = "lm", cores = 1L,
                            verbose = -1L, chains = 2L, warmup = 10L, iter = 20L,
                            bart_args = list(n.trees = 5L, keepTrees = TRUE)))
  expect_identical(seen$warned, "'bart_range' is not used for a binary response")
  expect_null(seen$value$range.bart)
  expect_null(seen$value$range.route)
  expect_false(any(grepl("range", suppressMessages(capture.output(print(summary(seen$value)))))))
  seen <- observe(stan4bart(yb ~ x + s + bart(u + v), d, bart_range = "mixed", cores = 1L,
                            verbose = -1L, chains = 1L, warmup = 10L, iter = 20L,
                            bart_args = list(n.trees = 5L)))
  expect_identical(seen$warned, character())
})

# offset_type -------------------------------------------------------------

# the model's own part holds what the forest's offset holds in the sweeps
if (has_lme4) local({
  control <- lme4::lmerControl(check.conv.grad = "ignore", check.conv.singular = "ignore",
                               check.conv.hess = "ignore")
  parts_of <- function(oracle) {
    effects <- lme4::ranef(oracle)$g
    effects <- effects[match(as.character(d$g), rownames(effects)), 1L]
    effects[is.na(effects)] <- 0
    list(fixed = (d$x - mean(d$x)) * lme4::fixef(oracle)[["x"]], random = effects)
  }
  # the offset is taken off the response except where it stands for the forest
  less <- parts_of(lme4::lmer(y ~ x + u + v + (1 | g), d[rows_in, ], weights = w, offset = o,
                              control = control))
  whole <- parts_of(lme4::lmer(y ~ x + u + v + (1 | g), d[rows_in, ], weights = w,
                               control = control))
  expected <- list(default = d$y - d$o - less$fixed - less$random,
                   fixef = d$y - d$o - less$random,
                   ranef = d$y - d$o - less$fixed,
                   bart = d$y - whole$fixed - whole$random)
  for (offset_type in names(expected)) {
    fit0 <- stan4bart(y ~ x + (1 | g) + bart(u + v), d, weights = w, offset = o,
                      offset_type = offset_type, iter = 0L)
    route <- route_of(fit0, "mixed")
    pair <- range(expected[[offset_type]][rows_in])
    expect_true(all(abs(route$range - pair) < 1e-6 * diff(pair)), info = offset_type)
    expect_equal(route$sigma, sd(expected[[offset_type]][rows_in]), tolerance = 1e-6,
                 info = offset_type)
    expect_null(route$fallback, info = offset_type)
  }
  # PRECONDITION: the four ranges differ
  expect_identical(anyDuplicated(vapply(expected, function(v) diff(range(v[rows_in])), 0)), 0L)
})
local({
  # "parametric": the offset stands for the whole parametric part, and nothing
  # is fitted
  fit0 <- stan4bart(mixed, d, weights = w, offset = o, offset_type = "parametric", iter = 0L)
  seen <- observe(route_of(fit0, "mixed", verbose = 1L, fitter = function(...) stop("fitted")))
  expect_identical(seen$value$range, own_range)
  expect_identical(seen$value[c("route", "fallback")], list(route = "mixed", fallback = NULL))
  expect_identical(seen$value$info$fit, "none")
  expect_identical(c(seen$warned, seen$printed), character())
  # with no offset the type does nothing
  fit0 <- stan4bart(plain, d, weights = w, offset_type = "parametric", iter = 0L)
  expect_identical(route_of(fit0, "lm")$info$fit, "linear")
})

# installing the range ----------------------------------------------------

# the made-up response holds the lower end once, on the first row of positive
# weight, and a sampler given it lands on the pair exactly
local({
  n <- 40L
  pair <- c(-3.25, 17.5)
  set.seed(26L)
  x <- matrix(runif(n * 2L), n)
  y <- 10 * x[, 1L] + rnorm(n)
  weights <- list(even = rep_len(c(1, 0), n), odd = rep_len(c(0, 1), n),
                  two = replace(numeric(n), c(17L, 31L), 2), absent = NULL)
  for (name in names(weights)) {
    w <- weights[[name]]
    first <- if (is.null(w)) 1L else which(w > 0)[1L]
    made <- stan4bart:::bart_range_vector(pair, w, n)
    expect_identical(which(made == pair[1L]), first, info = name)
    expect_identical(made[-first], rep_len(pair[2L], n - 1L), info = name)

    data <- dbarts::dbartsData(y ~ x)
    data@sigma <- 1
    if (!is.null(w)) data@weights <- w
    control <- dbarts::dbartsControl(n.chains = 1L, n.samples = 1L, n.burn = 0L, n.thin = 1L,
                                     n.threads = 1L, updateState = FALSE)
    spec <- dbarts::dbartsSpec(data = data, control = control, family = gaussian(sigma = fixed(1)))
    sampler <- new("dbartsSampler", spec$control, spec$model, spec$data)
    stan4bart:::install_bart_range(sampler, pair)
    leaf_prior <- sampler$getLeafPrior()
    expect_identical(leaf_prior$response.scale, pair[2L] - pair[1L], info = name)
    expect_identical(leaf_prior$response.shift, mean(pair), info = name)
    expect_identical(sampler$data@y, y, info = name)
  }
})

# every chain is put on the one range and holds it through warmup
local({
  held <- function(fit, pair, chains = 2L) {
    expect_identical(fit$range.bart, pair)
    expect_identical(state_pairs(fit), matrix(pair, 2L, chains))
    for (sampler in stan4bart:::getBartSampler(fit)) {
      leaf_prior <- sampler$getLeafPrior()
      expect_identical(leaf_prior$response.scale, pair[2L] - pair[1L])
      expect_true(abs(leaf_prior$response.shift - mean(pair)) <= 1e-12)
    }
  }
  # a pair far from the data's own, on weights with zeros and an offset
  pair <- c(-150.5, 300.25)
  held(fit_short(mixed, d, weights = w, offset = o, bart_range = pair), pair)
  # and an estimated one
  held(fit_short(mixed, d, weights = w, offset = o, bart_range = "lm"), route_of(mixed0, "lm")$range)
  if (has_lme4)
    held(fit_short(mixed, d, weights = w, offset = o), route_of(mixed0, "mixed")$range)
  if (at_home()) {
    held(fit_short(mixed, d, weights = w, offset = o, bart_range = pair, chains = 3L, cores = 2L),
         pair, 3L)
    if (has_lme4)
      held(fit_short(mixed, d, weights = w, offset = o, chains = 3L, cores = 2L),
           route_of(mixed0, "mixed")$range, 3L)
  }
})
