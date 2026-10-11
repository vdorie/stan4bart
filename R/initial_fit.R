# The initial fit, made before sampling from what stan4bart() has parsed. A
# continuous fit takes from it the response range of its BART component - the
# two numbers dbarts maps onto [-0.5, 0.5] and states its leaf prior between,
# one pair shared by every chain - and the values its first sweep starts from.
# A binary fit has no range and takes the starting values alone.

# 'bart_range' as stan4bart() takes it: the name of an estimate, or a pair of
# numbers that is used as given and never reordered.
resolve_bart_range <- function(bart_range)
{
  if (is.null(bart_range)) return("mixed")
  if (is.character(bart_range))
    return(match.arg(bart_range, c("mixed", "lm", "response")))
  if (!is.numeric(bart_range) || length(bart_range) != 2L || !all(is.finite(bart_range)) ||
      !(bart_range[1L] < bart_range[2L]) || !is.finite(bart_range[2L] - bart_range[1L]))
    stop("'bart_range' must be c(lower, upper), lower below upper, with a finite width",
         call. = FALSE)
  as.double(bart_range)
}

format_range <- function(x) vapply(x, format, "", digits = 6L)

# The rows an initial fit reads, the kept rows of positive weight, with their
# weights, and what the user's offset stands for: under an 'offset_type' other
# than "default" it replaces a part of the model, and with no offset the type
# does nothing. 'target' is what the forest is fitted to before the model's
# own part comes off: the response less the offset, or the response alone
# where the offset stands for the forest.
initial_fit_rows <- function(object)
{
  weights <- object$weights
  rows_in <- if (length(weights) > 0L) which(weights > 0) else seq_along(object$y)
  has_offset <- length(object$offset) > 0L && any(object$offset != 0)
  offset_type <- if (has_offset) object$offset_type else "default"
  offset <- if (has_offset && offset_type != "bart") as.vector(object$offset) else NULL
  list(rows_in = rows_in,
       w = if (length(weights) > 0L) as.double(weights[rows_in]) else rep_len(1, length(rows_in)),
       offset_type = offset_type, offset = offset,
       target = if (is.null(offset)) object$y else object$y - offset)
}

# The model's own part of an initial fit on every kept row: the fixed columns
# centered where the parametric sampler centers them times their coefficients,
# plus the predicted random effects through the whole design, and neither the
# intercept nor the forest's columns' part. Under an 'offset_type' other than
# "default" it holds what the forest's offset holds in the sweeps: the random
# part alone under "fixef", the fixed part alone under "ranef".
initial_fit_own <- function(object, columns, fit, offset_type)
{
  fixed <- if (length(columns$own) == 0L) 0 else
    drop(sweep(object$X[, columns$own_cols, drop = FALSE], 2L,
               object$X_means[columns$own_cols]) %*% fit$beta[columns$own])
  random <- if (is.null(fit[["b"]])) 0 else
    as.vector(Matrix::crossprod(object$reTrms$Zt, fit[["b"]]))
  numeric(length(object$y)) + switch(offset_type, fixef = random, ranef = fixed, fixed + random)
}

# The fixed columns of the initial fit on the rows in: an intercept, the
# model's own columns and, with 'forest', the forest's, in that order. A
# numeric or ordered forest column enters as the value the forest splits on.
# An unordered factor enters as indicators of its levels present among the
# rows in, less the first, unless each of those levels lies within one level
# of a grouping factor of the model: its indicators would take over that
# random effect, so it is left out. A pivoted QR at lm's tolerance, of the
# block as the fit weights it, then keeps the first of any dependent set,
# which puts the model's own columns ahead of the forest's; lme4 stops on a
# dependent column and so is never left to choose. 'own' indexes the kept
# columns that are the model's and 'own_cols' the columns of the model's
# design they are.
range_fixed_block <- function(object, rows_in, w, forest = TRUE)
{
  x_bart <- as.matrix(object$bartData@x)
  # dbarts marks an unordered factor's column with 1
  is_factor <- object$bartData@varTypes == 1L
  factor_levels <- attr(object$bartData@x, "factor.levels")
  added <- list()
  left_out <- character()
  for (j in seq_len(if (forest) ncol(x_bart) else 0L)) {
    name <- colnames(x_bart)[j]
    value <- x_bart[rows_in, j]
    if (!is_factor[j]) {
      added[[length(added) + 1L]] <- matrix(as.double(value), dimnames = list(NULL, name))
      next
    }
    nested <- vapply(object$reTrms$flist, function(group)
      all(tapply(as.integer(group)[rows_in], value, function(level) length(unique(level)) == 1L)),
      NA)
    if (any(nested)) {
      left_out <- c(left_out, name)
      next
    }
    present <- sort(unique(value))[-1L]
    if (length(present) == 0L) next
    indicators <- outer(value, present, "==") + 0
    colnames(indicators) <- paste0(name, factor_levels[[j]][present + 1L])
    added[[length(added) + 1L]] <- indicators
  }

  n_own <- NCOL(object$X)
  block <- cbind("(Intercept)" = rep_len(1, length(rows_in)),
                 object$X[rows_in, , drop = FALSE], do.call(cbind, added))
  decomposition <- qr(sqrt(w) * block, tol = 1e-7)
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  is_own <- keep > 1L & keep <= 1L + n_own
  list(X = block[, keep, drop = FALSE],
       own = which(is_own),
       own_cols = keep[is_own] - 1L,
       added = colnames(block)[keep[keep > 1L + n_own]],
       dropped = colnames(block)[setdiff(seq_len(ncol(block)), keep)],
       left.out = left_out)
}

# The initial fit of z on the fixed columns X with weights w, all on the rows
# in: a linear model, or with 'reTrms', the model's random-effect structure
# cut to those rows, a linear mixed model by restricted maximum likelihood.
# lme4 is handed the parsed pieces and not a formula; its constructor of a
# fitted model refuses a frame that carries none, so the estimates are read
# from the criterion's environment, as that constructor reads them. Warnings
# are the caller's to muffle: a singular fit or an optimizer's complaint
# leaves the estimates in use.
range_fit <- function(z, w, X, reTrms = NULL)
{
  if (is.null(reTrms)) {
    fit <- lm.wfit(X, z, w)
    return(list(beta = unname(fit$coefficients),
                sigma = sqrt(sum(w * fit$residuals^2) / fit$df.residual)))
  }
  control <- lme4::lmerControl(check.conv.grad     = "ignore",
                               check.conv.singular = "ignore",
                               check.conv.hess     = "ignore")
  frame <- model.frame(z ~ 1, data.frame(z = z, w = w), weights = w)
  criterion <- lme4::mkLmerDevfun(frame, X, reTrms, REML = TRUE, control = control)
  optimum <- lme4::optimizeLmer(criterion, optimizer = control$optimizer,
                                restart_edge = control$restart_edge,
                                boundary.tol = control$boundary.tol,
                                control = control$optCtrl, calc.derivs = FALSE,
                                use.last.params = control$use.last.params)
  criterion(optimum$par)
  state <- environment(criterion)
  beta <- state$pp$beta(1)
  list(beta = as.vector(beta), b = as.vector(state$pp$b(1)), theta = optimum$par,
       sigma = sqrt((state$resp$wrss() + state$pp$sqrL(1)) / (length(z) - length(beta))),
       message = optimum$message)
}

# The model's random-effect structure with its design cut to the rows in, for
# lme4. A level left with no row keeps its place and takes a predicted effect
# of zero. lme4 writes its estimates into the covariance factor and its
# parameters where they lie, so it is given copies: the fit's own must not
# change under it, or a second fit would start from the first one's end.
range_cut_reTrms <- function(reTrms, rows_in)
{
  reTrms$Zt <- reTrms$Zt[, rows_in, drop = FALSE]
  reTrms$flist[] <- lapply(reTrms$flist, function(group) group[rows_in])
  reTrms$Ztlist <- NULL
  reTrms$theta <- reTrms$theta + 0
  reTrms$Lambdat@x <- reTrms$Lambdat@x + 0
  reTrms
}

# The range of a continuous fit and the values its first sweep starts from,
# from what stan4bart() has parsed. The fit and the range read the rows in and
# no other. The range is that of what the forest is fitted to: the response,
# less the offset, less the model's own part of the initial fit, which holds
# the forest's columns as linear terms. That own part is the forest's first
# offset, on every kept row, and the spread of what is left on the rows in is
# the first residual standard deviation. An estimated range is never wider
# than the response's own, and where the fit cannot be computed the response's
# own range stands in with the start that goes with it: an own part of zero.
# Under an 'offset_type' of "parametric" the own part is zero and no fit is
# made.
bart_range_route <- function(object, bart_range, verbose = 0L, fitter = range_fit)
{
  verbose <- isTRUE(as.integer(verbose) > 0L)
  rows <- initial_fit_rows(object)
  rows_in <- rows$rows_in
  target <- rows$target
  n <- length(rows_in)

  supplied <- is.numeric(bart_range)
  route <- if (supplied) "supplied" else bart_range
  reTrms <- object$reTrms
  mixed <- route == "mixed" && length(reTrms) > 0L

  own <- numeric(length(target))
  coefficients <- NULL
  fallback <- NULL
  info <- list(n = n, p = NA_integer_, fit = "none", added = character(),
               dropped = character(), left.out = character(), theta = NULL,
               sigma = NULL, message = NULL)

  if (route != "response" && rows$offset_type != "parametric") {
    if (mixed && !lme4_is_used()) {
      fallback <- "lme4"
      warning("lme4 is not installed, so the BART response range is the response's own, ",
              "not the mixed model's; install lme4, or set 'bart_range'", call. = FALSE)
    } else {
      columns <- range_fixed_block(object, rows_in, rows$w)
      p <- ncol(columns$X)
      info$p <- p
      info[c("added", "dropped", "left.out")] <- columns[c("added", "dropped", "left.out")]
      # lme4's own refusal, which its builder of the criterion does not make
      n_effects <- if (mixed) vapply(seq_along(reTrms$cnms), function(i)
        length(reTrms$cnms[[i]]) *
          length(unique(reTrms$flist[[attr(reTrms$flist, "assign")[i]]][rows_in])), 0) else 0
      if (n - p < 1L || any(n_effects >= n)) {
        fallback <- "undefined"
      } else {
        if (verbose && !supplied) {
          cat("estimating the BART response range by a linear ",
              if (mixed) "mixed model (lme4)" else "model", " on ", n, " rows and ", p,
              " fixed columns; see 'bart_range'\n", sep = "")
          if (length(columns$added) == 0L)
            cat("no BART predictor enters the initial fit, so the range is that of the ",
                "model's own terms; see 'bart_range'\n", sep = "")
        }
        fit <- tryCatch(suppressMessages(suppressWarnings(
          if (mixed) fitter(target[rows_in], rows$w, columns$X, range_cut_reTrms(reTrms, rows_in))
          else fitter(target[rows_in], rows$w, columns$X))), error = function(e) e)
        fit <- checked_initial_fit(fit, p, if (mixed) nrow(reTrms$Zt) else 0L)
        if (inherits(fit, "error")) {
          fallback <- "failed"
          if (!supplied)
            warning("the initial fit for the BART response range failed (",
                    conditionMessage(fit), "); the response's own range is used; ",
                    "see 'bart_range'", call. = FALSE)
        } else {
          info$fit <- if (mixed) "lme4" else "linear"
          info[c("theta", "sigma", "message")] <- fit[c("theta", "sigma", "message")]
          coefficients <- setNames(fit$beta, colnames(columns$X))
          own <- initial_fit_own(object, columns, fit, rows$offset_type)
        }
      }
    }
  }

  range_own <- range(target[rows_in])
  no_width <- FALSE
  if (is.null(fallback) && info$fit != "none") {
    width <- diff(range((target - own)[rows_in]))
    # a response with no width of its own has a fitted range of none: what
    # rounding leaves of one is not a wider range
    if (diff(range_own) == 0 || width == 0) {
      fallback <- "undefined"
      no_width <- TRUE
    } else if (width > diff(range_own)) fallback <- "wider"
  }
  # a fallback takes the response's own range whole, starting values included
  if (!is.null(fallback)) own <- numeric(length(target))

  range_bart <- if (supplied) bart_range else range((target - own)[rows_in])
  # dbarts's own treatment of a constant response stands
  if (!(range_bart[2L] > range_bart[1L])) range_bart <- NULL

  if (supplied) {
    # a fallback of the start is not one of the range
    fallback <- NULL
    fitted_to <- range((target - own)[rows_in])
    if (bart_range[2L] < fitted_to[1L] || bart_range[1L] > fitted_to[2L])
      warning("'bart_range' (", paste(format_range(bart_range), collapse = ", "),
              ") lies wholly outside the values BART is fitted to (",
              paste(format_range(fitted_to), collapse = " to "), ")", call. = FALSE)
  } else if (verbose && !is.null(range_bart) && isTRUE(fallback %in% c("undefined", "wider"))) {
    cat("the BART response range is the response's own, ",
        paste(format_range(range_bart), collapse = " to "), ": ",
        if (fallback == "wider") "the fitted range was wider"
        else if (no_width) "the fitted range has no width"
        else paste("the initial fit is not defined on", n, "rows"),
        "; see 'bart_range'\n", sep = "")
  }

  sigma <- sd((target - own)[rows_in])
  if (!is.finite(sigma) || sigma == 0) sigma <- 1

  list(range = range_bart, route = route, fallback = fallback, info = info, own = own,
       sigma = sigma, coefficients = coefficients)
}

# The initial fit of a binary response y on the fixed columns X with an
# offset, all on the rows in: a probit regression, or with 'reTrms' a probit
# mixed model by lme4, built from the parsed pieces and optimized in the two
# stages glmer runs. 'usable' is FALSE where the regression reports, as
# glm.fit does, that it did not converge or fitted a probability of 0 or 1.
binary_fit <- function(y, offset, X, reTrms = NULL)
{
  family <- binomial(link = "probit")
  if (is.null(reTrms)) {
    fit <- glm.fit(X, y, offset = offset, family = family)
    eps <- 10 * .Machine$double.eps
    return(list(beta = unname(fit$coefficients),
                usable = fit$converged &&
                  !any(fit$fitted.values > 1 - eps | fit$fitted.values < eps)))
  }
  control <- lme4::glmerControl(check.conv.grad     = "ignore",
                                check.conv.singular = "ignore",
                                check.conv.hess     = "ignore")
  frame <- model.frame(y ~ 1, data.frame(y = y), offset = offset)
  # lme4's modular sequence as its help gives it (?lme4::modular), with the
  # stage of each optimization named as glmer names it: an lme4 before 2.0
  # takes an unnamed stage as the first and returns no coefficients. The
  # criterion's enclosure is the frame its builder is called from, this one,
  # and its first stage calls lme4's exported GHrule by name, so the name is
  # bound here.
  GHrule <- lme4::GHrule
  criterion <- lme4::mkGlmerDevfun(frame, X, reTrms, family = family, control = control)
  lme4::optimizeGlmer(criterion, optimizer = control$optimizer[[1L]], restart_edge = FALSE,
                      boundary.tol = 0, control = control$optCtrl, nAGQ = 0L, stage = 1,
                      calc.derivs = FALSE)
  criterion <- lme4::updateGlmerDevfun(criterion, reTrms, nAGQ = 1L)
  optimum <- lme4::optimizeGlmer(criterion, optimizer = control$optimizer[[2L]],
                                 restart_edge = control$restart_edge,
                                 boundary.tol = control$boundary.tol,
                                 control = control$optCtrl, nAGQ = 1L, stage = 2,
                                 calc.derivs = FALSE,
                                 use.last.params = control$use.last.params)
  criterion(optimum$par)
  state <- environment(criterion)
  # the second stage optimizes over the covariance parameters and then the
  # coefficients, which the criterion carries in the offset
  list(beta = optimum$par[-state$dpars], b = as.vector(state$pp$b(1)), usable = TRUE)
}

# An initial fit's estimates, or the error that stands for a failed fit: one
# the fitter raised, estimates that are not a coefficient for every fixed
# column and a predicted effect for every column of the random design, or one
# that is not finite. A failed fit is never a start.
checked_initial_fit <- function(fit, n_fixed, n_random)
{
  if (inherits(fit, "error")) return(fit)
  if (length(fit$beta) != n_fixed || length(fit[["b"]]) != n_random)
    return(simpleError("its estimates are not of the model's size"))
  if (!all(is.finite(c(fit$beta, fit[["b"]]))))
    return(simpleError("an estimate is not finite"))
  fit
}

# The values a binary fit's first sweep starts from: the model's own part of a
# probit fit of the response on the model's own columns, on the scale of the
# linear predictor. The fit reads the rows the 0/1 weights select. It is a
# mixed model where the model has random terms and lme4 is used, and else a
# regression on the fixed columns. Where no fit is to be had the start is zero.
binary_start <- function(object, verbose = 0L, fitter = binary_fit)
{
  rows <- initial_fit_rows(object)
  rows_in <- rows$rows_in
  own <- numeric(length(object$y))
  if (rows$offset_type == "parametric") return(own)

  reTrms <- object$reTrms
  mixed <- length(reTrms) > 0L && lme4_is_used()
  if (length(reTrms) > 0L && !mixed) {
    # where the offset stands for the fixed effects the start is the random
    # part alone, which a regression on the fixed columns does not have
    if (rows$offset_type == "fixef") {
      warning("lme4 is not installed, so the binary fit starts from zero: with an offset that ",
              "stands for the fixed effects its start is the mixed model's random effects; ",
              "install lme4 to start from them", call. = FALSE)
      return(own)
    }
    warning("lme4 is not installed, so the binary fit starts from a probit regression on ",
            "the fixed effects alone; install lme4 to start from the mixed model",
            call. = FALSE)
  }
  columns <- range_fixed_block(object, rows_in, rows$w, forest = FALSE)
  offset <- rows$offset[rows_in]
  fit <- tryCatch(suppressMessages(suppressWarnings(
    if (mixed) fitter(object$y[rows_in], offset, columns$X, range_cut_reTrms(reTrms, rows_in))
    else fitter(object$y[rows_in], offset, columns$X))), error = function(e) e)
  fit <- checked_initial_fit(fit, ncol(columns$X), if (mixed) nrow(reTrms$Zt) else 0L)
  if (inherits(fit, "error")) {
    warning("the initial fit for the binary response failed (", conditionMessage(fit),
            "); the fit starts from zero", call. = FALSE)
  } else if (!isTRUE(fit$usable)) {
    if (isTRUE(as.integer(verbose) > 0L))
      cat("the initial probit regression did not converge or fits a probability of 0 or 1, ",
          "so the fit starts from zero\n", sep = "")
  } else {
    own <- initial_fit_own(object, columns, fit, rows$offset_type)
  }
  own
}

# A response whose range, as dbarts derives it, is the pair exactly: the lower
# end on the first row of positive weight and the upper on every other row.
bart_range_vector <- function(pair, weights, n)
{
  rows_in <- if (length(weights) > 0L) which(weights > 0) else seq_len(n)
  result <- rep_len(pair[2L], n)
  result[rows_in[1L]] <- pair[1L]
  result
}

# Puts a sampler that holds no offset on the pair: dbarts derives the range
# from a response that spans it and is then given back its own response with
# the range held. The sampler must report the pair afterwards, its width
# exactly, or the chains would not share one prior. Its middle is formed here
# as dbarts forms it, width / 2 + lower, which no pair of a finite width
# overflows, and compared relative to the pair's larger end: the exact middle
# has roundings behind it that sum to at most one unit in that end's last
# place. 4 leaves room for a build that orders the sum otherwise and passes
# no other pair's middle.
install_bart_range <- function(sampler, pair)
{
  y <- sampler$data@y
  # updateState = FALSE explicitly: dbarts's setters default to following
  # control@updateState, and a user's bart_args updateState = TRUE must not
  # populate $state on a sampler that has not run or is about to be restored.
  sampler$setResponse(bart_range_vector(pair, sampler$data@weights, length(y)),
                      updateScale = TRUE, updateState = FALSE)
  sampler$setResponse(y, updateScale = FALSE, updateState = FALSE)
  leaf_prior <- sampler$getLeafPrior()
  width <- pair[2L] - pair[1L]
  tolerance <- 4 * .Machine$double.eps * max(abs(pair))
  if (!identical(leaf_prior$response.scale, width) ||
      !isTRUE(abs(leaf_prior$response.shift - (width * 0.5 + pair[1L])) <= tolerance))
    stop("internal error: the BART response range was not installed")
  invisible(sampler)
}
