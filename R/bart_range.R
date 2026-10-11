# The response range of the BART component: the two numbers dbarts maps onto
# [-0.5, 0.5] and states its leaf prior between. A continuous fit has one,
# fixed before sampling and shared by every chain; a binary fit has none.

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

# The fixed columns of the initial fit on the rows in: an intercept, the
# model's own columns and the forest's, in that order. A numeric or ordered
# forest column enters as the value the forest splits on. An unordered factor
# enters as indicators of its levels present among the rows in, less the
# first, unless each of those levels lies within one level of a grouping
# factor of the model: its indicators would take over that random effect, so
# it is left out. A pivoted QR at lm's tolerance, of the block as the fit
# weights it, then keeps the first of any dependent set, which puts the
# model's own columns ahead of the forest's; lme4 stops on a dependent column
# and so is never left to choose. 'own' indexes the kept columns that are the
# model's and 'own_cols' the columns of the model's design they are.
range_fixed_block <- function(object, rows_in, w)
{
  x_bart <- as.matrix(object$bartData@x)
  # dbarts marks an unordered factor's column with 1
  is_factor <- object$bartData@varTypes == 1L
  factor_levels <- attr(object$bartData@x, "factor.levels")
  added <- list()
  left_out <- character()
  for (j in seq_len(ncol(x_bart))) {
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
  reTrms$Ztlist <- NULL
  reTrms$theta <- reTrms$theta + 0
  reTrms$Lambdat@x <- reTrms$Lambdat@x + 0
  reTrms
}

# The range of a continuous fit and the values its first sweep starts from,
# from what stan4bart() has parsed. The rows in are the kept rows of positive
# weight; the fit and the range read no other row. The range is that of what
# the forest is fitted to: the response, less the offset, less the model's own
# part of the initial fit, which holds the forest's columns as linear terms.
# That own part - the fixed columns centered where the parametric sampler
# centers them times their coefficients, plus the predicted random effects,
# and neither the intercept nor the forest's columns' part - is the forest's
# first offset, on every kept row, and the spread of what is left on the rows
# in is the first residual standard deviation. An estimated range is never
# wider than the response's own, and where the fit cannot be computed the
# response's own range stands in with the start that goes with it: an own
# part of zero.
#
# Under an 'offset_type' other than "default" the own part holds what the
# forest's offset holds in the sweeps: the random part alone under "fixef",
# the fixed part alone under "ranef", both under "bart", where the offset is
# not taken off the response, and nothing under "parametric", where no fit is
# made.
bart_range_route <- function(object, bart_range, verbose = 0L, fitter = range_fit)
{
  verbose <- isTRUE(as.integer(verbose) > 0L)
  y <- object$y
  weights <- object$weights
  rows_in <- if (length(weights) > 0L) which(weights > 0) else seq_along(y)
  n <- length(rows_in)
  w <- if (length(weights) > 0L) as.double(weights[rows_in]) else rep_len(1, n)

  has_offset <- length(object$offset) > 0L && any(object$offset != 0)
  offset_type <- if (has_offset) object$offset_type else "default"
  target <- if (has_offset && offset_type != "bart") y - as.vector(object$offset) else y

  supplied <- is.numeric(bart_range)
  route <- if (supplied) "supplied" else bart_range
  reTrms <- object$reTrms
  mixed <- route == "mixed" && length(reTrms) > 0L

  own <- numeric(length(y))
  coefficients <- NULL
  fallback <- NULL
  info <- list(n = n, p = NA_integer_, fit = "none", added = character(),
               dropped = character(), left.out = character(), theta = NULL,
               sigma = NULL, message = NULL)

  if (route != "response" && offset_type != "parametric") {
    if (mixed && !lme4_is_used()) {
      fallback <- "lme4"
      warning("lme4 is not installed, so the BART response range is the response's own, ",
              "not the mixed model's; install lme4, or set 'bart_range'", call. = FALSE)
    } else {
      columns <- range_fixed_block(object, rows_in, w)
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
          if (mixed) fitter(target[rows_in], w, columns$X, range_cut_reTrms(reTrms, rows_in))
          else fitter(target[rows_in], w, columns$X))), error = function(e) e)
        if (!inherits(fit, "error") && !all(is.finite(c(fit$beta, fit$b))))
          fit <- simpleError("an estimate is not finite")
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
          fixed <- if (length(columns$own) == 0L) 0 else
            drop(sweep(object$X[, columns$own_cols, drop = FALSE], 2L,
                       object$X_means[columns$own_cols]) %*% fit$beta[columns$own])
          random <- if (mixed) as.vector(Matrix::crossprod(reTrms$Zt, fit$b)) else 0
          own <- own + switch(offset_type, fixef = random, ranef = fixed, fixed + random)
        }
      }
    }
  }

  range_own <- range(target[rows_in])
  if (is.null(fallback) && info$fit != "none") {
    width <- diff(range((target - own)[rows_in]))
    if (width == 0) fallback <- "undefined"
    else if (width > diff(range_own)) fallback <- "wider"
  }
  # a fallback takes the response's own range whole, starting values included
  if (!is.null(fallback)) own <- numeric(length(y))

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
        else paste("the initial fit is not defined on", n, "rows"),
        "; see 'bart_range'\n", sep = "")
  }

  sigma <- sd((target - own)[rows_in])
  if (!is.finite(sigma) || sigma == 0) sigma <- 1

  list(range = range_bart, route = route, fallback = fallback, info = info, own = own,
       sigma = sigma, coefficients = coefficients)
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
# exactly, or the chains would not share one prior.
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
  if (!identical(leaf_prior$response.scale, pair[2L] - pair[1L]) ||
      !isTRUE(abs(leaf_prior$response.shift - (pair[1L] + pair[2L]) / 2) <= 1e-12))
    stop("internal error: the BART response range was not installed")
  invisible(sampler)
}
