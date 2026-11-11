# End-of-fit mixing check. The split rank-normalized R-hat of Vehtari, Gelman,
# Simpson, Carpenter and Buerkner (2021), as posterior::rhat computes it, of
# the few quantities a fit is judged by: sigma, each grouping factor's
# standard deviations, and the in-sample RMSE of the expected value. Fixed
# effects are redundantly parameterized with the BART component and group
# effects are many, so neither is checked.

# R-hat above which summary() flags a quantity as not mixed.
mixing_rhat_bar <- 1.1

# Each chain cut in two, the middle draw of an odd length dropped.
split_halves <- function(x)
{
  n <- nrow(x)
  h <- n %/% 2L
  cbind(x[seq_len(h), , drop = FALSE], x[seq.int(n - h + 1L, n), , drop = FALSE])
}

# Plain R-hat of the columns of x, each a chain.
plain_rhat <- function(x)
{
  n <- nrow(x)
  w <- mean(apply(x, 2L, var))
  sqrt((n * var(colMeans(x)) / w + n - 1) / n)
}

normal_scores <- function(x)
{
  z <- qnorm((rank(x, ties.method = "average") - 3 / 8) / (length(x) + 1 / 4))
  dim(z) <- dim(x)
  z
}

# Split rank-normalized R-hat of a draws-by-chains matrix: the larger of the
# bulk value and the value for the draws folded about their median. NA when
# the draws cannot be compared: a missing or infinite value, a constant
# quantity, or fewer than four draws per chain.
split_rhat <- function(x)
{
  x <- as.matrix(x)
  if (nrow(x) < 4L || !all(is.finite(x)) || max(x) - min(x) < .Machine$double.eps)
    return(NA_real_)
  result <- max(plain_rhat(normal_scores(split_halves(x))),
                plain_rhat(normal_scores(split_halves(abs(x - median(x))))))
  # a folded column that is constant has no within-chain variance: 0 / 0
  if (is.nan(result)) NA_real_ else result
}

# The in-sample RMSE of the expected value against the response, per draw over
# the chains, the sums of squares taken in blocks of rows as fitted() streams
# them so that memory stays bounded.
fit_rmse_draws <- function(object, block_cells)
{
  M <- prod(bart_sample_chain_dims(object))
  y <- object$y
  sums <- numeric(M)
  n_rows <- length(y)
  streamed <- visit_ev_blocks(object, "ev", "train", TRUE, function(ev_b, idx, M) {
    resid <- matrix(ev_b, length(idx), M) - y[idx]
    sums <<- sums + colSums(resid^2)
  }, block_cells)
  if (!streamed) {
    # the rows extract() reports, those a component dropped missing
    ev <- extract(object, "ev", combine_chains = FALSE)
    resid <- ev - pad_kept_rows(object, y, NA_real_)
    dim(resid) <- c(dim(ev)[1L], M)
    sums <- colSums(resid^2, na.rm = TRUE)
    n_rows <- sum(!is.na(resid[, 1L]))
  }
  sqrt(sums / n_rows)
}

# R-hat of the per-draw RMSE. NA where it cannot be computed, and where
# computing it fails, as for want of memory: a check never loses a fit.
fit_rhat <- function(object, block_cells)
{
  tryCatch({
    n_draws <- bart_sample_chain_dims(object)
    split_rhat(matrix(fit_rmse_draws(object, block_cells), n_draws[1L], n_draws[2L]))
  }, error = function(e) NA_real_)
}

# Named vector of the R-hats of a fit's checked quantities: "sigma" unless
# the response is binary, "<factor> sd" (with the term appended for a factor
# with several) per grouping factor, and "fit".
fit_mixing <- function(object, block_cells = 2e6)
{
  result <- numeric()
  # keep_fits = FALSE keeps no draws to read
  if (is.null(object$stan)) return(result)
  if (!identical(object$family$family, "binomial"))
    result["sigma"] <- split_rhat(extract(object, "sigma", combine_chains = FALSE))

  if (any(startsWith(dimnames(object$stan)[[1L]], "theta_L."))) {
    Sigmas <- extract(object, "Sigma", combine_chains = FALSE)
    for (j in seq_along(Sigmas)) {
      Sigma <- Sigmas[[j]]
      terms <- dimnames(Sigma)[[1L]]
      for (k in seq_len(dim(Sigma)[1L])) {
        label <- paste0(names(Sigmas)[j], " sd", if (length(terms) > 1L) paste0(" ", terms[k]))
        result[label] <- split_rhat(matrix(sqrt(Sigma[k, k, , ]), dim(Sigma)[3L], dim(Sigma)[4L]))
      }
    }
  }

  # the response and the expected value
  if (!is.null(object$bart_train) || bart_needs_recompute(object, "train"))
    result["fit"] <- fit_rhat(object, block_cells)
  result
}

# The R-hat line of a summary, its quantities each followed by a star when over
# the bar, and the lines that explain what a star or an NA says.
format_mixing <- function(mixing, n_chains)
{
  flagged <- !is.na(mixing) & mixing > mixing_rhat_bar
  values <- paste0(ifelse(is.na(mixing), "NA", formatC(mixing, digits = 2L, format = "f")),
                   ifelse(flagged, "*", ""))
  c(paste0("R-hat: ", paste(names(mixing), values, collapse = ", ")),
    if (any(flagged))
      sprintf("* over %.1f: the %s disagree", mixing_rhat_bar,
              if (n_chains > 1L) "chains" else "halves of the chain"),
    if (anyNA(mixing)) "NA: could not be checked")
}

# The R-hats of an mvbart fit: each outcome's residual standard deviation (the
# square root of the covariance diagonal) and its in-sample RMSE, the forest's
# fit against the response, per draw.
mvbart_mixing <- function(Sigma, train, y, outcome_names)
{
  q <- length(outcome_names)
  d <- dim(Sigma)[3L:4L]
  sds <- numeric(q)
  names(sds) <- paste(outcome_names, "sd")
  fits <- numeric(q)
  names(fits) <- paste(outcome_names, "fit")
  for (k in seq_len(q)) {
    sds[k] <- split_rhat(matrix(sqrt(Sigma[k, k, , ]), d[1L], d[2L]))
    resid <- train[, k, , ] - y[, k]
    dim(resid) <- c(nrow(y), prod(d))
    fits[k] <- split_rhat(matrix(sqrt(colMeans(resid^2)), d[1L], d[2L]))
  }
  c(sds, fits)
}
