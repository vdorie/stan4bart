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

# Named vector of the R-hats of a fit's checked quantities: "sigma" unless
# the response is binary, "<factor> sd" (with the term appended for a factor
# with several) per grouping factor, and "fit".
fit_mixing <- function(object)
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

  # the response and the expected value, on the rows extract() reports
  if (!is.null(object$bart_train) || bart_needs_recompute(object, "train")) {
    ev <- extract(object, "ev", combine_chains = FALSE)
    n_draws <- dim(ev)[2L:3L]
    resid <- ev - pad_kept_rows(object, object$y, NA_real_)
    dim(resid) <- c(dim(ev)[1L], prod(n_draws))
    rmse <- sqrt(colMeans(resid^2, na.rm = TRUE))
    result["fit"] <- split_rhat(matrix(rmse, n_draws[1L], n_draws[2L]))
  }
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
