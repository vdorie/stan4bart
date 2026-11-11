\name{print.stan4bartFit}

\alias{print.stan4bartFit}
\alias{summary.stan4bartFit}
\alias{print.summary.stan4bartFit}

\title{Print and Summarize stan4bart Model Fits}
\description{
  Print a brief description of a fitted model, or compute and print a short
  summary of it. Both check the sampling-phase leapfrog-step diagnostic and
  warn when it suggests the model was under-warmed.
}
\usage{
\method{print}{stan4bartFit}(x, \dots)

\method{summary}{stan4bartFit}(object, \dots)

\method{print}{summary.stan4bartFit}(x, \dots)
}
\arguments{
  \item{x, object}{a fitted model resulting from a call to
  \code{\link{stan4bart}}, or for \code{print.summary.stan4bartFit}, the
  result of a call to \code{summary} on one.}
  \item{\dots}{unused, but provided to match the signatures of other
  generics.}
}
\details{
  Both \code{print} and \code{summary} examine
  \code{object$adaptation$mean_leapfrog}, the mean number of leapfrog
  (gradient evaluation) steps per transition during the sampling phase, per
  chain (see \code{\link{stan4bart}}). When the largest per-chain value
  exceeds 30, a one-line message is issued: \dQuote{parametric sampler tuning
  looks poor (mean leapfrog steps per transition = \dots); consider increasing
  warmup.} A well-tuned fit - one run with at least the default warmup, or
  with \code{warmup} at or above roughly 100 for large fits with many
  random-effect levels - should fall well under this value; the message is a
  cue to lengthen \code{warmup} rather than an error. No message is issued for
  a fit with no adaptation record, e.g. one run with \code{warmup = 0}.

  \code{summary} additionally returns an object of class
  \code{summary.stan4bartFit}, printed by \code{print.summary.stan4bartFit},
  giving the formula, family, chain/sample counts, the per-chain warmup
  and sampling mean leapfrog counts, and the fit's R-hats, which are the
  numbers of the mixing check (see \dQuote{Mixing check} in
  \code{\link{stan4bart}}): one compact line such as
  \code{R-hat: sigma 1.00, g sd 1.01, fit 1.00}, each value as it is, with no
  star, bar or note, and \code{NA: could not be checked}, with the reason
  when a computation failed, when a value is \code{NA}. \code{print} is
  unchanged.

  R-hat compares each chain's spread with the spread of the chains pooled, so
  a value above 1 says the chains disagree, not that pooling them is wrong. At
  the defaults BART's chains routinely disagree on \code{sigma}, on the
  standard deviations of the random effects, on the in-sample error and on a
  single row's prediction, so values above 1 there are expected, and no
  single bar says whether a fit is usable. Averaged estimands, such as the
  mean prediction over a subgroup or a contrast between two subgroups, mix far
  better, so read R-hat on the estimand you report. \code{extract} with
  \code{combine_chains = FALSE} gives each chain's draws of the expected
  value; the example computes R-hat of a contrast of subgroup means with them,
  by a short function of its own (\code{posterior::rhat} computes it too).
}
\value{
  \code{print.stan4bartFit} and \code{print.summary.stan4bartFit} return
  their argument, invisibly.

  \code{summary.stan4bartFit} returns a list of class
  \code{summary.stan4bartFit}, printed by \code{print.summary.stan4bartFit}.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link{stan4bart}}
}
\examples{
\donttest{
# R-hat of the estimand reported: the difference between the mean expected
# value of two subgroups, one value per draw and chain
set.seed(1)
n <- 200L
df <- data.frame(g = sample(10L, n, replace = TRUE), x1 = rnorm(n),
                 xb = runif(n))
df$y <- 3 * sin(pi * df$xb) + df$x1 + rnorm(df$g[1:10])[df$g] + rnorm(n)

fit <- stan4bart(y ~ x1 + bart(xb) + (1 | g), df, chains = 2L, iter = 300L,
                 warmup = 150L, verbose = -1, bart_args = list(n.trees = 25L))
summary(fit)

# observations by samples by chains
ev <- extract(fit, "ev", combine_chains = FALSE)
group <- df$xb > 0.5
contrast <- apply(ev[group, , , drop = FALSE], 2:3, mean) -
  apply(ev[!group, , , drop = FALSE], 2:3, mean)
# split R-hat of a samples-by-chains matrix: each chain cut in two, then the
# between-chain and within-chain variances compared (posterior::rhat adds rank
# normalization, which matters for heavy tails)
split.rhat <- function(x) {
  h <- nrow(x) \%/\% 2L
  x <- cbind(x[seq_len(h), ], x[h + seq_len(h), ])
  W <- mean(apply(x, 2L, var))
  sqrt(((h - 1) / h * W + var(colMeans(x))) / W)
}
split.rhat(contrast)
}
}
