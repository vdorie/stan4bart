\name{mvbart}

\alias{mvbart}
\alias{print.mvbartFit}
\alias{mvbartFit}

\title{Multivariate BART for Correlated Outcomes (Seemingly-Unrelated Regressions)}
\description{
  Jointly fits Bayesian additive regression tree (BART) models to two or more
  continuous outcomes that share a predictor set but have correlated residuals.
  Each outcome gets its own \dQuote{sum-of-trees} mean surface; the surfaces are
  coupled through a residual covariance matrix \eqn{\Sigma} that is given a
  conjugate inverse-Wishart prior and sampled jointly. This is the
  seemingly-unrelated regressions (SUR) analogue of BART: modeling the
  cross-outcome residual covariance yields a better-calibrated joint predictive
  distribution than fitting an independent BART per outcome whenever the outcomes
  are genuinely correlated.

  For \eqn{q} outcomes fits:
  \deqn{Y_{ik} = f_k(x_i) + e_{ik}, \quad e_i = (e_{i1}, \ldots, e_{iq}) \sim
    \mathrm{N}_q(0, \Sigma),}{Y_ik = f_k(x_i) + e_ik, e_i ~ N_q(0, Sigma),}
  where each \eqn{f_k} is a BART \link[dbarts:bart]{model} and \eqn{\Sigma} is the
  \eqn{q \times q}{q x q} residual covariance.
}
\usage{
mvbart(formula,
       data = NULL,
       subset,
       na.action = getOption("na.action", "na.omit"),
       test = NULL,
       n.samples = 1000L,
       n.burn = 1000L,
       n.chains = 4L,
       n.trees = 75L,
       prior = NULL,
       seed = NA_integer_,
       verbose = FALSE,
       bart_args = list())
}
\arguments{
  \item{formula}{a \code{\link{formula}} whose left-hand side is a multivariate
  response with at least two columns, e.g. \code{cbind(y1, y2) ~ x1 + x2}. The
  right-hand side gives the shared predictors used by every outcome's BART
  component; factors are expanded by \code{\link[dbarts]{dbarts}}.}
  \item{data}{an optional data frame containing the variables in the formula.}
  \item{subset, na.action}{optional components adjusting the model frame. See
  \code{\link{model.frame}}. \code{na.action} cannot be \code{"na.pass"}.}
  \item{test}{an optional data frame of held-out predictors. When present, the
  per-outcome forest predictions \eqn{f_k(x^*)} are stored for each retained
  draw.}
  \item{n.samples}{positive integer; number of retained posterior draws per
  chain.}
  \item{n.burn}{non-negative integer; number of burn-in sweeps discarded per
  chain.}
  \item{n.chains}{positive integer; number of independent chains, run
  sequentially.}
  \item{n.trees}{positive integer; number of trees in each outcome's forest.}
  \item{prior}{an optional named list specifying the inverse-Wishart prior on
  \eqn{\Sigma}, with elements \code{df} (degrees of freedom \eqn{\nu_0}) and
  \code{scale} (a \eqn{q \times q}{q x q} positive-definite scale matrix
  \eqn{\Psi_0}). Either may be omitted. The prior has mean
  \eqn{\Psi_0 / (\nu_0 - q - 1)}{Psi0 / (df - q - 1)}. Defaults to a
  weakly-informative prior: \code{df} \eqn{= q + 2} (the smallest integer giving a
  finite mean) and \code{scale} chosen so that the prior mean of \eqn{\Sigma}
  equals the diagonal matrix of marginal outcome variances.}
  \item{seed}{optional integer. If supplied, per-chain seeds are drawn
  deterministically from it and results are reproducible; the caller's RNG state
  is restored afterwards.}
  \item{verbose}{logical or integer; if a positive number, prints progress every
  \code{verbose} sweeps.}
  \item{bart_args}{an optional named list of further arguments forwarded to each
  outcome's BART sampler. Recognized \code{\link[dbarts]{dbartsControl}}
  arguments not already set by \code{mvbart} (such as \code{n.cuts} or
  \code{keepTrees}) are applied to the shared control; the end-node prior scale
  \code{k} is also honored.}
}
\details{
  Conditional on the other equations' current residuals, updating one outcome's
  forest is an ordinary univariate BART step: the working outcome is Gaussian
  with a per-observation offset (the conditional mean implied by \eqn{\Sigma})
  and a scalar conditional variance. These are injected into the outcome's
  \code{\link[dbarts]{dbarts}} sampler through its \code{setOffset} and
  \code{setSigma} methods; the sampler is created with a \code{fixed} residual
  prior so the engine never draws its own nugget. After the forest sweep,
  \eqn{\Sigma} is drawn from its inverse-Wishart full conditional,
  \eqn{\Sigma \mid \cdot \sim \mathrm{IW}(\nu_0 + n, \Psi_0 + E'E)}{Sigma | . ~
  IW(df + n, scale + E'E)}, where \eqn{E} is the matrix of structural residuals.
  The routine is built entirely on the public \code{dbarts} R interface and adds
  no compiled code.

  \subsection{Diagnostics}{
    Because the response passed to each sampler is fixed (only the offset moves),
    the creation-time \eqn{[-0.5, 0.5]} response rescaling used internally by
    BART stays anchored across sweeps. The returned \code{mMax} records, per
    retained draw, the largest per-sweep \eqn{\max_i |m_{ik}| / \mathrm{range}(Y_k)}{
    max_i |m_ik| / range(Y_k)}; values well below 1 indicate the fixed scale is
    well anchored. If the maximum drift exceeds 1 a warning is issued: for such
    heavy-tailed cases the recommended fallback is a longer warmup, or a
    warmup-only response rescale (re-anchoring the internal scale to include the
    offset range during burn-in only).
  }
  \subsection{Caveat: residual correlation bias with shared predictors}{
    When all outcomes share the same predictors, BART's flexibility introduces a
    small positive bias in the recovered residual correlation (on the order of
    \eqn{+0.01} to \eqn{+0.02} at moderate to high correlation), because
    finite-sample mean-fitting error is correlated across equations. This makes
    interval coverage of the true correlation mildly liberal. The bias shrinks
    toward zero as \eqn{n} grows and is inherent to the nonparametric mean model,
    not a defect of the sampler.
  }
  \subsection{Future extension}{
    \code{mvbart} is a pure-SUR sampler with no random-effect component.
    Combining the cross-outcome residual covariance modeled here with
    \code{\link{stan4bart}}'s multilevel machinery (SUR \eqn{\times}{x} random
    effects, i.e. a multivariate mixed model) is a natural extension that is not
    yet implemented.
  }
}
\value{
  Returns a list assigned class \code{mvbartFit} with components:
  \item{\code{Sigma}}{a \eqn{q \times q \times}{q x q x} \code{n.samples}
  \eqn{\times}{x} \code{n.chains} array of residual covariance draws, with the
  outcome names on the first two dimensions.}
  \item{\code{train}}{an \eqn{n \times q \times}{n x q x} \code{n.samples}
  \eqn{\times}{x} \code{n.chains} array of per-outcome forest fits
  \eqn{f_k(x_i)} on the training data.}
  \item{\code{test}}{when \code{test} was supplied, a
  \code{n.test} \eqn{\times q \times}{x q x} \code{n.samples} \eqn{\times}{x}
  \code{n.chains} array of per-outcome forest predictions \eqn{f_k(x^*)};
  otherwise \code{NULL}.}
  \item{\code{mMax}}{a \code{n.samples} \eqn{\times}{x} \code{n.chains} matrix of
  the per-sweep scale-drift diagnostic.}
  \item{\code{y}}{the response matrix.}
  \item{\code{outcome.names}}{character vector of outcome names.}
  \item{\code{n.samples}, \code{n.burn}, \code{n.chains}}{sampler settings.}
  \item{\code{prior}}{the inverse-Wishart prior actually used, a list with
  \code{df} and \code{scale}.}
  \item{\code{call}, \code{formula}, \code{terms}}{the matched call, formula, and
  model terms.}
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link{stan4bart}}, \code{\link[dbarts]{bart}},
  \code{\link[dbarts]{dbarts}}
}
\examples{
# two correlated outcomes sharing three predictors
set.seed(99)
n <- 200L
x <- matrix(runif(n * 3L), n, 3L)
colnames(x) <- paste0("x", 1:3)
B <- matrix(c(1.5, -1.0, 0.7, -0.6, 1.1, -0.9), 3L, 2L)
rho <- 0.7; s1 <- 1.0; s2 <- 1.3
Sigma <- matrix(c(s1^2, rho * s1 * s2, rho * s1 * s2, s2^2), 2L, 2L)
Y <- x \%*\% B + matrix(rnorm(n * 2L), n, 2L) \%*\% chol(Sigma)
df <- data.frame(y1 = Y[, 1L], y2 = Y[, 2L], x)

fit <- mvbart(cbind(y1, y2) ~ x1 + x2 + x3, data = df,
              n.samples = 200L, n.burn = 100L, n.chains = 1L,
              n.trees = 50L, seed = 1L)

# posterior of the residual correlation
Sig <- fit$Sigma
rho.draws <- Sig[1L, 2L, , ] / sqrt(Sig[1L, 1L, , ] * Sig[2L, 2L, , ])
c(mean = mean(rho.draws), quantile(rho.draws, c(0.025, 0.975)))
}
\keyword{models}
\keyword{regression}
