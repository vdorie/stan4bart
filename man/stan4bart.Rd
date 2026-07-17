\name{stan4bart}

\alias{stan4bart}

\title{Semiparametric Models Using a Gradient-Based Sampler and BART}
\description{
  This function fits semi-parametric linear and probit models that have a
  non-parametric, BART component and one or more of a parametric fixed
  effect (unmodeled coefficients), or a parametric random effect
  (modeled coefficients). If
  \eqn{f(x)} is a BART \dQuote{sum-of-trees} \link[dbarts:bart]{model}, fits:
\itemize{
  \item For continuous response variables:
  \deqn{Y \mid b \sim \mathrm{N}\left(f(X^b) + X^f\beta + Zb, \sigma^2\right) \\
  b \sim \mathrm{N}(0, \Sigma_b)}{Y | b ~ N(f(X^b) + X^f beta + Zb, sigma^2)
  b ~ N(0, Sigma_b) }
  \item For binary response variables:
  \deqn{P(Y = 1 \mid b) = \Phi\left(f(X^b) + X^f\beta + Zb\right) \\
  b \sim \mathrm{N}(0, \Sigma_b)}{P(Y = 1 | b) = \Phi(f(X^b) + X^f\beta + Zb)
  b ~ N(0, Sigma_b) }
}}
\usage{
stan4bart(formula,
          data = NULL,
          subset,
          weights,
          na.action = getOption("na.action", "na.omit"),
          offset,
          contrasts = NULL,
          test = NULL,
          treatment = NULL,
          offset_test = NULL,
          verbose = FALSE,
          iter = 2000L,
          warmup = iter \%/\% 2L,
          skip = 1L,
          chains = 4L,
          cores = getOption("mc.cores", 1L),
          refresh = max(iter \%/\% 10L, 1L),
          offset_type = c("default", "fixef", "ranef", "bart", "parametric"),
          seed = NA_integer_,
          keep_fits = TRUE,
          save_warmup = FALSE,
          store = c("fits", "trees"),
          callback = NULL,
          stan_args = NULL,
          bart_args = NULL)
}
\arguments{
  \item{formula}{a \code{\link{formula}} object, or one that can be coerced to
  that type. Terms on the right-hand-side of the formula that are encased in
  a symbolic call to \code{bart()} will be used to create the non-parametric
  component. Terms that use the \code{\link[lme4]{lmer}}-style grouping syntax
  will be added as parametric, hierarchical varying intercepts and slopes. All
  other terms will be added as fixed effects.}
  \item{data}{an optional data frame containing the variables in the formula.
  Its use is strongly encouraged.}
  \item{subset, weights, na.action, offset, contrasts}{optional components
  adjusting the constructed model matrices and otherwise changing the linear
  predictor. \code{na.action} cannot be \code{"na.pass"}. See
  \code{\link{lm}} and \code{\link{model.matrix.default}}.}
  \item{test}{an optional data frame to be used as test data. If present, the
  test predictions will be stored as the sampler runs and can be extracted
  later.}
  \item{treatment}{an optional symbol, that when present and refers to a binary
  variable, will be used to create a test data frame with the treatment variable
  set to its counterfactual. Only one of \code{test} and \code{treatment} can
  be supplied.}
  \item{offset_test}{optional vector which will be added to the test
  predictions.}
  \item{verbose}{a logical or integer. If \code{FALSE} or non-positive, runs
  quietly. Additional levels of information may be displayed for increasingly
  positive numbers, however a large number of diagnostics are suppressed when
  running multi-threaded. If negative, all diagnostic information is ignored.}
  \item{iter}{positive integer indicating the number of posterior samples to
  draw and return. Includes warmup samples.}
  \item{warmup}{non-negative integer indicating number of posterior samples to
  draw and throw away. Also controls the number of iterations for which
  adaptation is run. The default (half of \code{iter}) is ample for the
  parametric sampler's step size and mass tuning to converge. If \code{iter}
  is shortened, keep \code{warmup} at or above roughly 100 for large-\eqn{n}
  fits with many random-effect levels: under-warming leaves the tuning poor
  and can roughly double the per-iteration cost of the parametric step. The
  diagnostic is \code{fit$adaptation$mean_leapfrog} (mean leapfrog steps per
  transition during sampling; see \sQuote{Value} below) - a high value
  relative to a well-tuned fit indicates under-warming, and \code{print} and
  \code{summary} warn when it looks poor.}
  \item{skip}{one or two positive integers. Every \code{skip} sample will be
  kept, while every other sample will be discarded. If argument is length two,
  an attempt will be made to use he named element \code{"bart"} for BART and
  \code{"stan"} for Stan. If not named, BART is the first skip element and Stan
  is the second. This argument does not impact the number of \code{iter}s
  returned, unlike a conventional \dQuote{thinning} parameter.}
  \item{chains}{positive integer giving the number of Markov Chains to sample.}
  \item{cores}{positive integer giving the number of units of parallelization.
  Computation for each chain will be divide among the cores available. When 
  greater than one, verbose output within chains will not be available.}
  \item{refresh}{positive integer giving the frequency with which diagnostic
  information should be printed as the sampler runs. Only applies with
  \code{cores} (or \code{chains}) equal to 1.}
  \item{offset_type}{character; an experimental/testing feature that controls
  how \code{offset} is to be interpreted. When one of \code{"fixef"},
  \code{"ranef"}, or \code{"bart"}, the offset is used to replace that part of
  the model. When \code{"parametric"}, it replaces both of the fixed and random
  parametric components. Sampling is still done for these components and their
  draws are stored, however whenever they were present in the fit the
  supplied value is used instead.}
  \item{seed}{Optional integer specifying the desired pRNG \link[=set.seed]{seed}.
  It should not be needed when running single-threaded - calling
  \code{\link{set.seed}} will suffice. The primary use of \code{seed} is to
  obtain reproducible results when multi-threaded. See Reproducibility section
  below.}
  \item{keep_fits}{Logical that, when false, prevents the sampler from
  storing each draw. Intended to be used with \code{callback}.}
  \item{save_warmup}{Logical; when \code{FALSE} (the default) the warmup draws
  are not stored. In their place, per chain, the frozen sampler tuning (step
  size and diagonal inverse mass), a warmup-end position snapshot, and a
  thinned trace of the monitored scalars are recorded in the \code{adaptation}
  component of the fit; this is what the 2026 Stan/mixed-model toolchain does
  and cuts the fit object roughly in half at scale (see the storage-policy
  note in \code{docs/design/walnuts.md}). When \code{TRUE}, the full per-draw
  warmup is stored under \code{warmup} as before, and
  \code{include_warmup = TRUE} accessors work; on a default fit those
  accessors raise an informative error.}
  \item{store}{character; how the fitted BART component is retained.
  \code{"fits"} (the default) stores the full \eqn{n \times}{n x}
  draws \eqn{\times}{x} chains \code{bart_train} (and \code{bart_test})
  blocks, as in prior releases. \code{"trees"} instead stores only the kept
  trees and recomputes those blocks on demand through the \code{dbarts}
  predict path; it implies \code{keepTrees} (and errors if \code{bart_args}
  contradicts it with \code{keepTrees = FALSE}). This is an opt-in memory/time
  trade: at rest the fit no longer carries the \code{bart_train} block, which is
  94-98\% of a large fit object (\code{n >= 10000}); in exchange every reader of
  the BART fits pays a recompute cost on access. \code{fitted} streams the
  posterior mean in row blocks, so its peak memory stays bounded, but
  \code{extract} must materialize the whole block (unavoidable) - roughly 15 s
  per call at \code{n = 10000} with 4 chains x 1000 draws, and over a minute at
  \code{n = 50000}. Storing trees costs essentially nothing during sampling. The
  saving only exists above a crossover near \code{n = 550}: below it the kept
  trees cost as much as the fits they replace and add latency, so \code{"fits"}
  remains the default and the recommended choice at small \code{n}. Under
  \code{"trees"} the \code{bart_train}/\code{bart_test} components are absent;
  read the BART fits through \code{\link[=stan4bart-generics]{extract}},
  \code{\link[=stan4bart-generics]{fitted}}, or
  \code{\link[=stan4bart-generics]{predict}} instead.}
  \item{callback}{A function that will be called at each iteration, accepting
  three arguments: \code{yhat.train}, \code{yhat.test}, \code{stan_pars}.
  See details for more information.}
  \item{stan_args}{optional list, specifying further arguments to the
  parametric (gradient-sampled) component. See details below.}
  \item{bart_args}{optional list, specifying further arguments to BART. See
  details below.}
}
\details{
  Fits a Bayesian \dQuote{mixed effect} model with a non-parametric Bayesian
  Additive Regression Trees (BART) component. For continuous responses:

  \itemize{
    \item \eqn{Y_i \mid b \sim \mathrm{N}\left(f(X^b_i) + X^f_i\beta +
    Z_ib_{g[i]}, \sigma^2\right)}{Y_i | b ~ N(f(X^b_i) + X^f_i\beta + Z_ib_{g[i]},
    \sigma^2)}
    \item \eqn{b_j \sim \mathrm{N}(0, \Sigma_b)}{b_j ~ N(0, \Sigma_b)}
  }

  where \eqn{b_j} are the \dQuote{random effects} - random intercepts and
  slopes - that correspond to group \eqn{j}, \eqn{g[i]} is a mapping from
  individual \eqn{i} to its group index,  \eqn{f} - a BART sum-of-trees model,
  \eqn{X^b} are predictors used in the BART model, \eqn{X^f} are predictors
  in a parametric, linear \dQuote{fixed effect} component, \eqn{Z} is the 
  design matrix for the random intercept and slopes, and \eqn{sigma} and
  \eqn{Sigma_b} are variance components.
  
  Binary outcome models are obtained by assuming a latent variable that has
  the above distribution, and that the observed response is 1 when that
  variable is positive and 0 otherwise. The response variable marginally has
  the distribution:
  
  \itemize{
    \item \deqn{P(Y_i = 1 \mid b) = \Phi\left(f(X^b_i) + X^f_i\beta +
    Z_ib_{g[i]}\right)}{P(Y_i = 1 | b) = \Phi(f(X_i^b) + X^f_i\beta +
    Z_ib_{g[i]})}
  }
  
  where \eqn{\Phi} is the cumulative distribution function of the standard
  normal distribution.

  
  \subsection{Terminology}{
    As \code{stan4bart} fits a Bayesian model, essentially all components are
    \dQuote{modeled}. Furthermore, as it has two first-level, non-hierarchical
    components, \dQuote{fixed} effects are ambiguous. Thus we adopt:
    \itemize{
      \item \dQuote{fixed} - refers only to the parametric, linear, individual
      level mean component, \eqn{X^f\beta}; these are \dQuote{unmodeled
      coefficients} in other contexts
      \item \dQuote{random} - refers only to the parametric, linear,
      hierarchical mean component, \eqn{Zb}; these are \dQuote{modeled
      coefficients} in other contexts
      \item \dQuote{bart} - refers only to the nonparametric, individual level
      mean component, \eqn{f(X^b)}
    }
  }
  \subsection{Model Specification}{
    Model specification occurs in the \code{formula} object, according to the
    following rules:
    
    \itemize{
      \item variables or terms specified inside a pseudo-call to \code{bart}
      are used for the \dQuote{bart} component, e.g. \code{y ~ bart(x_1 + x_2)}
      \item variables or terms specified according to \code{\link[lme4]{lmer}}
      syntax are used for the \dQuote{random} effect component, e.g.
      \code{y ~ (1 | g_1) + (1 + x_3 | g_1)}
      \item remaining variables not inside a \code{bart} or \dQuote{bars}
      construct are used for the \dQuote{fixed} effect component; e.g.
      \code{y ~ x_4}
    }
    
    All three components can be present in a single model, however are
    \code{bart} part must present. If you wish to fit a model without one, use
    \code{stan_glmer} in the \code{rstanarm} package instead.

    Factor variables are handled differently in each component. In the
    \code{fixed} and \code{random} components they use the usual
    \code{\link{model.matrix}} contrasts and \code{lme4} grouping-factor
    machinery, respectively. In the \code{bart} component a factor is encoded
    as \link[dbarts]{dbarts}'s categorical splits: a single column whose tree
    prior chooses level subsets to split on, rather than one indicator column
    per level. Consequently \code{extract(type = "varcount")} reports one row
    per factor, named by the factor's variable name, not one row per level.
    A factor level in \code{newdata} (to \code{predict}) or \code{test} that
    was not observed when fitting is an error for a \code{bart} factor; this
    is unlike the \code{random} component, where an unseen grouping level is
    handled per \code{sample_new_levels} in \code{\link{extract}}.
  }
  \subsection{Additional Arguments}{
    The \code{stan_args} and \code{bart_args} arguments to \code{stan4bart} can
    be used to pass further arguments to the parametric and BART components
    respectively. \code{stan_args} is named for historical reasons - the
    parametric component was originally sampled with Stan's NUTS; it is now
    sampled with a gradient-based WALNUTS sampler driven by a hand-derived
    analytic log-posterior gradient, and \code{stan_args} configures that
    sampler and its priors. \code{bart_args} is similar to
    \code{\link[dbarts]{bart}}, but not identical, as \code{stan4bart}
    constructs its own model internally.

    Prior- and model-related \code{stan_args} include:
    \itemize{
      \item \code{prior_covariance} - always a \code{decov()} prior on the
      random effect covariance.
      \item \code{prior} - the fixed effect coefficient prior; supports
      \code{normal}, \code{student_t}, and \code{cauchy}. The shrinkage
      families (\code{hs}, \code{hs_plus}, \code{lasso}, \code{laplace},
      \code{product_normal}) are not supported by the gradient-based sampler
      and raise an error if supplied.
      \item \code{prior_aux} - the residual sd prior (continuous responses
      only); supports \code{normal}, \code{student_t}, \code{cauchy}, and
      \code{exponential}.
      \item \code{prior_intercept} - accepted for backwards compatibility, but
      a warning is raised and it is dropped; BART already carries the
      intercept.
      \item \code{QR} - whether to use a QR decomposition of the fixed effect
      design matrix; sampler-agnostic.
      \item \code{init_r} - half-width of the uniform interval the sampler's
      initial unconstrained position is drawn from.
      \item \code{adapt_delta} - the WALNUTS step-size acceptance-rate target, a
      value in \eqn{(0, 1)} (default \code{0.8}); the analog of Stan's
      \code{adapt_delta}. A higher value targets a smaller step size, so warmup
      spends more gradient evaluations per transition but tracks the posterior
      geometry more finely; a lower value takes larger, cheaper steps.
      \item \code{save_raw_parameters} - logical, default \code{FALSE}. When
      \code{FALSE} the stored parametric draw is only the transformed block
      every consumer reads (\code{beta}, \code{b}, \code{theta_L}, and
      \code{aux}) plus the two live diagnostics (\code{lp__}, \code{stepsize__});
      the raw unconstrained rows (\code{z_beta}, \code{z_b}, \code{z_T},
      \code{rho}, \code{zeta}, \code{tau}, \code{aux_unscaled}) and the five
      constant-zero placeholder diagnostic rows are dropped. Set \code{TRUE} to
      restore the full historical Stan-layout rows, e.g. as funnel coordinates
      for forensics on the decov covariance prior.
    }
    See the help page for \code{stan_glmer} in the \code{rstanarm} package for
    the prior constructors (\code{normal}, \code{student_t}, \code{cauchy},
    \code{exponential}, \code{decov}).

    The remaining \code{stan_args} - \code{adapt_gamma},
    \code{adapt_kappa}, \code{adapt_t0}, \code{adapt_init_buffer},
    \code{adapt_term_buffer}, \code{adapt_window}, \code{stepsize},
    \code{stepsize_jitter}, and \code{max_treedepth} - configured Stan's NUTS
    adaptation and have no analog under WALNUTS. They are still accepted so
    that existing scripts do not error, but are ignored, and a warning naming
    each supplied argument is raised at fit time; they are deprecated and will
    be removed in a future release.

    BART arguments include:
    \itemize{
      \item further arguments to \code{\link[dbarts]{dbartsControl}} that are not
      specified by \code{stan4bart}, such as \code{keepTrees} or
      \code{n.trees}; keeping trees can be costly in terms of memory, but is
      required to use \code{\link[=stan4bart-generics]{predict}}
    }
  }
  \subsection{Reproducibility}{
    Behavior differs when running multi- and single-threaded, as the pseudo random
    number generators (pRNG) used by R are not thread safe. When single-threaded,
    R's built-in generator is used; if set at the start, \code{\link{.Random.seed}}
    will be used and its value updated as samples are drawn. When multi-threaded,
    the default behavior is draw new random seeds for each thread using the clock
    and use thread-specific pRNGs.
    
    This behavior can be modified by setting \code{seed}. For the single-threaded
    case, that seed will be installed and the existing seed replaced at the end,
    if applicable. For multi-threaded runs, the seeds for threads are drawn
    sequentially using the supplied seed, and will not change the state of 
    R's built-in generator.
    
    Consequently, the \code{seed} argument should not be needed when running
    single-threaded - \code{\link{set.seed}} will suffice. When multi-threaded,
    \code{seed} can be used to obtain reproducible results.
  }
  \subsection{Callbacks}{
    Callbacks can be used to avoid expensive memory allocation, aggregating
    results as the sampler proceeds instead of waiting until the end. A callback
    funtion must accept the arguments:
    \itemize{
      \item \code{yhat.train} - BART predictions of the expected value of the
      training data, conditioned on the Stan parameters
      \item \code{yhat.test} - when applicable, the same values as above but
      for the test data; \code{NULL} otherwise
      \item \code{stan_pars} - named vector of parametric-component samples,
      including fixed and random effects and variance parameters
    }
    It is expected that the callback will return a vector of the same length
    each time. If not, invalid memory writes will likely result. If the
    result of the callback has names, they will be added to the result.
  }
  \subsection{Serialization}{
    At present, \code{stan4bart} models cannot be safely \code{\link{save}}d
    and \code{\link{load}ed} in a way that the sampler can be restarted. This
    feature may be added in the future.
  }
}
\value{
  Returns a list assigned class \code{stan4bartFit}. Has components below, some
  of which will be \code{NULL} if not applicable.
  
  Input values:
  \item{\code{y}}{response vector}
  \item{\code{weights}}{weights vector or null}
  \item{\code{offset}}{offset vector or null}
  \item{\code{frame}}{joint model frame for all components}
  \item{\code{formula}}{formula used to specify the model}
  \item{\code{na.action}}{supplied na.action}
  \item{\code{call}}{original call}
  
  Stored data:
  \item{\code{bartData}}{\code{\link[dbarts:dbartsData]{data}} object used for BART
  component}
  \item{\code{X}}{fixed effect design matrix or NULL}
  \item{\code{X_means}}{column means of fixed effect design matrix when
  appropriate}
  \item{\code{reTrms}}{random effect \dQuote{terms} object when applicable, as
  used by \code{\link[lme4]{lmer}}}
  \item{\code{test}}{named list when applicable, having components \code{X} and
  \code{reTrms}; test data for BART is added to the \code{bartData} result}
  \item{\code{treatment}}{treatment vector, when applicable}
  
  Results, better accessed using \code{\link[=stan4bart-generics]{extract}}:
  \item{\code{bart_train}}{samples of individual posterior predictions for BART
  component; absent (recomputed from the kept trees on access) when
  \code{store = "trees"}}
  \item{\code{bart_test}}{predicted test values for BART component, when
  applicable; likewise absent under \code{store = "trees"}}
  \item{\code{bart_varcount}}{BART variable counts}
  \item{\code{sigma}}{samples of residual standard error; not present for
  binary outcomes}
  \item{\code{k}}{samples of the end-node sensitivity parameter; only present when
  it is modeled}
  \item{\code{ranef}}{samples of random effects, or modeled coefficients; will
  be a named list, with effects for each grouping factor}
  \item{\code{Sigma}}{samples of covariance of random effects; also a named list
  with one element for each grouping factor}
  \item{\code{fixef}}{samples of the fixed effects, or unmodeled coefficients}
  \item{\code{callback}}{samples returned by an optional callback function}
  
  Other items:
  \item{\code{warmup}}{a list of warmup samples, containing the same objects in
  the results subsection; present only when \code{save_warmup = TRUE}}
  \item{\code{adaptation}}{present when warmup was run and \code{save_warmup} is
  \code{FALSE}. A list with, per chain, the frozen \code{step_size}, the diagonal
  \code{inv_mass}, the mean number of leapfrog (gradient evaluation) steps per
  transition during the sampling phase (\code{mean_leapfrog}) and during
  warmup (\code{mean_leapfrog_warmup}) - the runtime-relevant tuning
  diagnostic; a high \code{mean_leapfrog} relative to a well-tuned fit
  indicates under-warming (see the \code{warmup} argument above) - a
  warmup-end \code{snapshot} of the parametric draw, and a \code{trace} (a
  thinned warmup trace of the monitored scalars: the parametric \code{stan}
  rows and the BART \code{sigma}) for a qualitative where-did-it-stabilize
  plot}
  \item{\code{sampler.bart}}{external points to BART samplers; used only for
  \code{\link[=stan4bart-generics]{predict}} when \code{keepTrees} is \code{TRUE}}
  \item{\code{range.bart}}{internal scale used by BART samplers, used by 
  \code{\link[=stan4bart-generics]{predict}} when \code{keepTrees} is \code{TRUE}}
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link[dbarts]{bart}}, \code{\link[lme4]{lmer}}, and \code{stan_glmer} in the 
  \code{rstanarm} package
}
\examples{
# simulate data (extension of Friedman MARS paper)
# x consists of 10 variables, only first 5 matter
# x_4 is linear
f <- function(x)
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
      10 * x[,4] + 5 * x[,5]

set.seed(99)
sigma <- 1.0

n <- 100
n.g.1 <- 5L
n.g.2 <- 8L

# sample observation level covariates and calculate marginal mean
x <- matrix(runif(n * 10), n, 10)
mu.bart  <- f(x) - 10 * x[,4]
mu.fixef <- 10 * x[,4]

# varying intercepts and slopes for first grouping factor
g.1 <- sample(n.g.1, n, replace = TRUE)
Sigma.b.1 <- matrix(c(1.5^2, .2, .2, 1^2), 2)
b.1 <- matrix(rnorm(2 * n.g.1), n.g.1) \%*\% chol(Sigma.b.1)

# varying intercepts for second grouping factor
g.2 <- sample(n.g.2, n, replace = TRUE)
Sigma.b.2 <- as.matrix(1.2)
b.2 <- rnorm(n.g.2, 0, sqrt(Sigma.b.2))

mu.ranef <- b.1[g.1,1] + x[,4] * b.1[g.1,2] + b.2[g.2]

y <- mu.bart + mu.fixef + mu.ranef + rnorm(n, 0, sigma)

df <- data.frame(y, x, g.1, g.2)

fit <- stan4bart(
    formula = y ~
        X4 +                         # linear component ("fixef")
        (1 + X4 | g.1) + (1 | g.2) + # multilevel ("ranef")
        bart(. - g.1 - g.2 - X4),    # use bart for other variables
    verbose = -1, # suppress ALL output
    # low numbers for illustration
    data = df,
    chains = 1, iter = 10, bart_args = list(n.trees = 5))

# posterior means of individual expected values
y.hat <- fitted(fit)

# posterior means of the random effects
ranef.hat <- fitted(fit, type = "ranef")
}
