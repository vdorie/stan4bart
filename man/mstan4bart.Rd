\name{mstan4bart}

\alias{mstan4bart}

\title{Multilevel Models Using Stan and BART}
\description{
  This function fits hierarchical/multilevel models that have a non-parametric,
  BART component and one or more of a parametric fixed effect (unmodeled
  coefficients), or a parametric random effect (modeled coefficients). If
  \eqn{f(x)} is a BART \dQuote{sum-of-trees} \link[=bart]{model}, fits:
\itemize{
  \item For continuous response variables:
  \deqn{Y \mid b \sim \mathrm{N}\left(f(X^b) + X^f\beta + Zb, \sigma^2\right) \\
  b \sim \mathrm{N}(0, \Sigma_b)}{Y | b ~ N(f(X^b) + X^f beta + Zb, sigma^2)
  b ~ N(0, Sigma_b) }
  \item Binary response coming soon.
}}
\usage{
mstan4bart(formula,
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
           thin = 1L,
           chains = 4L,
           cores = getOption("mc.cores", 1L),
           refresh = max(iter \%/\% 10L, 1L),
           offset_type = c("default", "fixef", "ranef", "bart", "parametric"),
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
  predictor. See \code{\link{lm}} and \code{\link{model.matrix.default}}.}
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
  running multi-threaded.}
  \item{iter}{positive integer indicating the number of posterior samples to
  draw and return. Includes warmup samples.}
  \item{warmup}{non-negative integer indicating number of posterior samples to
  draw and throw away. Also controls the number of iterations for which
  adaptation is run.}
  \item{thin}{one or two positive integer. Every \code{thin} sample will be
  kept, while every other sample will be discarded. If is length two, will
  attempt to use the named element \code{"bart"} for BART and \code{"stan"} for
  Stan. If not named, BART is the first thin element and Stan is the second.}
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
  \item{stan_args}{optional list, specifying further arguments to Stan. See
  details below.}
  \item{bart_args}{optional list, specifying further arguments to BART. See
  details below.}
}
\details{
  Fits a Bayesian \dQuote{mixed effect} model with a non-parametric Bayesian
  Additive Regression Trees (BART) component. For continuous responses:

  \itemize{
    \item \eqn{Y_i \mid b \sim \mathrm{N}\left(f(X^b_i) + X^f\beta +
    Z_ib_{g[i]}, \sigma^2\right)}{Y_i | b ~ N(f(X^b_i) + X^f\beta + Z_ib_{g[i]},
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
  
  \subsection{Terminology}{
    As \code{mstan4bart} fits a Bayesian model, essentially all components are
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
    \code{\link[rstanarm]{stan_glmer}} instead.
  }
  \subsection{Additional Arguments}{
    The \code{stan_args} and \code{bart_args} arguments to \code{mstan4bart} can
    be used to pass further arguments to \code{stan} and \code{bart}
    respectively. These are similar to the functions \code{\link[rstan]{stan}}
    and \code{\link{bart}}, but not identical as \code{mstan4bart} constructs its
    own model internally.
    
    Stan arguments include:
    \itemize{
      \item \code{prior_covariance} - see \code{\link[rstanarm]{stan_glmer}}
      \item \code{prior}, \code{prior_intercept}, \code{prior_aux} - see
      \code{\link[rstanarm]{stan_glm}}
      \item \code{init_r}, \code{adapt_gamma}, \code{adapt_delta},
      \code{adapt_kappa} - see \code{\link[rstan]{stan}}
    }
    
    BART arguments include:
    \itemize{
      \item further arguments to \code{\link{dbartsControl}} that are not
      specified by \code{mstan4bart}, such as \code{keepTrees} or
      \code{n.trees}; keeping trees can be costly in terms of memory, but is
      required to use \code{\link[=mstan4bart-generics]{predict}}
    }
  }
}
\value{
  Returns a list assigned class \code{mstan4bartFit}. Has components below, some
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
  \item{\code{bartData}}{\code{\link[=dbartsData]{data}} object used for BART
  component}
  \item{\code{X}}{fixed effect design matrix or NULL}
  \item{\code{X_means}}{column means of fixed effect design matrix when
  appropriate}
  \item{\code{reTrms}}{random effect \dQuote{terms} object when applicable, as
  used by \code{\link[lme4]{lmer}}}
  \item{\code{test}}{named list when applicable, having components \code{X} and
  \code{reTrms}; test data for BART is added to the \code{bartData} result}
  \item{\code{treatment}}{treatment vector, when applicable}
  
  Results, better accessed using \code{\link[=mstan4bart-generics]{extract}}:
  \item{\code{bart_train}}{samples of individual posterior predictions for BART
  component}
  \item{\code{bart_test}}{predicted test values for BART component, when
  applicable}
  \item{\code{bart_varcount}}{BART variable counts}
  \item{\code{sigma}}{samples of residual standard error}
  \item{\code{ranef}}{samples of random effects, or modeled coefficients; will
  be a named list, with effects for each grouping factor}
  \item{\code{Sigma}}{samples of covariance of random effects; also a named list
  with one element for each grouping factor}
  \item{\code{fixef}}{samples of the fixed effects, or unmodeled coefficients}
  
  Other items:
  \item{\code{warmup}}{a list of warmup samples, containing the same objects in
  the results subsection}
  \item{\code{sampler.bart}}{external points to BART samplers; used only for
  \code{\link[=mstan4bart-generics]{predict}} when \code{keepTrees} is \code{TRUE}}
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link{bart}}, \code{\link[rstan]{stan}}, \code{\link[lme4]{lmer}}
}
\examples{\donttest{
# Todo
}}
