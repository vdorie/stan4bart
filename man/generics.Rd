\name{stan4bart-generics}

\alias{extract}
\alias{stan4bart-generics}
\alias{extract.stan4bartFit}
\alias{fitted.stan4bartFit}
\alias{predict.stan4bartFit}

\title{Generic Functions for stan4bart Model Fits}
\description{
  Commonly expected utility functions to derive useful quantities from fitted
  models.
}
\usage{
\method{extract}{stan4bartFit}(
  object,
  type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
           "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan",
           "trees", "callback"),
  sample = c("train", "test"),
  combine_chains = TRUE,
  sample_new_levels = TRUE,
  include_warmup = FALSE,
  \dots)
  
\method{fitted}{stan4bartFit}(
  object,
  type = c("ev", "ppd", "fixef", "indiv.fixef", "ranef", "indiv.ranef",
           "indiv.bart", "sigma", "Sigma", "k", "varcount", "stan",
           "callback"),
  sample = c("train", "test"),
  sample_new_levels = TRUE,
  \dots)

\method{predict}{stan4bartFit}(
  object, newdata, offset,
  type = c("ev", "ppd", "indiv.fixef", "indiv.ranef", "indiv.bart"),
  combine_chains = TRUE,
  sample_new_levels = TRUE,
  \dots)
}
\arguments{
  \item{object}{a fitted model resulting from a call to
  \code{\link{stan4bart}}.}
  \item{type}{a character vector; one of the options listed below.}
  \item{sample}{one of \code{"train"} or \code{"test"}, indicating if the
  training or test data frames should be used.}
  \item{combine_chains}{logical controlling if chain information should be
  discarded and the result returned as a matrix instead of an array.}
  \item{sample_new_levels}{logical; if \code{TRUE}, levels out of the training
  sample will have random effects drawn from their posterior predictive
  distribution. If \code{FALSE}, their random effects will be fixed to 0.}
  \item{include_warmup}{logical or \code{"only"}; when \code{TRUE}/\code{FALSE},
  warmup samples will or will not be included in the result respectively. When
  \code{"only"}, only the warmup samples will be returned.}
  \item{newdata}{data frame for making out of sample predictions.}
  \item{offset}{optional vector which will be added to test predictors.}
  \item{\dots}{not currently in use, but provided to match signatures of
  other generics.}
}
\details{
  \code{extract} is used to obtain raw samples using the training or test data,
  \code{fitted} averages those samples, and \code{predict} operates on data
  not available at the time of fitting. Note: \code{predict} requires that the
  model be fit with \code{args_bart = list(keepTrees = TRUE)}.
  
  \subsection{Return type}{
    The type argument accepts:
    \itemize{
      \item \code{"ev"} - the individual level expected value, that is draws
      from \eqn{E[Y \mid X^b, X^f, Z] \mid Y = f(X^b) + X^f\beta + Zb
      \mid Y}{E[Y | X^b, X^f, Z] | Y = f(X^b) + X^f \beta + Zb
      \mid Y} where the expectation is with respect to the posterior
      distribution of the parameters given the data
      \item \code{"ppd"} - draws from the individual level posterior predictive
      distribution, generally speaking adding noise to the result for
      \code{"ev"} or simulating new Bernoulli trials.
      \item \code{"fixef"} - draws from the posterior of the fixed effects
      (also known as the \dQuote{unmodeled} coefficients),
      \eqn{\beta \mid Y}{\beta | Y}
      \item \code{"indiv.fixef"} - draws from the posterior distribution of the
      individual level mean component deriving from the fixed effects,
      \eqn{X^f\beta}
      \item \code{"ranef"} - the random effects, varying intercepts and slopes,
      or \dQuote{modeled} coefficients, \eqn{b}; \eqn{b} has substantial
      structure that is represented as the returned value, where coefficients
      are reported within their grouping factors
      \item \code{"indiv.ranef"} - individual level mean component deriving
      from the random effects, \eqn{Zb}
      \item \code{"indiv.bart"} - individual level mean component deriving
      from the BART model, \eqn{f(X^b)}
      \item \code{"sigma"} - for continuous responses, the residual standard
      error
      \item \code{"Sigma"} - when applicable, the covariance matrices of the
      random effects
      \item \code{"stan"} - raw matrix or array of Stan sampled transformed
      parameters.
      \item \code{"trees"} - a data frame of flatted trees; see the subsection
      on extracted trees in \code{\link[dbarts]{bart}} and note that stan4bart variable
      names can be found in the \code{bartData@x} element of a fitted stan4bart
      model
      \item \code{"callback"} - if a callback function was provided while
      fitting, the results of that for each sample
    }
  }
}
\value{
  \code{extract} and \code{predict} return either arrays of dimensions equal to
  \code{n.observations x n.samples x n.chains} when \code{combine_chains} is
  \code{FALSE}, or matrices of dimensions equal to
  \code{n.observations x (n.samples * n.chains)} when \code{combine_chains} is
  \code{TRUE}.
  
  \code{fitted} returns a vector of the appropriate length by averaging the
  result of a call to \code{extract}.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}

