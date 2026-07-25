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
  giving the formula, family, chain/sample counts, and the per-chain warmup
  and sampling mean leapfrog counts.
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
