\name{NEWS}
\title{stan4bart News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 0.0-14}{
  \subsection{UPGRADING}{
    \itemize{
      \item Requires \pkg{dbarts} (>= 1.0-0) and a C++20 compiler. Results
            differ from earlier versions, including under a fixed seed.
      \item Stan is no longer used. The parametric component is drawn with the
            WALNUTS sampler, so the name \code{stan_args} is kept but
            configures that sampler. The NUTS tuning arguments
            (\code{adapt_gamma}, \code{adapt_kappa}, \code{adapt_t0},
            \code{adapt_init_buffer}, \code{adapt_term_buffer},
            \code{adapt_window}, \code{stepsize}, \code{stepsize_jitter},
            \code{max_treedepth}) are ignored with a warning. The shrinkage
            priors \code{hs}, \code{hs_plus}, \code{lasso}, \code{laplace} and
            \code{product_normal} are no longer supported and \code{prior_intercept}
            is dropped with a warning.
      \item \code{bart_args} follows \pkg{dbarts} 1.0-0: \code{node.prior} is
            now \code{leaf.prior} (an error naming the new spelling),
            \code{rngSeed} is now \code{seed} (a warning), and the residual
            model arguments (\code{sigma}, \code{resid.prior},
            \code{resid.dist}) are refused. A name in \code{bart_args} that
            nothing recognizes is now an error; before it was silently ignored.
      \item Warmup draws are no longer stored by default. A per-chain summary
            of the sampler tuning is returned as \code{adaptation} instead;
            use \code{save_warmup = TRUE} to keep the draws.
      \item For a binary response, \code{weights} that are all 0 or 1 now select
            the rows in use and any other weights are an error.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{mvbart} fits multivariate BART to correlated continuous
            outcomes with a shared predictor set.
      \item \code{store = "trees"} keeps the trees rather than the full matrix
            of BART fits, which shrinks large fit objects considerably;
            \code{extract}, \code{fitted} and \code{predict} work either way.
      \item \code{print}, \code{summary} and \code{residuals} methods for
            fitted models; the first two warn when the sampler looks
            under-warmed.
      \item \code{skip} accepts a \code{"stan"} element to take several
            parametric transitions per BART sweep, which improves mixing of the
            random effect standard deviations.
      \item Factors in the \code{bart} component are split by level subsets
            rather than indicator columns, so \code{extract(type = "varcount")}
            reports one row per factor.
    }
  }
}
\section{CHANGES IN VERSION 0.0-8}{
  \subsection{BUG FIXES}{
    \itemize{
      \item No longer uses \code{rpath} on Window64 arm. PR thanks to
            Andrew Johnson.
    }
  }
}
\section{CHANGES IN VERSION 0.0-7}{
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{extract} argument now works correctly with varcount. Bug
            report thanks to Joshua Bon.
      \item \code{predict} now works with single chains. Bug report thanks to
            github user Pentaonia (Louis).
    }
  }
}
\section{CHANGES IN VERSION 0.0-4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item An experimental implementation of callbacks has been added as
            has the option to avoid writing out the full set of Stan and
            BART samples. This allows a model to be fit just saving samples
            of the quantity of interested.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{treated} argument now works correctly with logical vectors.
    }
  }
}
\section{CHANGES IN VERSION 0.0-3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{extract} now accepts as a type \code{"trees"}, which
            allows for easier inspection of models fit with
            \code{bart_args = list("keepTrees" = TRUE)}.
    }
  }
}
\section{CHANGES IN VERSION 0.0-2}{
  \itemize{
    \item Initial public release.
  }
}

