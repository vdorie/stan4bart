\name{NEWS}
\title{stan4bart News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 0.0-13.1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Compiles and fits under Eigen 5, which a coming release of
            \pkg{RcppEigen} brings; the package failed to compile there,
            and once compiled crashed when fitting. Report and fix thanks
            to Dirk Eddelbuettel.
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

