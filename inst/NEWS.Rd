\name{NEWS}
\title{stan4bart News}
\encoding{UTF-8}
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

