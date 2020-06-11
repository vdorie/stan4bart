#ifndef STAN_UTIL_HPP
#define STAN_UTIL_HPP

#include <ext/Rinternals.h>

namespace stan4bart {

  struct double_writer;
  
  SEXP createStanResultsExpr(const double_writer& sample_writer);

}

#endif // BART_UTIL_HPP
