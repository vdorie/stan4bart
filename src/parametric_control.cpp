/// \file parametric_control.cpp
/// \brief De-Stanned R-marshaling seam for the parametric conditional's
///        sampler: parse the sampler control list, print it, and re-emit the
///        stored draws as the R `stan` results object.
///
/// These three routines used to live in stan_sampler.cpp alongside the Stan
/// model marshaling; they carry no Stan dependency (only rc/misc + the shared
/// double_writer), so they survive the Stan deletion here. The result object's
/// shape and row names are UNCHANGED - the R generics consume it by name.

#include "parametric_sampler.hpp"  // ParametricControl, double_writer, decls

#include <ext/Rinternals.h>
#include <R_ext/Print.h>  // Rprintf

#include <misc/stddef.h>  // size_t
#include <cstring>        // memcpy, strerror

#include <misc/alloca.h>  // misc_stackAllocate
#include <misc/string.h>  // misc_str_matchAllInArray

#include <rc/util.h>
#include <rc/bounds.h>

namespace {

// The control vocabulary the R side passes. The NUTS-specific knobs
// (adapt_*, stepsize*, max_treedepth) have no faithful WALNUTS analog and are
// accepted-but-ignored (parsed, never consumed); only seed, init_r, and skip
// reach the WALNUTS sampler.
const char* const controlNames[] = {
  "seed",
  "init_r",
  "skip",
  "adapt_gamma",
  "adapt_delta",
  "adapt_kappa",
  "adapt_init_buffer",
  "adapt_term_buffer",
  "adapt_window",
  "adapt_t0",
  "stepsize",
  "stepsize_jitter",
  "max_treedepth"
};

}  // anonymous namespace

namespace stan4bart {

// The rc bounds vocabulary ORs distinct enum types (RC_VALUE | RC_GEQ, ...),
// which C++20 deprecates; suppress it exactly as the deleted stan_sampler.cpp
// did around this same code.
#ifdef __clang__
#  if __has_warning("-Wenum-enum-conversion")
#    define SUPPRESS_ENUM_CONVERSION_WARNING 1
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wenum-enum-conversion"
#  endif
#endif

void initializeStanControlFromExpression(StanControl& control, SEXP controlExpr)
{
  SEXP controlNamesExpr = rc_getNames(controlExpr);
  if (Rf_isNull(controlNamesExpr))
    Rf_error("names for stanControl object cannot be NULL");

  size_t numControlNames = sizeof(controlNames) / sizeof(controlNames[0]);
  size_t* matchPos = misc_stackAllocate(numControlNames, size_t);

  size_t numInputControlNames = rc_getLength(controlNamesExpr);
  const char** inputControlNames = misc_stackAllocate(numInputControlNames, const char*);
  for (size_t i = 0; i < numInputControlNames; ++i)
    inputControlNames[i] = CHAR(STRING_ELT(controlNamesExpr, i));

  int errorCode = misc_str_matchAllInArray(controlNames, numControlNames, inputControlNames, numInputControlNames, matchPos);
  misc_stackFree(inputControlNames);

  if (errorCode != 0) {
    misc_stackFree(matchPos);
    Rf_error("error matching names: %s", std::strerror(errorCode));
  }

  if (matchPos[0] == MISC_STR_NO_MATCH)
    Rf_error("stanControl requires 'seed' to be specified");

  control.random_seed = static_cast<unsigned int>(rc_getInt0(VECTOR_ELT(controlExpr, matchPos[0]), "seed"));
  control.init_radius = rc_getDoubleAt(controlExpr, matchPos[1], "init_r",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 2.0, RC_END);
  control.skip = rc_getIntAt(controlExpr, matchPos[2], "skip",
    RC_VALUE | RC_GT, 0,
    RC_NA | RC_YES, RC_END);
  control.adapt_gamma = rc_getDoubleAt(controlExpr, matchPos[3], "adapt_gamma",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 0.05, RC_END);
  control.adapt_delta = rc_getDoubleAt(controlExpr, matchPos[4], "adapt_delta",
    RC_VALUE | RC_GT, 0.0,
    RC_VALUE | RC_LT, 1.0,
    RC_VALUE | RC_DEFAULT, 0.8, RC_END);
  control.adapt_kappa = rc_getDoubleAt(controlExpr, matchPos[5], "adapt_kappa",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 0.75, RC_END);
  control.adapt_init_buffer = static_cast<unsigned int>(rc_getIntAt(controlExpr, matchPos[6], "adapt_init_buffer",
    RC_VALUE | RC_DEFAULT, 75u, RC_END));
  control.adapt_term_buffer = static_cast<unsigned int>(rc_getIntAt(controlExpr, matchPos[7], "adapt_term_buffer",
    RC_VALUE | RC_DEFAULT, 50u, RC_END));
  control.adapt_window = static_cast<unsigned int>(rc_getIntAt(controlExpr, matchPos[8], "adapt_window",
    RC_VALUE | RC_DEFAULT, 25u, RC_END));
  control.adapt_t0 = rc_getDoubleAt(controlExpr, matchPos[9], "adapt_t0",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 10.0, RC_END);
  control.stepsize = rc_getDoubleAt(controlExpr, matchPos[10], "stepsize",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 1.0, RC_END);
  control.stepsize_jitter = rc_getDoubleAt(controlExpr, matchPos[11], "stepsize_jitter",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_LEQ, 1.0,
    RC_VALUE | RC_DEFAULT, 0.0, RC_END);
  control.max_treedepth = rc_getIntAt(controlExpr, matchPos[12], "max_treedepth",
    RC_VALUE | RC_GEQ, 0,
    RC_VALUE | RC_DEFAULT, 10, RC_END);

  misc_stackFree(matchPos);
}

#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#  undef SUPPRESS_ENUM_CONVERSION_WARNING
#endif

SEXP createStanResultsExpr(const double_writer& sample_writer)
{
  SEXP resultExpr = PROTECT(rc_newReal(sample_writer.num_pars * sample_writer.num_samples));
  rc_setDims(resultExpr, static_cast<int>(sample_writer.num_pars), static_cast<int>(sample_writer.num_samples), -1);
  std::memcpy(REAL(resultExpr), sample_writer.x_base, sample_writer.num_pars * sample_writer.num_samples * sizeof(double));

  SEXP dimnamesExpr = PROTECT(rc_newList(2));
  SET_VECTOR_ELT(dimnamesExpr, 0, rc_newCharacter(sample_writer.num_pars));
  SET_VECTOR_ELT(dimnamesExpr, 1, R_NilValue);

  rc_setDimNames(resultExpr, dimnamesExpr);

  SEXP namesExpr = VECTOR_ELT(dimnamesExpr, 0);
  for (size_t i = 0; i < sample_writer.num_pars; ++i)
    SET_STRING_ELT(namesExpr, i, Rf_mkChar(sample_writer.names[i].c_str()));

  UNPROTECT(2);

  return resultExpr;
}

void printStanControl(const StanControl& control)
{
  Rprintf("  seed: %u\n"
          "  init_r: %f\n"
          "  skip: %d\n"
          "  adapt_gamma: %f\n"
          "  adapt_delta: %f\n"
          "  adapt_kappa: %f\n"
          "  adapt_init_buffer: %u\n"
          "  adapt_term_buffer: %u\n"
          "  adapt_window: %u\n"
          "  adapt_t0: %f\n"
          "  stepsize: %f\n"
          "  stepsize_jitter: %f\n"
          "  max_treedepth: %d\n",
          control.random_seed, control.init_radius,
          control.skip, control.adapt_gamma, control.adapt_delta,
          control.adapt_kappa, control.adapt_init_buffer, control.adapt_term_buffer,
          control.adapt_window, control.adapt_t0, control.stepsize,
          control.stepsize_jitter, control.max_treedepth);
}

}  // namespace stan4bart
