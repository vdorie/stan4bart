/// \file logdensity_export.cpp
/// \brief Internal .Call gate surface for the WALNUTS parametric target.
///
/// One dot-prefixed entry, taking the marshaled data list (the same named list
/// the sampler consumes) plus an unconstrained position:
///   .stan4bart_logdensity_grad(data, par) -> list(value, gradient)
///       evaluates the hand-written ParametricModel functor.
/// This is the finite-difference correctness gate; it carries no man page. The
/// former Stan log_prob cross-check oracle went with the deleted Stan model.

#include <ext/Rinternals.h>

#include <cstring>
#include <vector>

#include "parametric_model.hpp"
#include "parametric_model_io.hpp"    // buildParametricModel (shared with walnuts_sampler.cpp)

extern "C" {

SEXP stan4bart_logdensity_grad(SEXP dataExpr, SEXP parExpr) {
  stan4bart::ParametricModel pm = stan4bart::buildParametricModel(dataExpr);

  R_xlen_t np = XLENGTH(parExpr);
  if (static_cast<int>(np) != pm.dim())
    Rf_error("par length %d does not match model dimension %d", static_cast<int>(np), pm.dim());

  Eigen::VectorXd par(np);
  for (R_xlen_t i = 0; i < np; ++i) par(i) = REAL(parExpr)[i];

  double value = 0.0;
  Eigen::VectorXd grad(pm.dim());
  pm(par, value, grad);

  SEXP valueExpr = PROTECT(Rf_ScalarReal(value));
  SEXP gradExpr = PROTECT(Rf_allocVector(REALSXP, pm.dim()));
  for (int i = 0; i < pm.dim(); ++i) REAL(gradExpr)[i] = grad(i);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, valueExpr);
  SET_VECTOR_ELT(out, 1, gradExpr);
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("value"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("gradient"));
  Rf_setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(4);
  return out;
}

}  // extern "C"
