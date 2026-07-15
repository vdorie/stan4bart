/// \file logdensity_export.cpp
/// \brief Internal .Call gate surface for the WALNUTS parametric target (C1).
///
/// Two dot-prefixed entries, both taking the marshaled Stan data list (the
/// same named list continuous_model consumes) plus an unconstrained position:
///   .stan4bart_logdensity_grad(data, par) -> list(value, gradient)
///       evaluates the hand-written ParametricModel functor.
///   .stan4bart_stan_logdensity(data, par) -> numeric
///       evaluates Stan's own continuous_model::log_prob<false, true> at the
///       SAME point (the near-exact oracle for the target cross-check).
/// These are the correctness gates; they carry no man pages.

#include <ext/Rinternals.h>

#include <cstring>
#include <vector>

#include "stan_sampler.hpp"           // createStanModelFromExpression, StanModel
#include <stan/model/model_base.hpp>  // model_base::log_prob_jacobian (the vtable slot)

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

SEXP stan4bart_stan_logdensity(SEXP dataExpr, SEXP parExpr) {
  // Construct Stan's own continuous_model (heavy machinery lives in
  // stan_sampler.o) and evaluate its log_prob<false, true> via the base-class
  // virtual log_prob_jacobian. Going through the vtable avoids re-including the
  // generated continuous.hpp here, which would duplicate its namespace-scope
  // globals (e.g. profiles__) at link time. continuous_model derives from
  // stan::model::model_base through a single-inheritance chain, so its
  // model_base subobject is at offset 0 and the cast is address-preserving.
  stan4bart::StanModel* model = stan4bart::createStanModelFromExpression(dataExpr);
  stan::model::model_base* mb = reinterpret_cast<stan::model::model_base*>(model);

  double value = R_NaReal;
  try {
    R_xlen_t np = XLENGTH(parExpr);
    Eigen::VectorXd par(np);
    for (R_xlen_t i = 0; i < np; ++i) par(i) = REAL(parExpr)[i];
    value = mb->log_prob_jacobian(par, nullptr);
  } catch (const std::exception& e) {
    stan4bart::deleteStanModel(model);
    Rf_error("stan log_prob failed: %s", e.what());
  }
  stan4bart::deleteStanModel(model);
  return Rf_ScalarReal(value);
}

}  // extern "C"
