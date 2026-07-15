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

namespace {

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_d;

SEXP getElt(SEXP list, const char* name) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_isNull(names)) Rf_error("stan data list has no names");
  R_xlen_t n = XLENGTH(list);
  for (R_xlen_t i = 0; i < n; ++i)
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(list, i);
  Rf_error("stan data field '%s' missing", name);
  return R_NilValue;  // unreachable
}

int getIntElt(SEXP list, const char* name) {
  SEXP e = getElt(list, name);
  if (TYPEOF(e) == INTSXP) return INTEGER(e)[0];
  if (TYPEOF(e) == REALSXP) return static_cast<int>(REAL(e)[0]);
  Rf_error("stan data field '%s' not scalar-numeric", name);
  return 0;
}

double getDoubleElt(SEXP list, const char* name) {
  SEXP e = getElt(list, name);
  if (TYPEOF(e) == REALSXP) return REAL(e)[0];
  if (TYPEOF(e) == INTSXP) return static_cast<double>(INTEGER(e)[0]);
  Rf_error("stan data field '%s' not scalar-numeric", name);
  return 0.0;
}

vector_d getEigenVec(SEXP list, const char* name) {
  SEXP e = getElt(list, name);
  R_xlen_t len = XLENGTH(e);
  vector_d out(len);
  if (TYPEOF(e) == REALSXP)      for (R_xlen_t i = 0; i < len; ++i) out(i) = REAL(e)[i];
  else if (TYPEOF(e) == INTSXP)  for (R_xlen_t i = 0; i < len; ++i) out(i) = INTEGER(e)[i];
  else Rf_error("stan data field '%s' not numeric", name);
  return out;
}

std::vector<int> getIntVec(SEXP list, const char* name) {
  SEXP e = getElt(list, name);
  R_xlen_t len = XLENGTH(e);
  std::vector<int> out(len);
  if (TYPEOF(e) == INTSXP)       for (R_xlen_t i = 0; i < len; ++i) out[i] = INTEGER(e)[i];
  else if (TYPEOF(e) == REALSXP) for (R_xlen_t i = 0; i < len; ++i) out[i] = static_cast<int>(REAL(e)[i]);
  else Rf_error("stan data field '%s' not integer", name);
  return out;
}

matrix_d getEigenMat(SEXP list, const char* name) {
  SEXP e = getElt(list, name);
  SEXP dims = Rf_getAttrib(e, R_DimSymbol);
  if (Rf_isNull(dims) || XLENGTH(dims) != 2) Rf_error("stan data field '%s' not a matrix", name);
  int nr = INTEGER(dims)[0], nc = INTEGER(dims)[1];
  matrix_d out(nr, nc);
  R_xlen_t off = 0;
  for (int c = 0; c < nc; ++c) for (int r = 0; r < nr; ++r) out(r, c) = REAL(e)[off++];
  return out;
}

// Build a ParametricModel from the marshaled Stan data list, replicating the
// transformed-data derivations (len_z_T, len_rho, delta) from init.cpp.
stan4bart::ParametricModel buildModel(SEXP dataExpr) {
  stan4bart::ParametricModel pm;

  pm.N = getIntElt(dataExpr, "N");
  pm.K = getIntElt(dataExpr, "K");
  pm.t = getIntElt(dataExpr, "t");
  pm.q = getIntElt(dataExpr, "q");
  pm.is_binary   = getIntElt(dataExpr, "is_binary") != 0;
  pm.has_weights = getIntElt(dataExpr, "has_weights") != 0;

  pm.X = getEigenMat(dataExpr, "X");
  pm.y_ = getEigenVec(dataExpr, "y");
  pm.offset_ = getEigenVec(dataExpr, "offset_");
  if (pm.has_weights) pm.weights = getEigenVec(dataExpr, "weights");

  pm.prior_dist  = getIntElt(dataExpr, "prior_dist");
  pm.prior_scale = getEigenVec(dataExpr, "prior_scale");
  pm.prior_mean  = getEigenVec(dataExpr, "prior_mean");
  pm.prior_df    = getEigenVec(dataExpr, "prior_df");

  if (!pm.is_binary) {
    pm.prior_dist_for_aux  = getIntElt(dataExpr, "prior_dist_for_aux");
    pm.prior_scale_for_aux = getDoubleElt(dataExpr, "prior_scale_for_aux");
    pm.prior_mean_for_aux  = getDoubleElt(dataExpr, "prior_mean_for_aux");
    pm.prior_df_for_aux    = getDoubleElt(dataExpr, "prior_df_for_aux");
  }

  pm.p = getIntVec(dataExpr, "p");
  pm.l = getIntVec(dataExpr, "l");
  pm.len_theta_L       = getIntElt(dataExpr, "len_theta_L");
  pm.len_concentration = getIntElt(dataExpr, "len_concentration");
  if (pm.t > 0) {
    pm.re_scale  = getEigenVec(dataExpr, "scale");
    pm.tau_shape = getEigenVec(dataExpr, "shape");
    vector_d reg = getEigenVec(dataExpr, "regularization");
    pm.regularization = std::vector<double>(reg.data(), reg.data() + reg.size());
  }

  // transformed-data derivations (init.cpp:157-182)
  int len_var_group = 0;
  for (int i = 0; i < pm.t; ++i) len_var_group += pm.p[i];
  pm.len_rho = len_var_group - pm.t;

  vector_d concentration = pm.len_concentration > 0 ? getEigenVec(dataExpr, "concentration") : vector_d();
  pm.len_z_T = 0;
  pm.delta.clear();
  for (int i = 0; i < pm.t; ++i) {
    if (pm.p[i] > 1)
      for (int j = 0; j < pm.p[i]; ++j) pm.delta.push_back(concentration[j]);
    for (int j = 2; j < pm.p[i]; ++j) pm.len_z_T += pm.p[i] - 1;
  }

  pm.Zw = getEigenVec(dataExpr, "w");
  pm.Zv = getIntVec(dataExpr, "v");
  pm.Zu = getIntVec(dataExpr, "u");

  pm.finalize();
  return pm;
}

}  // namespace

extern "C" {

SEXP stan4bart_logdensity_grad(SEXP dataExpr, SEXP parExpr) {
  stan4bart::ParametricModel pm = buildModel(dataExpr);

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
