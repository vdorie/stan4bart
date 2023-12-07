#include "stan_sampler.hpp"

#include <misc/stddef.h> // size_t
#include <cstring>       // memcpy, sterror

#include <misc/alloca.h> // misc_stackAllocate
#include <misc/string.h> // misc_str_matchAllInArray

#include <rc/util.h>
#include <rc/bounds.h>

#include <rstan/io/r_ostream.hpp>

#include <stan/services/util/create_unit_e_diag_inv_metric.hpp>

#include "stan_files/continuous.hpp" // should include Eigen

#include "double_writer.hpp"

#define TREE_DEPTH_OFFSET 1
#define DIVERGENT_TRANSITION_OFFSET 3

#if __cplusplus < 201112L
#  if defined(_WIN64) || SIZEOF_SIZE_T == 8
#    define SIZE_T_SPECIFIER "%lu"
#  else
#    define SIZE_T_SPECIFIER "%u"
#  endif
#else
#  define SIZE_T_SPECIFIER "%zu"
#endif

namespace {

std::vector<int> getIntVector(SEXP x) {
  return std::vector<int>(INTEGER(x), INTEGER(x) + rc_getLength(x));
}
std::vector<double> getDoubleVector(SEXP x) {
  return std::vector<double>(REAL(x), REAL(x) + rc_getLength(x));
}

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_d;

vector_d getEigenVector(SEXP x) {
  size_t len = rc_getLength(x);
  vector_d result(len);
  for (size_t i = 0; i < len; ++i)
    result(i) = REAL(x)[i];
  
  return result;
}

matrix_d getEigenMatrix(SEXP x) {
  int* dims = INTEGER(rc_getDims(x));
  matrix_d result(dims[0], dims[1]);
  size_t offset = 0;
  for (int col = 0; col < dims[1]; ++col) {
    for (int row = 0; row < dims[0]; ++row) {
      result(row, col) = REAL(x)[offset++];
    }
  }
  
  return result;
}

const char* const dataNames[] = {
  "N", "K", "X", "len_y", "lb_y", "ub_y", "y",
  "has_intercept", "is_binary",
  "prior_dist", "prior_dist_for_intercept", "prior_dist_for_aux",
  "has_weights", "weights", "offset_",
  "prior_scale", "prior_scale_for_intercept", "prior_scale_for_aux",
  "prior_mean", "prior_mean_for_intercept", "prior_mean_for_aux",
  "prior_df", "prior_df_for_intercept", "prior_df_for_aux",
  "global_prior_df", "global_prior_scale", "slab_df", "slab_scale", "num_normals",
  "t", "p", "l", "q", "len_theta_L",
  "shape", "scale", "len_concentration", "concentration",
  "len_regularization", "regularization",
  "num_non_zero", "w", "v", "u"
};

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

std::ostream nullout(nullptr);

} // anonymous namespace

namespace stan4bart {

#ifdef __clang__
#  if __has_warning("-Wenum-enum-conversion")
#    define SUPPRESS_ENUM_CONVERSION_WARNING 1
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wenum-enum-conversion"
#  endif
#endif

StanModel* createStanModelFromExpression(SEXP dataExpr)
{
  SEXP inputDataNamesExpr = rc_getNames(dataExpr);
  if (Rf_isNull(inputDataNamesExpr))
    Rf_error("names for dataExpr object cannot be NULL");
  
  size_t numDataNames = sizeof(dataNames) / sizeof(dataNames[0]);
  size_t* matchPos = misc_stackAllocate(numDataNames, size_t);
  
  size_t numInputDataNames = rc_getLength(inputDataNamesExpr);
  const char** inputDataNames = misc_stackAllocate(numInputDataNames, const char*);
  for (size_t i = 0; i < numInputDataNames; ++i)
    inputDataNames[i] = CHAR(STRING_ELT(inputDataNamesExpr, i));
  
  int errorCode = misc_str_matchAllInArray(dataNames, numDataNames, inputDataNames, numInputDataNames, matchPos);
  misc_stackFree(inputDataNames);
  
  if (errorCode != 0) {
    misc_stackFree(matchPos);
    Rf_error("error matching stan_data names: %s", std::strerror(errorCode));
  }
  
  for (size_t i = 0; i < numDataNames; ++i) {
    if (matchPos[i] == MISC_STR_NO_MATCH) {
      misc_stackFree(matchPos);
      Rf_error("mismatched stan data name: '%s' missing", dataNames[i]);
    }
  }
  
  // transformed data
  int prior_dist = rc_getIntAt(dataExpr, matchPos[9], "prior_dist", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 7, RC_END);
  int hs = 0;
  if (prior_dist <= 2) hs = 0;
  else if (prior_dist == 3) hs = 2;
  else if (prior_dist == 4) hs = 4;
  
  int t = rc_getIntAt(dataExpr, matchPos[29], "t", RC_VALUE | RC_GEQ, 0, RC_END);
  
  SEXP pExpr = VECTOR_ELT(dataExpr, matchPos[30]);
  if (t != static_cast<int>(rc_getLength(pExpr))) {
    misc_stackFree(matchPos);
    Rf_error("length of p must be equal to %d", t);
  }
  const int* p_int = INTEGER(pExpr);
  
  int len_var_group = 0;
  for (size_t i = 0; i < rc_getLength(pExpr); ++i) {
    if (p_int[i] == R_NaInt || p_int[i] < 1) {
      misc_stackFree(matchPos);
      Rf_error("p[" SIZE_T_SPECIFIER "] NA or less than 1", i + 1);
    }
    len_var_group += p_int[i];
  }
  int len_rho = len_var_group - t;
  
  int len_z_T = 0;
  int pos = 0;
  int len_concentration = rc_getIntAt(dataExpr, matchPos[36], "len_concentration", RC_VALUE | RC_GEQ, 0, RC_END);
  std::vector<double> delta_v(len_concentration);
  SEXP concentrationExpr = VECTOR_ELT(dataExpr, matchPos[37]);
  if (len_concentration > 0)
    rc_assertDoubleConstraints(concentrationExpr, "concentration", RC_VALUE | RC_GEQ, 0.0, RC_END);
  double* concentration = REAL(concentrationExpr);
  for (int i = 0; i < t; ++i) {
    if (p_int[i] > 1) {
      for (int j = 0; j < p_int[i]; ++j)
        delta_v[pos++] = concentration[j];
    }
    for (int j = 2; j < p_int[i]; ++j)
      len_z_T += p_int[i] - 1;
  }
  
  const int* l = INTEGER(VECTOR_ELT(dataExpr, matchPos[31]));
  for (int i = 0; i < t; ++i) {
    if (l[i] < 1) {
      misc_stackFree(matchPos);
      Rf_error("l[%d] less than 1", i + 1);
    }
  }
  
  /* Rprintf("name at pos: %s\n", inputDataNames[matchPos[11]]);
  Rprintf("prior mean for aux: %.16f\n", REAL(VECTOR_ELT(dataExpr, matchPos[11]))[0]);
  Rprintf("value >= 0.0: %s, >= -0.0: %s, >= -1.0e-16: %s\n", REAL(VECTOR_ELT(dataExpr, matchPos[11]))[0] >= 0.0 ? "true" : "false",
                                                    REAL(VECTOR_ELT(dataExpr, matchPos[11]))[0] >= -0.0 ? "true" : "false",
                                                    REAL(VECTOR_ELT(dataExpr, matchPos[11]))[0] >= -1.0e-16 ? "true" : "false"); */
  continuous_model_namespace::continuous_model* result = new continuous_model_namespace::continuous_model(
    rc_getIntAt(  dataExpr, matchPos[ 0], "N",     RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getIntAt(  dataExpr, matchPos[ 1], "K",     RC_VALUE | RC_GEQ, 0, RC_END),
    getEigenMatrix(VECTOR_ELT(dataExpr, matchPos[ 2])), // X
    rc_getIntAt(  dataExpr, matchPos[ 3], "len_y", RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getDouble0(VECTOR_ELT(dataExpr, matchPos[ 4]), "lb_y"),
    rc_getDouble0(VECTOR_ELT(dataExpr, matchPos[ 5]), "ub_y"),
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[ 6])), // y
    rc_getIntAt(  dataExpr, matchPos[ 7], "has_intercept", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 1, RC_END),
    rc_getIntAt(  dataExpr, matchPos[ 8], "is_binary", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 1, RC_END),
    prior_dist, // 9
    rc_getIntAt(  dataExpr, matchPos[10], "prior_dist_for_intercept", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 2, RC_END),
    rc_getIntAt(  dataExpr, matchPos[11], "prior_dist_for_aux", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 3, RC_END),
    rc_getIntAt(  dataExpr, matchPos[12], "has_weights",        RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 1, RC_END),
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[13])), // weights
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[14])), // offset_
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[15])), // prior_scale
    rc_getDoubleAt(dataExpr, matchPos[16], "prior_scale_for_intercept", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[17], "prior_scale_for_aux",       RC_VALUE | RC_GEQ, 0.0, RC_END),
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[18])), // prior_mean
    rc_getDoubleAt(dataExpr, matchPos[19], "prior_mean_for_intercept", RC_NA | RC_NO, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[20], "prior_mean_for_aux",  RC_VALUE | RC_GEQ, 0.0, RC_END),
    getEigenVector(VECTOR_ELT(dataExpr, matchPos[21])), // prior_df
    rc_getDoubleAt(dataExpr, matchPos[22], "prior_df_for_intercept", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[23], "prior_df_for_aux", RC_VALUE | RC_GEQ, 0.0, RC_END),
    
    rc_getDoubleAt(dataExpr, matchPos[24], "global_prior_df", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[25], "global_prior_scale", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[26], "slab_df", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[27], "slab_scale", RC_VALUE | RC_GEQ, 0.0, RC_END),
    getIntVector(VECTOR_ELT(dataExpr, matchPos[28])), // num_normals
    t, // 29
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[30])), // p
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[31])), // l
    rc_getIntAt(dataExpr, matchPos[32], "q",           RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getIntAt(dataExpr, matchPos[33], "len_theta_L", RC_VALUE | RC_GEQ, 0, RC_END),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[34])), // shape
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[35])), // scale
    len_concentration, // 36
    getDoubleVector(VECTOR_ELT(dataExpr, matchPos[37])), // concentration
    rc_getIntAt(dataExpr, matchPos[38], "len_regularization", RC_VALUE | RC_GEQ, 0, RC_END),
    getDoubleVector(VECTOR_ELT(dataExpr, matchPos[39])), // regularization
    rc_getIntAt(dataExpr, matchPos[40], "num_non_zero", RC_VALUE | RC_GEQ, 0, RC_END),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[41])), // w
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[42])), // v
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[43])), // u
    len_z_T,
    len_var_group,
    len_rho,
    pos,
    delta_v,
    hs
  );
  
  // check some extra constraints
  size_t n = static_cast<size_t>(INTEGER(VECTOR_ELT(dataExpr, matchPos[0]))[0]);
  size_t k = static_cast<size_t>(INTEGER(VECTOR_ELT(dataExpr, matchPos[1]))[0]);
  
  const int* x_dims = Rf_isNull(rc_getDims(VECTOR_ELT(dataExpr, matchPos[2]))) ? NULL : INTEGER(rc_getDims(VECTOR_ELT(dataExpr, matchPos[2])));
  
  double lb = REAL(VECTOR_ELT(dataExpr, matchPos[4]))[0];
  double ub = REAL(VECTOR_ELT(dataExpr, matchPos[5]))[0];
  const double* y = REAL(VECTOR_ELT(dataExpr, matchPos[6]));
  size_t y_len = rc_getLength(VECTOR_ELT(dataExpr, matchPos[6]));
  
  bool has_weights = rc_getIntAt(  dataExpr, matchPos[12], "has_weights",        RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 1, RC_END);
  size_t weights_len = rc_getLength(VECTOR_ELT(dataExpr, matchPos[13]));
  
  size_t offset_len = rc_getLength(VECTOR_ELT(dataExpr, matchPos[14]));
  
  SEXP prior_scaleExpr = VECTOR_ELT(dataExpr, matchPos[15]);
  SEXP prior_meanExpr  = VECTOR_ELT(dataExpr, matchPos[18]);
  SEXP prior_dfExpr    = VECTOR_ELT(dataExpr, matchPos[21]);
  
  SEXP num_normalsExpr = VECTOR_ELT(dataExpr, matchPos[28]);

  
  int numNonZero = rc_getInt0(VECTOR_ELT(dataExpr, matchPos[40]), "num_non_zero");
  const int* v = INTEGER(VECTOR_ELT(dataExpr, matchPos[42]));
  SEXP uExpr = VECTOR_ELT(dataExpr, matchPos[43]);
  const int* u = INTEGER(uExpr);
  
  int maxU = static_cast<int>(rc_getLength(VECTOR_ELT(dataExpr, matchPos[41]))) + 1;
  int q = INTEGER(VECTOR_ELT(dataExpr, matchPos[32]))[0];
  
  misc_stackFree(matchPos);
  
  if (x_dims == NULL) {
    delete result;
    Rf_error("x is not a matrix");
  }
  if (static_cast<size_t>(x_dims[0]) != n || static_cast<size_t>(x_dims[1]) != k) {
    delete result;
    Rf_error("x dim mismatch: got [%d, %d], expected [" SIZE_T_SPECIFIER ", " SIZE_T_SPECIFIER "]", x_dims[0], x_dims[1], n, k);
  }
  
  if (y_len != n) {
    delete result;
    Rf_error("y length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, y_len, n);
  }

  for (size_t i = 0; i < n; ++i) {
    if (y[i] < lb || y[i] > ub) {
      delete result;
      Rf_error("y[" SIZE_T_SPECIFIER "] out of [%f, %f] range", i + 1, lb, ub);
    }
  }
  if (rc_getLength(uExpr) != n + 1) {
    delete result;
    Rf_error("length of u must be equal to number of observations + 1");
  }
  for (size_t i = 0; i < n + 1; ++i) {
    if (u[i] < 0 || u[i] > maxU) {
      delete result;
      Rf_error("u[" SIZE_T_SPECIFIER "] out of [0, %d] range", i + 1, maxU);
    }
  }
  
  if (has_weights ? (weights_len != n) : (has_weights != 0)) {
    delete result;
    Rf_error("weights length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, weights_len, has_weights ? n : static_cast<size_t>(0));
  }
  
  if (offset_len != n) {
    delete result;
    Rf_error("offset length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, offset_len, n);
  }
  
  if (rc_getLength(prior_scaleExpr) != k) {
    delete result;
    Rf_error("prior_scale length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, rc_getLength(prior_scaleExpr), k);
  }
  for (size_t i = 0; i < k; ++i) {
    if (REAL(prior_scaleExpr)[i] < 0.0) {
      delete result;
      Rf_error("prior_scale[" SIZE_T_SPECIFIER "] out of [0, Inf) range", i + 1);
    }
  }
  if (rc_getLength(prior_meanExpr) != k) {
    delete result;
    Rf_error("prior_mean length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, rc_getLength(prior_meanExpr), k);
  }
  for (size_t i = 0; i < k; ++i) {
    if (REAL(prior_meanExpr)[i] < 0.0) {
      delete result;
      Rf_error("prior_mean[" SIZE_T_SPECIFIER "] out of [0, Inf) range", i + 1);
    }
  }
  if (rc_getLength(prior_dfExpr) != k) {
    delete result;
    Rf_error("prior_df length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, rc_getLength(prior_dfExpr), k);
  }
  for (size_t i = 0; i < k; ++i) {
    if (REAL(prior_dfExpr)[i] < 0.0) {
      delete result;
      Rf_error("prior_df[" SIZE_T_SPECIFIER "] out of [0, Inf) range", i + 1);
    }
  }
  
  if (prior_dist == 7) {
    if (rc_getLength(num_normalsExpr) != k) {
      delete result;
      Rf_error("num_normals length mismatch: got " SIZE_T_SPECIFIER ", expected " SIZE_T_SPECIFIER, rc_getLength(num_normalsExpr), k);
    }
    for (size_t i = 0; i < k; ++i) {
      if (INTEGER(num_normalsExpr)[i] < 2) {
        delete result;
        Rf_error("num_normals[" SIZE_T_SPECIFIER "] out of [2, Inf) range", i);
      }
    }
  } else if (rc_getLength(num_normalsExpr) != 0) {
    delete result;
    Rf_error("prior_df length mismatch: got " SIZE_T_SPECIFIER ", expected 0", rc_getLength(num_normalsExpr));
  }
  
  for (size_t i = 0; i < static_cast<size_t>(numNonZero); ++i) {
    if (v[i] < 0 || v[i] > q - 1) {
      delete result;
      Rf_error("v[" SIZE_T_SPECIFIER "] out of [0, %d] range", i + 1, q - 1);
    }
  }
  
  return result;
}

#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic pop
#endif

void deleteStanModel(StanModel* model) {
  delete model;
}

#ifdef SUPPRESS_ENUM_CONVERSION_WARNING
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wenum-enum-conversion"
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
#endif

StanSampler::StanSampler(StanModel& stanModel, const StanControl& stanControl, int chain_id, int num_warmup, int verbose) :
  c_out(verbose >  0 ? &rstan::io::rcout : &nullout),
  c_err(verbose >= 0 ? &rstan::io::rcerr : &nullout),
  logger(*c_out, *c_out, *c_out, *c_err, *c_err),
  diagnostic_writer(diagnostic_stream, "# "),
  init_context_ptr(new stan::io::empty_var_context()),
  init_writer("init"),
  sample_writer("sample"),
  dmp(stan::services::util::create_unit_e_diag_inv_metric(stanModel.num_params_r())),
  unit_e_metric(dmp),
  sampler(NULL)
{
  stanModel.constrained_param_names(constrained_param_names);
  stan::mcmc::sample::get_sample_param_names(sample_names);
  // start NUTS block
  sampler_names.resize(5);
  sampler_names[0] = "stepsize__";
  sampler_names[TREE_DEPTH_OFFSET] = "treedepth__";
  sampler_names[2] = "n_leapfrog__";
  sampler_names[DIVERGENT_TRANSITION_OFFSET] = "divergent__";
  sampler_names[4] = "energy__";
  sample_writer_offset = sample_names.size() + sampler_names.size();
  
  init_writer.resize(stanModel.num_params_r(), 1);
  
  num_pars = static_cast<int>(sample_names.size() + sampler_names.size() + constrained_param_names.size());
  // sample_writer.resize(sample_names.size() + sampler_names.size() + constrained_param_names.size(), num_warmup + num_iter);
  
  try {
  sampler = new stan4bart::interruptable_sampler<StanModel>(
    stanModel,
    *init_context_ptr,
    unit_e_metric,
    stanControl.random_seed,
    chain_id,
    stanControl.init_radius,
    num_warmup,
    stanControl.skip,
    stanControl.stepsize,
    stanControl.stepsize_jitter,
    stanControl.max_treedepth,
    stanControl.adapt_delta,
    stanControl.adapt_gamma,
    stanControl.adapt_kappa,
    stanControl.adapt_t0,
    stanControl.adapt_init_buffer,
    stanControl.adapt_term_buffer,
    stanControl.adapt_window,
    interrupt, 
    logger,
    init_writer,
    sample_writer,
    diagnostic_writer);
  } catch (const std::bad_alloc& e) {
    Rprintf("bad alloc: %s", e.what());
    throw e;
  }
}

void StanSampler::setVerbose(int level)
{
  c_out = level >  0 ? &rstan::io::rcout : &nullout;
  c_err = level >= 0 ? &rstan::io::rcerr : &nullout;
  new (&logger) stan::callbacks::stream_logger(*c_out, *c_out, *c_out, *c_err, *c_err);
}

void setStanOffset(StanModel& model, const double* offset)
{
  model.set_offset(offset);
}

void setResponse(StanModel& model, const double* response)
{
  model.set_response(response);
}

void StanSampler::getParametricMean(const StanModel& model, double* result) const
{
  model.get_parametric_mean(sample_writer.x_curr + sample_writer_offset, result);
}

void StanSampler::getParametricMean(const StanModel& model, double* result,
                                    bool includeFixed, bool includeRandom) const
{
  model.get_parametric_mean(sample_writer.x_curr + sample_writer_offset, result,
                            includeFixed, includeRandom);
}

double StanSampler::getSigma(const StanModel& model) const
{
  return model.get_aux(sample_writer.x_curr + sample_writer_offset);
}

void StanSampler::copyOutParameters(double* result, int offset) const
{
  std::memcpy(result, const_cast<const double*>(sample_writer.x_curr) + offset * num_pars, num_pars * sizeof(double));
}

void StanSampler::run(bool isWarmup)
{
  sampler->run(isWarmup);
}

/* 
bool StanSampler::isDivergentTransition() const {
  return sample_writer.x_curr[DIVERGENT_TRANSITION_OFFSET] == 1.0;
}

int StanSampler::getTreeDepth() const {
  return static_cast<int>(sample_writer.x_curr[TREE_DEPTH_OFFSET] + 0.5);
}
*/

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


void printStanModel(const StanModel* model)
{
  model->print(rstan::io::rcout);
}


}

