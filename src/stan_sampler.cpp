#include "stan_sampler.hpp"

#include <misc/stddef.h> // size_t
#include <cstring>       // memcpy, sterror

#include <misc/alloca.h> // misc_stackAllocate
#include <misc/string.h> // misc_str_matchAllInArray

#include <rc/util.h>
#include <rc/bounds.h>

#include <stan/services/util/create_unit_e_diag_inv_metric.hpp>

#if defined(__GNUC__) && (\
  (!defined(__clang__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6))) || \
  ( defined(__clang__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 7))))
#  define SUPPRESS_DIAGNOSTIC 1
#endif

#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wunused-local-typedef"
#endif
#include "stan_files/continuous.hpp" // should include Eigen
#ifdef SUPPRESS_DIAGNOSTIC
#  pragma GCC diagnostic pop
#endif

#include "double_writer.hpp"

namespace {

std::vector<int> getIntVector(SEXP x) {
  return std::vector<int>(INTEGER(x), INTEGER(x) + rc_getLength(x));
}
std::vector<double> getDoubleVector(SEXP x) {
  return std::vector<double>(REAL(x), REAL(x) + rc_getLength(x));
}

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;

vector_d getEigenVector(SEXP x) {
  size_t len = rc_getLength(x);
  vector_d result(len);
  for (size_t i = 0; i < len; ++i)
    result(i) = REAL(x)[i];
  
  return result;
}

const char* const dataNames[] = {
  "N", "len_y", "lb_y", "ub_y", "y",
  "prior_dist_for_aux", "has_weights", "weights",
  "offset_", "prior_scale_for_aux", "prior_mean_for_aux",
  "prior_df_for_aux", "t", "p", "l", "q", "len_theta_L",
  "shape", "scale", "len_concentration", "concentration",
  "len_regularization", "regularization", "num_non_zero",
  "w", "v", "u"
};

const char* const argNames[] = {
  "seed",
  "init_r",
  "thin",
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
  int t = rc_getIntAt(dataExpr, matchPos[12], "t", RC_VALUE | RC_GEQ, 0, RC_END);
  
  int len_var_group = 0;
  SEXP pExpr = VECTOR_ELT(dataExpr, matchPos[13]);
  rc_assertIntConstraints(pExpr, "p", RC_VALUE | RC_GEQ, 1, RC_END);
  const int* p_int = INTEGER(pExpr);
  for (size_t i = 0; i < rc_getLength(pExpr); ++i) {
    if (p_int[i] < 1) {
      misc_stackFree(matchPos);
      Rf_error("p[%lu] less than 1", i + 1);
    }
    len_var_group += p_int[i];
  }
  int len_rho = len_var_group - t;
  
  int len_z_T = 0;
  int pos = 0;
  int len_concentration = rc_getIntAt(dataExpr, matchPos[19], "len_concentration", RC_VALUE | RC_GEQ, 0, RC_END);
  std::vector<double> delta_v(len_concentration);
  SEXP concentrationExpr = VECTOR_ELT(dataExpr, matchPos[20]);
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
  
  const int* l = INTEGER(VECTOR_ELT(dataExpr, matchPos[14]));
  for (int i = 0; i < t; ++i) {
    if (l[i] < 1) {
      misc_stackFree(matchPos);
      Rf_error("l[%d] less than 1", i + 1);
    }
  }
  
  /* Rprintf("name at pos: %s\n", inputDataNames[matchPos[10]]);
  Rprintf("prior mean for aux: %.16f\n", REAL(VECTOR_ELT(dataExpr, matchPos[10]))[0]);
  Rprintf("value >= 0.0: %s, >= -0.0: %s, >= -1.0e-16: %s\n", REAL(VECTOR_ELT(dataExpr, matchPos[10]))[0] >= 0.0 ? "true" : "false",
                                                    REAL(VECTOR_ELT(dataExpr, matchPos[10]))[0] >= -0.0 ? "true" : "false",
                                                    REAL(VECTOR_ELT(dataExpr, matchPos[10]))[0] >= -1.0e-16 ? "true" : "false"); */
  model_continuous_namespace::model_continuous* result = new model_continuous_namespace::model_continuous(
    rc_getIntAt(  dataExpr, matchPos[ 0], "N",     RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getIntAt(  dataExpr, matchPos[ 1], "len_y", RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getDouble0(VECTOR_ELT(dataExpr, matchPos[ 2]), "lb_y"),
    rc_getDouble0(VECTOR_ELT(dataExpr, matchPos[ 3]), "ub_y"),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[ 4])), // y
    rc_getIntAt(  dataExpr, matchPos[ 5], "prior_dist_for_aux", RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 3, RC_END),
    rc_getIntAt(  dataExpr, matchPos[ 6], "has_weights",        RC_VALUE | RC_GEQ, 0, RC_VALUE | RC_LEQ, 1, RC_END),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[ 7])), // weights
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[ 8])), // offset_
    rc_getDoubleAt(dataExpr, matchPos[ 9], "prior_scale_for_aux", RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[10], "prior_mean_for_aux",  RC_VALUE | RC_GEQ, 0.0, RC_END),
    rc_getDoubleAt(dataExpr, matchPos[11], "prior_df_for_aux",    RC_VALUE | RC_GEQ, 0.0, RC_END),
    t,
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[13])), // p
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[14])), // l
    rc_getIntAt(dataExpr, matchPos[15], "q",           RC_VALUE | RC_GEQ, 0, RC_END),
    rc_getIntAt(dataExpr, matchPos[16], "len_theta_L", RC_VALUE | RC_GEQ, 0, RC_END),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[17])), // shape
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[18])), // scale
    len_concentration,
    getDoubleVector(VECTOR_ELT(dataExpr, matchPos[20])), // concentration
    rc_getIntAt(dataExpr, matchPos[21], "len_regularization", RC_VALUE | RC_GEQ, 0, RC_END),
    getDoubleVector(VECTOR_ELT(dataExpr, matchPos[22])), // regularization
    rc_getIntAt(dataExpr, matchPos[23], "num_non_zero", RC_VALUE | RC_GEQ, 0, RC_END),
    getEigenVector( VECTOR_ELT(dataExpr, matchPos[24])), // w
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[25])), // v
    getIntVector(   VECTOR_ELT(dataExpr, matchPos[26])), // u
    len_z_T,
    len_var_group,
    len_rho,
    pos,
    delta_v
  );
  
  // check some extra constraints
  size_t n = static_cast<size_t>(INTEGER(VECTOR_ELT(dataExpr, matchPos[0]))[0]);
  const double* y = REAL(VECTOR_ELT(dataExpr, matchPos[4]));
  double lb = REAL(VECTOR_ELT(dataExpr, matchPos[2]))[0];
  double ub = REAL(VECTOR_ELT(dataExpr, matchPos[3]))[0];
  
  int numNonZero = rc_getInt0(VECTOR_ELT(dataExpr, matchPos[23]), "num_non_zero");
  const int* v = INTEGER(VECTOR_ELT(dataExpr, matchPos[25]));
  const int* u = INTEGER(VECTOR_ELT(dataExpr, matchPos[26]));
  
  int maxU = rc_getLength(VECTOR_ELT(dataExpr, matchPos[24])) + 1;
  int q = INTEGER(VECTOR_ELT(dataExpr, matchPos[15]))[0];
  
  misc_stackFree(matchPos);
  
  for (size_t i = 0; i < n; ++i) {
    if (y[i] < lb || y[i] > ub) {
      delete result;
      Rf_error("y[%lu] out of [%f, %f] range", i + 1, lb, ub);
    }
    if (u[i] < 0 || u[i] > maxU) {
      delete result;
      Rf_error("u[%lu] out of [0, %d] range", i + 1, maxU);
    }
  }
  if (u[n] < 0 || u[n] > maxU) {
    delete result;
    Rf_error("u[%lu] out of [0, %d] range", n + 1, maxU);
  }
  for (int i = 0; i < numNonZero; ++i) {
    if (v[i] < 0 || v[i] > q - 1) {
      delete result;
      Rf_error("v[%d] out of [0, %d] range", i + 1, q - 1);
    }
  }
  
  return result;
}

void deleteStanModel(StanModel* model) {
  delete model;
}

void initializeStanArgsFromExpression(StanArgs& args, SEXP argsExpr)
{
  SEXP argNamesExpr = rc_getNames(argsExpr);
  if (Rf_isNull(argNamesExpr))
    Rf_error("names for stanArgs object cannot be NULL");
  
  size_t numArgNames = sizeof(argNames) / sizeof(argNames[0]);
  size_t* matchPos = misc_stackAllocate(numArgNames, size_t);
  
  size_t numInputArgNames = rc_getLength(argNamesExpr);
  const char** inputArgNames = misc_stackAllocate(numInputArgNames, const char*);
  for (size_t i = 0; i < numInputArgNames; ++i)
    inputArgNames[i] = CHAR(STRING_ELT(argNamesExpr, i));
  
  int errorCode = misc_str_matchAllInArray(argNames, numArgNames, inputArgNames, numInputArgNames, matchPos);
  misc_stackFree(inputArgNames);
  
  if (errorCode != 0) {
    misc_stackFree(matchPos);
    Rf_error("error matching names: %s", std::strerror(errorCode));
  }
  
  if (matchPos[0] == MISC_STR_NO_MATCH)
    Rf_error("stanArgs requires 'seed' to be specified");
  
  args.random_seed = static_cast<unsigned int>(rc_getInt0(VECTOR_ELT(argsExpr, matchPos[0]), "seed"));
  args.init_radius = rc_getDoubleAt(argsExpr, matchPos[1], "init_r",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 2.0, RC_END);
  args.thin = rc_getIntAt(argsExpr, matchPos[2], "thin",
    RC_VALUE | RC_GT, 0,
    RC_NA | RC_YES, RC_END);
  args.adapt_gamma = rc_getDoubleAt(argsExpr, matchPos[3], "adapt_gamma",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 0.05, RC_END);
  args.adapt_delta = rc_getDoubleAt(argsExpr, matchPos[4], "adapt_delta",
    RC_VALUE | RC_GT, 0.0,
    RC_VALUE | RC_LT, 1.0,
    RC_VALUE | RC_DEFAULT, 0.8, RC_END);
  args.adapt_kappa = rc_getDoubleAt(argsExpr, matchPos[5], "adapt_kappa",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 0.75, RC_END);
  args.adapt_init_buffer = static_cast<unsigned int>(rc_getIntAt(argsExpr, matchPos[6], "adapt_init_buffer",
    RC_VALUE | RC_DEFAULT, 75u, RC_END));
  args.adapt_term_buffer = static_cast<unsigned int>(rc_getIntAt(argsExpr, matchPos[7], "adapt_term_buffer",
    RC_VALUE | RC_DEFAULT, 50u, RC_END));
  args.adapt_window = static_cast<unsigned int>(rc_getIntAt(argsExpr, matchPos[8], "adapt_window",
    RC_VALUE | RC_DEFAULT, 25u, RC_END));
  args.adapt_t0 = rc_getDoubleAt(argsExpr, matchPos[9], "adapt_t0",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 10.0, RC_END);
  args.stepsize = rc_getDoubleAt(argsExpr, matchPos[10], "stepsize",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_DEFAULT, 1.0, RC_END);
  args.stepsize_jitter = rc_getDoubleAt(argsExpr, matchPos[11], "stepsize_jitter",
    RC_VALUE | RC_GEQ, 0.0,
    RC_VALUE | RC_LEQ, 1.0,
    RC_VALUE | RC_DEFAULT, 0.0, RC_END);
  args.max_treedepth = rc_getIntAt(argsExpr, matchPos[12], "max_treedepth",
    RC_VALUE | RC_GEQ, 0,
    RC_VALUE | RC_DEFAULT, 10, RC_END);
  
  misc_stackFree(matchPos);
}

StanSampler::StanSampler(StanModel& stanModel, const StanArgs& stanArgs, int chain_id, int num_warmup) :
  logger(nullout, nullout, nullout, nullout, nullout),
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
  sampler_names[1] = "treedepth__";
  sampler_names[2] = "n_leapfrog__";
  sampler_names[3] = "divergent__";
  sampler_names[4] = "energy__";
  sample_writer_offset = sample_names.size() + sampler_names.size();
  
  init_writer.resize(stanModel.num_params_r(), 1);
  
  num_pars = sample_names.size() + sampler_names.size() + constrained_param_names.size();
  // sample_writer.resize(sample_names.size() + sampler_names.size() + constrained_param_names.size(), num_warmup + num_iter);
  
  sampler = new stan4bart::interruptable_sampler<StanModel>(
    stanModel,
    *init_context_ptr,
    unit_e_metric,
    stanArgs.random_seed,
    chain_id,
    stanArgs.init_radius,
    num_warmup,
    stanArgs.thin,
    stanArgs.stepsize,
    stanArgs.stepsize_jitter,
    stanArgs.max_treedepth,
    stanArgs.adapt_delta,
    stanArgs.adapt_gamma,
    stanArgs.adapt_kappa,
    stanArgs.adapt_t0,
    stanArgs.adapt_init_buffer,
    stanArgs.adapt_term_buffer,
    stanArgs.adapt_window,
    interrupt, 
    logger,
    init_writer,
    sample_writer,
    diagnostic_writer);
}

void setStanOffset(StanModel& model, const double* offset)
{
  model.set_offset(offset);
}

void getZb(const StanSampler& sampler, const StanModel& model, double* Zb)
{
  model.get_Zb(sampler.sample_writer.x_curr + sampler.sample_writer_offset, Zb);
}

double getSigma(const StanSampler& sampler, const StanModel& model)
{
  return model.get_aux(sampler.sample_writer.x_curr + sampler.sample_writer_offset);
}

void StanSampler::run(bool isWarmup)
{
  sampler->run(isWarmup);
}

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

void printStanArgs(const StanArgs& args)
{
  Rprintf("  seed: %u\n"
          "  init_r: %f\n"
          "  thin: %d\n"
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
          args.random_seed, args.init_radius,
          args.thin, args.adapt_gamma, args.adapt_delta,
          args.adapt_kappa, args.adapt_init_buffer, args.adapt_term_buffer,
          args.adapt_window, args.adapt_t0, args.stepsize,
          args.stepsize_jitter, args.max_treedepth);
}

}

