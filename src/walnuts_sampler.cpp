/// \file walnuts_sampler.cpp
/// \brief WALNUTS wiring for the continuous parametric conditional (C2), plus
///        the shared ParametricModel builder consumed by the gradient gate.

#include "walnuts_sampler.hpp"
#include "parametric_model_io.hpp"

#include <cstring>   // memcpy
#include <optional>
#include <random>
#include <string>
#include <vector>

// The WALNUTS + Eigen headers trip warnings we neither own nor can fix; silence
// them just around these includes (pattern from bairrtt_types.h / stan4bart's
// interruptable_sampler.hpp).
#if (defined(__clang__) && (__clang_major__ > 3 || (__clang_major__ == 3 && __clang_minor__ >= 7))) || \
    (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#  define S4B_WALNUTS_SUPPRESS_DIAGNOSTIC 1
#endif

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 1
#ifdef S4B_WALNUTS_SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wunknown-pragmas"
#    pragma clang diagnostic ignored "-Wunused-variable"
#    pragma clang diagnostic ignored "-Wunused-parameter"
#    pragma clang diagnostic ignored "-Wsign-compare"
#    pragma clang diagnostic ignored "-Wignored-qualifiers"
#    pragma clang diagnostic ignored "-Wshorten-64-to-32"
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wunknown-pragmas"
#    pragma GCC diagnostic ignored "-Wunused-variable"
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#    pragma GCC diagnostic ignored "-Wsign-compare"
#    pragma GCC diagnostic ignored "-Wignored-qualifiers"
#  endif
#endif

#include <walnuts/adaptive_walnuts.hpp>
#include <walnuts/config.hpp>
#include <walnuts/walnuts.hpp>

#ifdef S4B_WALNUTS_SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif

namespace stan4bart {

// ---- shared ParametricModel builder ---------------------------------------
// Reads the marshaled Stan data list (a named R list, the same one
// continuous_model consumes) and replicates init.cpp's transformed-data
// derivations (len_z_T, len_rho, delta).

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

}  // namespace

ParametricModel buildParametricModel(SEXP dataExpr) {
  ParametricModel pm;

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

// ---- the WALNUTS sampler ----------------------------------------------------

// Minimal ChainHandler (bairrtt's LatestDraw): captures the most recent draw
// and the latest tuning out of WALNUTS. On_sample fires per post-warmup draw;
// on_warmup per warmup draw; on_warmup_complete at freeze().
namespace {
struct LatestDraw {
  Eigen::VectorXd position;
  double lp = 0.0;
  double step_size = 0.0;
  Eigen::VectorXd inv_mass;
  void on_sample(const Eigen::VectorXd& p, double l) { position = p; lp = l; }
  void on_warmup(const Eigen::VectorXd& p, double l, double s, const Eigen::VectorXd& m) {
    position = p; lp = l; step_size = s; inv_mass = m;
  }
  void on_warmup_complete(double s, const Eigen::VectorXd& m) { step_size = s; inv_mass = m; }
};

using Adapter = walnuts::AdaptiveWalnuts<ParametricModel, std::mt19937_64, LatestDraw>;
using Sampler = walnuts::WalnutsSampler<ParametricModel, std::mt19937_64, LatestDraw>;
}  // namespace

struct WalnutsSampler::Impl {
  ParametricModel model;       // holds X, y_, offset_, ... ; refreshed in place
  LatestDraw handler;          // captures each draw / latest tuning
  std::mt19937_64 rng;
  Eigen::VectorXd position;
  Eigen::VectorXd grad;        // scratch for the per-draw re-emission
  double step_size = 0.1;

  // WALNUTS holds these two by const reference; they must outlive the adapter,
  // hence they are declared before (destroyed after) the adapter/sampler.
  walnuts::WarmupConfig   warmup_cfg   = walnuts::WarmupConfigBuilder().build();
  walnuts::SamplingConfig sampling_cfg = walnuts::SamplingConfigBuilder().build();

  std::optional<Adapter> adapter;   // live during warmup
  std::optional<Sampler> sampler;   // live after freeze()
};

WalnutsSampler::WalnutsSampler(SEXP dataExpr, unsigned int random_seed,
                               double init_radius, int num_warmup) {
  impl_ = new Impl();
  impl_->model = buildParametricModel(dataExpr);
  const int dim = impl_->model.dim();

  // Seed the WALNUTS rng from the per-chain seed exactly as Stan's was threaded
  // (control.stan$seed, drawn deterministically from the master seed), so
  // same seed -> identical draws (test-05-rng).
  impl_->rng.seed(static_cast<std::mt19937_64::result_type>(random_seed));

  // Initial position: uniform(-init_radius, init_radius) per coordinate from
  // the seeded rng, mirroring Stan's init_r random start (also keeps z_T off
  // the onion's dot_self == 0 singularity).
  impl_->position.resize(dim);
  std::uniform_real_distribution<double> unif(-init_radius, init_radius);
  for (int i = 0; i < dim; ++i) impl_->position[i] = unif(impl_->rng);
  impl_->grad.resize(dim);

  impl_->warmup_cfg = walnuts::WarmupConfigBuilder()
                          .min_max_iter(num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1,
                                        num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1)
                          .build();
  impl_->sampling_cfg = walnuts::SamplingConfigBuilder().build();

  // Stan-identical row layout + names, so the R side consumes by name unchanged.
  sample_names = {"lp__", "accept_stat__"};
  sampler_names = {"stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"};
  constrained_param_names = impl_->model.constrainedParamNames();
  sample_writer_offset = sample_names.size() + sampler_names.size();
  num_pars = static_cast<int>(sample_writer_offset + constrained_param_names.size());
  sample_writer.num_pars = num_pars;
  {
    std::vector<std::string> all;
    all.reserve(static_cast<size_t>(num_pars));
    all.insert(all.end(), sample_names.begin(), sample_names.end());
    all.insert(all.end(), sampler_names.begin(), sampler_names.end());
    all.insert(all.end(), constrained_param_names.begin(), constrained_param_names.end());
    sample_writer.names = std::move(all);
  }

  // Start adaptation: identity initial mass (Stan's unit_e initial metric),
  // user step size (Adam adapts from here). Config objects are Impl members,
  // referenced by the adapter, so they outlive it.
  walnuts::InitChainConfig init(impl_->step_size, impl_->position,
                                Eigen::VectorXd::Ones(dim));
  impl_->adapter.emplace(impl_->rng, impl_->handler, impl_->model, init,
                         impl_->warmup_cfg, impl_->sampling_cfg);
}

WalnutsSampler::~WalnutsSampler() { delete impl_; }

void WalnutsSampler::run(bool isWarmup) {
  if (isWarmup)
    (*impl_->adapter)();   // one adapting transition
  else
    (*impl_->sampler)();   // one fixed-tuning transition
  impl_->position = impl_->handler.position;

  // Re-emit the constrained draw into the current writer row, Stan layout.
  double* row = sample_writer.x_curr;
  double logp = 0.0;
  impl_->model.eval(impl_->position, logp, impl_->grad, row + sample_writer_offset);
  row[0] = logp;                      // lp__
  row[1] = 1.0;                       // accept_stat__ (WALNUTS reports no analog)
  row[2] = impl_->handler.step_size;  // stepsize__
  row[3] = 0.0;                       // treedepth__
  row[4] = 0.0;                       // n_leapfrog__
  row[5] = 0.0;                       // divergent__ (WALNUTS delivers no divergences)
  row[6] = 0.0;                       // energy__
}

void WalnutsSampler::freeze() {
  impl_->sampler.emplace(impl_->adapter->sampler());  // AdaptiveWalnuts -> WalnutsSampler
  impl_->adapter.reset();
}

void WalnutsSampler::getParametricMean(double* result) const {
  impl_->model.parametricMean(sample_writer.x_curr + sample_writer_offset, result, true, true);
}

void WalnutsSampler::getParametricMean(double* result, bool includeFixed,
                                       bool includeRandom) const {
  impl_->model.parametricMean(sample_writer.x_curr + sample_writer_offset, result,
                              includeFixed, includeRandom);
}

double WalnutsSampler::getSigma() const {
  return impl_->model.getAux(sample_writer.x_curr + sample_writer_offset);
}

void WalnutsSampler::setOffset(const double* offset) {
  std::memcpy(impl_->model.offset_.data(), offset,
              static_cast<size_t>(impl_->model.N) * sizeof(double));
}

void WalnutsSampler::setResponse(const double* y) {
  std::memcpy(impl_->model.y_.data(), y,
              static_cast<size_t>(impl_->model.N) * sizeof(double));
}

void WalnutsSampler::setVerbose(int /*level*/) {}

}  // namespace stan4bart
