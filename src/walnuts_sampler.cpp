/// \file walnuts_sampler.cpp
/// \brief WALNUTS wiring for the parametric conditional, plus
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
                               double init_radius, int num_warmup,
                               double step_accept_rate_target, bool save_raw) {
  impl_ = new Impl();
  impl_->model = buildParametricModel(dataExpr);
  impl_->model.save_raw = save_raw;
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

  // step_accept_rate_target is the Adam acceptance-rate target (adapt_delta);
  // its 0.8 default matches WarmupConfigBuilder's own, so an unset adapt_delta
  // builds a bit-identical config to the historical fixed-target path.
  impl_->warmup_cfg = walnuts::WarmupConfigBuilder()
                          .min_max_iter(num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1,
                                        num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1)
                          .step_accept_rate_target(step_accept_rate_target)
                          .build();
  impl_->sampling_cfg = walnuts::SamplingConfigBuilder().build();

  // Row layout + names. By default only the two LIVE diagnostics (lp__ and
  // stepsize__) lead each row; save_raw restores the full Stan-identical
  // header (lp__, accept_stat__, then the five constant-zero placeholder
  // sampler rows) alongside the raw constrained block.
  if (save_raw) {
    sample_names = {"lp__", "accept_stat__"};
    sampler_names = {"stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"};
  } else {
    sample_names = {"lp__"};
    sampler_names = {"stepsize__"};
  }
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

  // Seed adaptation's initial mass from upstream's callable Nutpie heuristic
  // instead of the historical identity metric: mass = (1 - s)|grad| + s
  // (config.hpp InitConfigBuilder::masses, one gradient eval), evaluated at the
  // chain's initial unconstrained position. That position is the
  // uniform(-init_radius, init_radius) draw, NOT a centered zero: unconstrained
  // zero is exactly the onion dot_self == 0 singularity (parametric_model.hpp
  // make_theta_L, sf = sqrt(rho / D) with D == 0) that the random start is
  // chosen to avoid, so seeding at zero would throw and fall back for every
  // nc >= 3 block. The actual init point is finite for every reachable model
  // and stays reproducible (the draw is a deterministic function of the chain
  // seed). The smoothing s reuses upstream's own mass additive-smoothing
  // default (WarmupConfig::mass_additive_smoothing, 1e-5) - the same
  // (1 - s)x + s interpolation the online mass estimator applies - so a
  // near-zero gradient coordinate keeps a positive floor rather than a
  // degenerate zero mass. mass_init_count (4) regularizes the whole warmup mass
  // estimate toward this seed, so it helps well past t = 0.
  //
  // The initial STEP stays at the conservative 0.1 constant; upstream's step
  // probe (adapt_step_build) is deliberately NOT run. Probing against the
  // seeded mass at the construction-time target - conditioned on the initial
  // BART offset / probit latents, both far from where the Gibbs sweeps settle -
  // returns steps 2-4x Adam's equilibrium, an oversized start measured to lock
  // ~15% of binary-tier chains into a full-sampling-phase rejection state
  // (ess ~= 2; 5/32 chains across four seeds vs 0/32 unseeded, 0/24 with
  // mass-only). Adam recovers the step from 0.1 within tens of warmup
  // transitions, so the probe's upside is small and its tail risk is not.
  //
  // The seed's gradient eval runs HERE, at construction, before the first
  // warmup transition, so it folds into the warmup phase's
  // mean_leapfrog_warmup (evals_warmup counts it) and leaves the
  // sampling-phase mean_leapfrog clean (the eval counters assembled into
  // fit$adaptation by stan4bart_fit.R).
  //
  // Robust fallback: the model throws std::domain_error on a non-finite log
  // density or gradient (parametric_model.hpp), matching NoExceptLogpGrad's
  // divergence semantics. A poisoned seed eval falls back to the historical
  // identity mass rather than seeding from a non-finite value.
  //
  // The seed only moves the STARTING point of adaptation; the freeze ->
  // fixed-draws lifecycle is untouched. Config objects are Impl members,
  // referenced by the adapter, so they outlive it.
  const double mass_smoothing = impl_->warmup_cfg.mass_additive_smoothing();
  std::optional<walnuts::InitChainConfig> init;
  try {
    walnuts::InitConfig seeded =
        walnuts::InitConfigBuilder(1u, static_cast<std::size_t>(dim))
            .positions(impl_->position)
            .masses(impl_->model, mass_smoothing)   // one gradient eval
            .build();
    init.emplace(impl_->step_size, impl_->position, seeded.mass(0u));
  } catch (...) {
    init.emplace(impl_->step_size, impl_->position, Eigen::VectorXd::Ones(dim));
  }
  impl_->adapter.emplace(impl_->rng, impl_->handler, impl_->model, *init,
                         impl_->warmup_cfg, impl_->sampling_cfg);
}

WalnutsSampler::~WalnutsSampler() { delete impl_; }

void WalnutsSampler::run(bool isWarmup) {
  if (isWarmup)
    (*impl_->adapter)();   // one adapting transition
  else
    (*impl_->sampler)();   // one fixed-tuning transition
  impl_->position = impl_->handler.position;

  // Re-emit the constrained draw into the current writer row.
  double* row = sample_writer.x_curr;
  double logp = 0.0;
  impl_->model.eval(impl_->position, logp, impl_->grad, row + sample_writer_offset);
  if (impl_->model.save_raw) {
    row[0] = logp;                      // lp__
    row[1] = 1.0;                       // accept_stat__ (WALNUTS reports no analog)
    row[2] = impl_->handler.step_size;  // stepsize__
    row[3] = 0.0;                       // treedepth__
    row[4] = 0.0;                       // n_leapfrog__
    row[5] = 0.0;                       // divergent__ (WALNUTS delivers no divergences)
    row[6] = 0.0;                       // energy__
  } else {
    row[0] = logp;                      // lp__
    row[1] = impl_->handler.step_size;  // stepsize__
  }
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

// The handler captured these at freeze() (AdaptiveWalnuts::sampler() fires
// on_warmup_complete with the frozen step size + diagonal inverse mass); valid
// only once disengageAdaptation has run.
double WalnutsSampler::getStepSize() const { return impl_->handler.step_size; }

int WalnutsSampler::getAdaptDim() const { return impl_->model.dim(); }

void WalnutsSampler::getInvMass(double* out) const {
  const Eigen::VectorXd& m = impl_->handler.inv_mass;
  const int n = impl_->model.dim();
  // m.size() == n after freeze(); the guard is purely defensive.
  for (int i = 0; i < n; ++i)
    out[i] = i < m.size() ? m[i] : 0.0;
}

long long WalnutsSampler::getEvalCount() const { return impl_->model.evalCount(); }

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
