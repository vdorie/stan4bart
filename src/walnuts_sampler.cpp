/// \file walnuts_sampler.cpp
/// \brief WALNUTS wiring for the parametric conditional, plus
///        the shared ParametricModel builder consumed by the gradient gate.

#include "walnuts_sampler.hpp"
#include "parametric_model_io.hpp"

#include <cmath>     // exp, log, sqrt, isfinite
#include <cstring>   // memcpy
#include <exception>
#include <limits>
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
  // Required by the SampleHandler concept. The vendored NoExceptLogpGrad
  // reports a thrown log density here and then returns -inf with a zero
  // gradient, which is the divergence the model's own domain_error already
  // means; nothing further is recorded, and this must not throw because the
  // caller is noexcept.
  void on_logp_exception(const Eigen::VectorXd&, const std::exception&) noexcept {}
};

using Adapter = walnutpie::AdaptiveWalnuts<ParametricModel, std::mt19937_64, LatestDraw>;
using Sampler = walnutpie::WalnutsSampler<ParametricModel, std::mt19937_64, LatestDraw>;

/// \brief One draw from a scalar log-concave density by Neal's (2003)
///        stepping-out slice sampler.
///
/// `logDensity` need be correct only up to an additive constant and must be
/// finite at `x0`, where its value is `logp0`. `width` must not depend on the
/// current position: the interval it seeds is what makes the draw reversible.
/// Log-concavity bounds both loops; their caps are guards, not policy.
template <class F, class RNG>
double sliceDraw(const F& logDensity, double x0, double logp0, double width,
                 RNG& rng) {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::exponential_distribution<double> expo(1.0);
  const double level = logp0 - expo(rng);

  // Figure 3: a window of the given width placed uniformly around x0, stepped
  // out at most m times with the budget split at random between the two sides.
  const int m = 64;
  double lo = x0 - width * unif(rng);
  double hi = lo + width;
  int left = static_cast<int>(m * unif(rng));
  int right = m - 1 - left;
  while (left > 0 && logDensity(lo) > level) { lo -= width; --left; }
  while (right > 0 && logDensity(hi) > level) { hi += width; --right; }

  // Figure 5: shrink until a point clears the slice, each rejection replacing
  // the endpoint on its own side of x0 so that x0 stays inside the interval.
  for (int iter = 0; iter < 200; ++iter) {
    const double x = lo + unif(rng) * (hi - lo);
    if (logDensity(x) > level) return x;
    if (x < x0) lo = x; else hi = x;
  }
  return x0;  // the interval underflowed to a point; stay put
}
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
  walnutpie::WarmupConfig   warmup_cfg   = walnutpie::WarmupConfigBuilder().build();
  walnutpie::SamplingConfig sampling_cfg = walnutpie::SamplingConfigBuilder().build();

  std::optional<Adapter> adapter;   // live during warmup
  std::optional<Sampler> sampler;   // live after freeze()

  // WALNUTS caches the log density and its gradient at its current position
  // from one transition to the next. This package moves the target between
  // transitions - the BART fit becomes the parametric offset, and for binary
  // responses the latents become the response - so that cache goes stale on
  // every sweep and the next transition would score its initial state under
  // the previous sweep's target. Set here, cleared by one refresh at the top
  // of the next run(): one evaluation per sweep rather than one per
  // transition, so skip > 1 pays it once.
  bool logp_stale = false;

  // The tuning the frozen sampler was built with. WALNUTS holds its position
  // privately, so the ridge move's setter is a rebuild at the moved position,
  // and a rebuild has to restate every tuning value. Three are readable off
  // the sampler and two off this package's own SamplingConfig, but
  // min_micro_steps is estimated during warmup and freeze() drops the adapter
  // that holds it, so it must be captured before the adapter goes.
  Eigen::VectorXd frozen_inv_mass;
  double frozen_macro_time = 0.0;
  double frozen_max_error = 0.0;
  std::size_t frozen_max_nuts_depth = 0;
  std::size_t frozen_max_step_halvings = 0;
  std::size_t frozen_min_micro_steps = 0;
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
  impl_->warmup_cfg = walnutpie::WarmupConfigBuilder()
                          .min_max_iter(num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1,
                                        num_warmup > 0 ? static_cast<std::size_t>(num_warmup) : 1)
                          .step_accept_rate_target(step_accept_rate_target)
                          .build();
  impl_->sampling_cfg = walnutpie::SamplingConfigBuilder().build();

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
  std::optional<walnutpie::InitChainConfig> init;
  try {
    walnutpie::InitConfig seeded =
        walnutpie::InitConfigBuilder(1u, static_cast<std::size_t>(dim))
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
  if (impl_->logp_stale) {
    if (isWarmup)
      impl_->adapter->refresh_logp_grad();
    else
      impl_->sampler->refresh_logp_grad();
    impl_->logp_stale = false;
  }
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

/// \brief One slice draw of every random-effect block's scale along the curve
///        that holds the linear predictor fixed, then a rebuild of the frozen
///        sampler at the moved position.
///
/// THE CURVE. Every block's Cholesky factor is homogeneous of degree one in
/// that block's own scale. In ParametricModel::eval the scale is
/// s_i = tau_i * re_scale_i * dispersion, and every entry of T_i is a multiple
/// of s_i: a scalar block is T = s_i; a correlated block has trace = nc s_i^2,
/// so sd_c = sqrt(pi_c * trace) = s_i sqrt(nc pi_c), and each onion entry is a
/// fixed function of rho, zeta and z_T times one sd. With b_level = T_i
/// z_b_level, the map
///
///     tau_i -> tau_i e^u,   z_b(block i) -> z_b(block i) e^-u
///
/// therefore leaves b, the linear predictor, and the likelihood exactly where
/// they were. That holds for every reachable block shape - several
/// random-effect terms (each block has its own tau and its own z_b segment,
/// and the moves are separate conditionals), correlated slopes under the decov
/// structure, and the binary family, whose latents enter only as the response
/// of that same invariant likelihood. re_scale and dispersion never enter the
/// conditional: they multiply tau_i and z_b in one product.
///
/// THE CONDITIONAL. Sample along the curve as a Gibbs step in the chart that
/// holds its invariant fixed. Write x = log tau_i (which IS the unconstrained
/// coordinate the sampler carries) and v = z_b(block i) in R^{q_i}, where
/// q_i = p_i * l_i counts every coordinate of every level of the block, not
/// the coordinates per level. Change variables (x, v) -> (x, w) with
/// w = e^x v = tau_i v, the quantity the curve preserves. The Jacobian of
/// v -> w at fixed x is e^{q_i x} times the identity, so the density in the new
/// chart carries a factor e^{-q_i x}. The only terms of eval()'s log density
/// that involve x or v are the standardized effects' N(0, 1) prior, tau's
/// Gamma(shape_i, 1) prior, and the log Jacobian of tau = e^x; the likelihood,
/// beta, rho, zeta, z_T and aux are all invariant. Collecting them,
///
///     log p(x | w, rest) = (shape_i - q_i) x - e^x - (A_i / 2) e^{-2x},
///     A_i = ||w||^2 = tau_i^2 ||z_b(block i)||^2,
///
/// with (shape_i - 1) x - e^x the Gamma prior, +x its Jacobian, -q_i x the
/// chart's, and the last term the N(0, 1) prior at v = e^{-x} w. A_i is
/// invariant along the curve, so it is computed once from the current state.
/// The second derivative is -e^x - 2 A_i e^{-2x} < 0, so the conditional is
/// log-concave on all of R with no condition on shape or on q_i, and a
/// stepping-out slice sampler draws it exactly: O(q_i) to form A_i, then a
/// handful of scalar evaluations and no gradient.
///
/// STATIONARITY. The move conditions on nothing the likelihood supplies, so it
/// leaves the parametric conditional invariant wherever in the sweep it is
/// applied. It runs at the top of the sweep, after the BART draw has set the
/// offset (and, for binary responses, the latent response) and before the
/// WALNUTS transitions, because that is where the rebuild's own log density
/// and gradient evaluation is the one the stale cache already owed.
void WalnutsSampler::ridgeMove() {
  Impl& impl = *impl_;
  const ParametricModel& model = impl.model;
  if (model.t == 0 || !impl.sampler) return;

  Eigen::VectorXd& theta = impl.position;
  bool moved = false;
  int b_off = 0;
  for (int i = 0; i < model.t; ++i) {
    const int q_i = model.p[i] * model.l[i];
    if (q_i <= 0) continue;
    double* z = theta.data() + model.off_z_b + b_off;
    b_off += q_i;

    double zz = 0.0;
    for (int j = 0; j < q_i; ++j) zz += z[j] * z[j];
    if (!(zz > 0.0)) continue;  // the curve degenerates at z_b == 0

    const double x0 = theta[model.off_tau + i];
    const double logA = 2.0 * x0 + std::log(zz);
    const double coef = model.tau_shape[i] - static_cast<double>(q_i);
    // Overflow in either exponential is the density underflowing to zero in
    // one tail or the other, which is what -inf means to the slice.
    const auto logDensity = [coef, logA](double x) {
      const double v = coef * x - std::exp(x) - 0.5 * std::exp(logA - 2.0 * x);
      return std::isfinite(v) ? v : -std::numeric_limits<double>::infinity();
    };
    const double logp0 = logDensity(x0);
    if (!std::isfinite(logp0)) continue;

    // -d^2/dx^2 = e^x + 2 A_i e^{-2x} is 2 q_i at the mode up to the scale's
    // own size, so a few posterior sds is a few over sqrt(q_i). The width is a
    // function of the block's geometry only: one that read the current
    // position would cost the draw its reversibility.
    const double width = 2.0 / std::sqrt(static_cast<double>(q_i));
    const double x1 = sliceDraw(logDensity, x0, logp0, width, impl.rng);
    if (x1 == x0) continue;

    theta[model.off_tau + i] = x1;
    const double shrink = std::exp(x0 - x1);
    for (int j = 0; j < q_i; ++j) z[j] *= shrink;
    moved = true;
  }
  if (!moved) return;

  // WALNUTS holds its position privately, so the setter is a rebuild at the
  // moved position with the tuning captured at freeze(). The base generator is
  // held by reference and continues where it was; the rebuilt sampler's own
  // normal_distribution is fresh, and it caches a spare variate, so a sweep
  // that rebuilds is not bitwise what a persistent sampler would have drawn.
  impl.sampler.emplace(impl.rng, impl.handler, impl.model, theta,
                       impl.frozen_inv_mass, impl.frozen_macro_time,
                       impl.frozen_max_nuts_depth,
                       impl.frozen_max_step_halvings,
                       impl.frozen_min_micro_steps, impl.frozen_max_error);
  // That constructor evaluates the log density and gradient at the handed
  // position under the live target, which is the refresh this sweep owed.
  impl.logp_stale = false;
}

void WalnutsSampler::freeze() {
  // Read before the handoff: sampler() drops the adapter's min-micro estimator
  // and WalnutsSampler exposes no getter for the value it was handed.
  impl_->frozen_min_micro_steps = impl_->adapter->min_micro_steps();
  impl_->sampler.emplace(impl_->adapter->sampler());  // AdaptiveWalnuts -> WalnutsSampler
  impl_->adapter.reset();
  impl_->frozen_inv_mass = impl_->sampler->inverse_mass_matrix_diagonal();
  impl_->frozen_macro_time = impl_->sampler->macro_time();
  impl_->frozen_max_error = impl_->sampler->max_error();
  impl_->frozen_max_nuts_depth = impl_->sampling_cfg.max_trajectory_doublings();
  impl_->frozen_max_step_halvings = impl_->sampling_cfg.max_step_halvings();
  // The frozen sampler's constructor evaluates the log density and gradient at
  // the handed-over position, so it starts from the live target.
  impl_->logp_stale = false;
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
  impl_->logp_stale = true;
}

void WalnutsSampler::setResponse(const double* y) {
  std::memcpy(impl_->model.y_.data(), y,
              static_cast<size_t>(impl_->model.N) * sizeof(double));
  impl_->logp_stale = true;
}

void WalnutsSampler::setVerbose(int /*level*/) {}

}  // namespace stan4bart
