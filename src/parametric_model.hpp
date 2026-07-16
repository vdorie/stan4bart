#ifndef STAN4BART_PARAMETRIC_MODEL_HPP
#define STAN4BART_PARAMETRIC_MODEL_HPP

/// \file parametric_model.hpp
/// \brief Hand-derived log-posterior and reverse-mode gradient for the
///        parametric conditional of stan4bart's Gaussian linear mixed model.
///
/// This is the WALNUTS target that replaces the autodiff'd continuous.stan
/// model. It reproduces continuous.stan's log-density (including its quirks:
/// the weighted normalizer's N-not-sum(w), the decov "onion" prior on the RE
/// covariance, and the aux->dispersion->RE-scale coupling) for the reachable
/// prior scope: Gaussian likelihood + {none, normal, student_t/cauchy} beta +
/// {none, normal, student_t, exponential} aux + decov random effects. The
/// shrinkage families (hs, hs_plus, laplace, lasso, product_normal) are NOT
/// implemented; construction throws on their prior_dist codes.
///
/// The gradient is a hand-written reverse-mode adjoint mirroring the forward
/// loops (NOT a general autodiff), gated by finite differences and by a
/// cross-check against Stan's own log_prob (see tests/testthat/test-12).
///
/// Free (unconstrained) position vector, in continuous.stan's `parameters`
/// declaration order (the same order Stan's log_prob reads):
///   [ z_beta(K), z_b(q), z_T(len_z_T), rho_free(len_rho),
///     zeta_free(len_concentration), tau_free(t), aux_unscaled_free(!binary) ]

#include <cassert>    // assert
#include <cmath>      // std::exp, std::log, std::sqrt
#include <cstring>    // std::memcpy
#include <stdexcept>  // std::invalid_argument
#include <string>     // std::string, std::to_string
#include <vector>

#include <Eigen/Dense>

namespace stan4bart {

/// \brief Cornish-Fisher polynomial mapping a standard-normal deviate to a
///        student-t deviate (continuous.stan:146-158). Used for the
///        student_t / cauchy fixed-effect coefficient prior.
inline double CFt(double z, double df) {
  const double z2 = z * z, z3 = z2 * z, z5 = z2 * z3, z7 = z2 * z5, z9 = z2 * z7;
  const double df2 = df * df, df3 = df2 * df, df4 = df2 * df2;
  return z + (z3 + z) / (4.0 * df) + (5.0 * z5 + 16.0 * z3 + 3.0 * z) / (96.0 * df2)
         + (3.0 * z7 + 19.0 * z5 + 17.0 * z3 - 15.0 * z) / (384.0 * df3)
         + (79.0 * z9 + 776.0 * z7 + 1482.0 * z5 - 1920.0 * z3 - 945.0 * z) / (92160.0 * df4);
}

/// \brief d CFt / dz, the reverse-mode local derivative of the CFt transform.
inline double dCFt(double z, double df) {
  const double z2 = z * z, z4 = z2 * z2, z6 = z2 * z4, z8 = z2 * z6;
  const double df2 = df * df, df3 = df2 * df, df4 = df2 * df2;
  return 1.0 + (3.0 * z2 + 1.0) / (4.0 * df)
         + (25.0 * z4 + 48.0 * z2 + 3.0) / (96.0 * df2)
         + (21.0 * z6 + 95.0 * z4 + 51.0 * z2 - 15.0) / (384.0 * df3)
         + (711.0 * z8 + 5432.0 * z6 + 7410.0 * z4 - 5760.0 * z2 - 945.0) / (92160.0 * df4);
}

/// \brief The parametric target. Conditioning data (X, y, weights, offset_,
///        CSR Z, prior constants, block geometry) are held by value so the
///        sampler can refresh y_ / offset_ in place between Gibbs transitions. Matches the
///        walnuts::LogpGrad concept: void operator()(const VectorXd&, double&,
///        VectorXd&) const.
struct ParametricModel {
  // --- dimensions ---
  int N = 0, K = 0, q = 0, t = 0;
  bool is_binary = false;
  bool has_weights = false;

  /// \brief When false (the storage default), eval() emits and
  ///        constrainedParamNames() names ONLY the transformed block consumers
  ///        read (aux, beta, b, theta_L); the raw unconstrained rows
  ///        (z_beta/z_b/z_T/rho/zeta/tau/aux_unscaled) are dropped. When true
  ///        (save_raw_parameters opt-in, funnel forensics) the full Stan
  ///        write_array block is emitted exactly as historically. See
  ///        docs/design/walnuts.md, storage policy.
  bool save_raw = false;

  // --- data (y_ and offset_ mutable for the per-sweep Gibbs refresh) ---
  Eigen::MatrixXd X;        // N x K, centered fixed-effect design
  Eigen::VectorXd y_;       // N
  Eigen::VectorXd offset_;  // N
  Eigen::VectorXd weights;  // N (empty when !has_weights)

  // --- fixed-effect coefficient prior ---
  int prior_dist = 1;             // 0 none, 1 normal, 2 student_t/cauchy
  Eigen::VectorXd prior_scale;    // K
  Eigen::VectorXd prior_mean;     // K
  Eigen::VectorXd prior_df;       // K (student_t/cauchy df)

  // --- residual sd (aux) prior ---
  int prior_dist_for_aux = 3;     // 0 none, 1 normal, 2 student_t, 3 exponential
  double prior_scale_for_aux = 1.0;
  double prior_mean_for_aux = 0.0;
  double prior_df_for_aux = 0.0;

  // --- random-effect (decov) geometry ---
  std::vector<int> p;             // t: correlated coefs per block
  std::vector<int> l;             // t: levels per block
  int len_theta_L = 0, len_z_T = 0, len_rho = 0, len_concentration = 0;
  Eigen::VectorXd re_scale;       // t: decov `scale`
  Eigen::VectorXd tau_shape;      // t: gamma shape for tau
  std::vector<double> delta;      // len_concentration: gamma shape for zeta
  std::vector<double> regularization;  // per block with p_i > 1

  // --- sparse Z (CSR, 0-based; eta += Z b) ---
  Eigen::VectorXd Zw;             // num_non_zero nonzero values
  std::vector<int> Zv;            // num_non_zero column indices
  std::vector<int> Zu;            // N + 1 row pointers

  // --- free-vector segment offsets (filled by finalize()) ---
  int off_z_beta = 0, off_z_b = 0, off_z_T = 0, off_rho = 0;
  int off_zeta = 0, off_tau = 0, off_aux = 0, dim_ = 0;

  // Beta-onion prior shapes for each rho entry (len_rho), precomputed.
  Eigen::VectorXd rho_a1, rho_a2;

  /// \brief Per-grouping-factor forward cache + reverse-pass adjoints for the
  ///        make_theta_L / make_b construction. Held across eval() calls so its
  ///        buffers are reused rather than reallocated (block geometry - nc,
  ///        level_count - is fixed by finalize(), so every buffer's size is
  ///        call-invariant). g_T/g_sd/g_pi are the reverse-pass scratch (g_T was
  ///        formerly a separate per-call std::vector<MatrixXd>).
  struct Block {
    int nc, level_count, b_off, tau_idx, rho_off, zeta_off;
    double s, trace, Zsum;
    Eigen::VectorXd pi, sd;
    Eigen::MatrixXd T;
    std::vector<double> sf, Dm;      // per onion row m (index m in 2..nc-1)
    std::vector<int> zT_row_start;   // per onion row m: start index into z_T
    Eigen::MatrixXd g_T;             // reverse adjoint of T
    Eigen::VectorXd g_sd, g_pi;      // reverse-pass per-block scratch
  };

  /// \brief Reusable eval() scratch. THREAD-SAFETY: one ParametricModel lives
  ///        per chain (walnuts_sampler.cpp Impl::model); WALNUTS holds it by
  ///        const reference (util.hpp NoExceptLogpGrad's std::cref), and eval()
  ///        is never entered concurrently on a single instance - chains run as
  ///        separate R processes (stan4bart_fit.R makeCluster), the per-chain
  ///        WALNUTS classes (AdaptiveWalnuts/WalnutsSampler) spawn no threads,
  ///        and the leapfrog evaluations and the re-emission eval run in
  ///        sequence. These mutable buffers are therefore reused across calls
  ///        (buffer reuse only; the arithmetic and its order are unchanged).
  ///        Sizes are call-invariant; buffers are sized on first use.
  mutable Eigen::VectorXd rho_, zeta_, tau_;              // len_rho, len_conc, t
  mutable Eigen::VectorXd beta_, dbeta_dz_;               // K
  mutable Eigen::VectorXd b_, eta_, g_eta_, g_beta_;      // q, N, N, K
  mutable Eigen::VectorXd g_b_, g_z_b_;                   // q
  mutable Eigen::VectorXd g_rho_, g_zeta_, g_tau_chain_, g_zT_chain_;
  mutable std::vector<Block> blocks_;                     // t

  /// \brief Count of gradient evaluations on the leapfrog hot path - eval()
  ///        calls with constrained == nullptr, i.e. every operator() the WALNUTS
  ///        integrator makes, one per leapfrog micro-step. The once-per-stored-
  ///        draw re-emission eval (constrained != nullptr) is excluded, so
  ///        (evals accrued over a phase) / (transitions in that phase) is the
  ///        mean leapfrog steps per transition - the quantity that IS the
  ///        parametric-step runtime. mutable, reached through std::cref like
  ///        the scratch buffers; one model per chain and chains in separate
  ///        processes make it race-free. A single increment on a member
  ///        disjoint from the Eigen scratch, so the __restrict hot-loop
  ///        codegen is unperturbed and draws stay bitwise identical.
  mutable long long eval_count_ = 0;

  /// \brief The cumulative leapfrog-hot-path eval count (see eval_count_).
  long long evalCount() const { return eval_count_; }

  /// \brief Validate the reachable scope and precompute segment offsets and
  ///        the rho Beta shapes. Call once after populating the data members.
  void finalize() {
    if (prior_dist != 0 && prior_dist != 1 && prior_dist != 2)
      throw std::invalid_argument(
          "parametric_model: unsupported beta prior_dist (only none/normal/"
          "student_t/cauchy are implemented under WALNUTS)");
    if (!is_binary &&
        prior_dist_for_aux != 0 && prior_dist_for_aux != 1 &&
        prior_dist_for_aux != 2 && prior_dist_for_aux != 3)
      throw std::invalid_argument("parametric_model: unsupported aux prior_dist");

    off_z_beta = 0;
    off_z_b    = off_z_beta + K;
    off_z_T    = off_z_b + q;
    off_rho    = off_z_T + len_z_T;
    off_zeta   = off_rho + len_rho;
    off_tau    = off_zeta + len_concentration;
    off_aux    = off_tau + t;
    dim_       = off_aux + (is_binary ? 0 : 1);

    // decov_lp's per-block Beta shapes for rho (continuous.stan:104-118).
    rho_a1 = Eigen::VectorXd::Zero(len_rho);
    rho_a2 = Eigen::VectorXd::Zero(len_rho);
    int pos_reg = 0, pos_rho = 0;
    for (int i = 0; i < t; ++i) {
      const int nc = p[i];
      if (nc > 1) {
        double nu = regularization[pos_reg++] + 0.5 * (nc - 2);
        rho_a1[pos_rho] = nu;
        rho_a2[pos_rho] = nu;
        for (int j = 2; j <= nc - 1; ++j) {  // 1-based j, local index j-1
          nu -= 0.5;
          rho_a1[pos_rho + (j - 1)] = 0.5 * j;
          rho_a2[pos_rho + (j - 1)] = nu;
        }
        pos_rho += nc - 1;
      }
    }
  }

  int dim() const { return dim_; }

  /// \brief The raw/unconstrained prefix (z_beta, z_b, z_T, rho, zeta, tau,
  ///        [aux_unscaled]) that precedes the transformed block; stored only
  ///        under save_raw.
  int rawConstrainedDim() const {
    return K + q + len_z_T + len_rho + len_concentration + t + (is_binary ? 0 : 1);
  }
  /// \brief The transformed block every consumer reads: [aux], beta, b, theta_L.
  int transformedConstrainedDim() const {
    return (is_binary ? 0 : 1) + K + q + len_theta_L;
  }

  /// \brief Length of the constrained parameter array emitted per stored draw.
  ///        Under save_raw this is the full Stan write_array block
  ///        (z_beta, z_b, z_T, rho, zeta, tau, [aux_unscaled, aux], beta, b,
  ///        theta_L); by default only the transformed block ([aux], beta, b,
  ///        theta_L) is stored.
  int constrainedDim() const {
    return (save_raw ? rawConstrainedDim() : 0) + transformedConstrainedDim();
  }

  /// \brief Index of beta.1 within the constrained array (b follows at +K).
  int betaConstrainedOffset() const {
    return (save_raw ? rawConstrainedDim() : 0) + (is_binary ? 0 : 1);
  }
  /// \brief Index of aux.1 (the residual sd) within the constrained array;
  ///        valid only when !is_binary. aux leads the transformed block.
  int auxConstrainedOffset() const {
    return save_raw ? rawConstrainedDim() : 0;
  }

  /// \brief The Stan-identical constrained parameter names, in write_array
  ///        order, for the reachable scope. Mirrors continuous.hpp's
  ///        constrained_param_names() so the R side consumes rows by name
  ///        exactly as under the Stan sampler.
  std::vector<std::string> constrainedParamNames() const {
    std::vector<std::string> names;
    names.reserve(static_cast<size_t>(constrainedDim()));
    auto emit = [&](const char* base, int n) {
      for (int i = 1; i <= n; ++i)
        names.emplace_back(std::string(base) + "." + std::to_string(i));
    };
    if (save_raw) {
      emit("z_beta", K);
      emit("z_b", q);
      emit("z_T", len_z_T);
      emit("rho", len_rho);
      emit("zeta", len_concentration);
      emit("tau", t);
      if (!is_binary) emit("aux_unscaled", 1);
    }
    if (!is_binary) emit("aux", 1);
    emit("beta", K);
    emit("b", q);
    emit("theta_L", len_theta_L);
    return names;
  }

  /// \brief The parametric mean X*beta + Z*b (offset-free), read from a
  ///        constrained array. Mirrors continuous.hpp::get_parametric_mean.
  void parametricMean(const double* constrained, double* result,
                      bool includeFixed, bool includeRandom) const {
    const double* betap = constrained + betaConstrainedOffset();
    const double* bp = betap + K;
    Eigen::VectorXd eta = Eigen::VectorXd::Zero(N);
    if (includeFixed && K > 0) {
      Eigen::Map<const Eigen::VectorXd> betaVec(betap, K);
      eta.noalias() += X * betaVec;
    }
    if (includeRandom && t > 0) {
      for (int i = 0; i < N; ++i) {
        double acc = 0.0;
        for (int k = Zu[i]; k < Zu[i + 1]; ++k) acc += Zw[k] * bp[Zv[k]];
        eta[i] += acc;
      }
    }
    std::memcpy(result, eta.data(), static_cast<size_t>(N) * sizeof(double));
  }

  /// \brief The residual sd (aux) read from a constrained array.
  double getAux(const double* constrained) const {
    return constrained[auxConstrainedOffset()];
  }

  /// \brief Evaluate the log-posterior and its gradient at the unconstrained
  ///        point `par`. Matches the walnuts::LogpGrad concept.
  void operator()(const Eigen::VectorXd& par, double& logp, Eigen::VectorXd& grad) const {
    eval(par, logp, grad, nullptr);
  }

  /// \brief The log-posterior, its gradient, and (when `constrained` is
  ///        non-null) the constrained parameter array in Stan write_array
  ///        order. C2's WALNUTS wrapper calls this once per stored draw to
  ///        re-emit beta/b/theta_L/aux with Stan-identical names; the sampling
  ///        hot path calls operator() (constrained == nullptr) so no re-emission
  ///        cost is paid per leapfrog. Parameter-independent constants are
  ///        dropped (so the value matches Stan's log_prob up to a constant).
  void eval(const Eigen::VectorXd& par, double& logp, Eigen::VectorXd& grad,
            double* constrained) const {
    // Count only the leapfrog hot path (operator() -> constrained == nullptr);
    // the re-emission eval carries a constrained buffer and is not a leapfrog.
    if (constrained == nullptr) ++eval_count_;
    assert(par.size() == dim_ && "eval: position vector size must equal dim()");
    if (grad.size() != dim_) grad.resize(dim_);
    grad.setZero();
    logp = 0.0;

    const double* z_beta_p  = par.data() + off_z_beta;
    const double* z_b_p     = par.data() + off_z_b;
    const double* z_T_p     = par.data() + off_z_T;
    const double* rho_f_p   = par.data() + off_rho;
    const double* zeta_f_p  = par.data() + off_zeta;
    const double* tau_f_p   = par.data() + off_tau;

    // --- constrain rho / zeta / tau / aux ---
    // (reused scratch bound to local names; arithmetic body unchanged)
    rho_.resize(len_rho); zeta_.resize(len_concentration); tau_.resize(t);
    Eigen::VectorXd& rho = rho_;
    Eigen::VectorXd& zeta = zeta_;
    Eigen::VectorXd& tau = tau_;
    for (int j = 0; j < len_rho; ++j)          rho[j]  = 1.0 / (1.0 + std::exp(-rho_f_p[j]));
    for (int j = 0; j < len_concentration; ++j) zeta[j] = std::exp(zeta_f_p[j]);
    for (int i = 0; i < t; ++i)                tau[i]  = std::exp(tau_f_p[i]);

    double au = 1.0, aux = 1.0, dispersion = 1.0, actual_aux = 1.0, daux_dau = 1.0;
    if (!is_binary) {
      au = std::exp(par[off_aux]);
      if (prior_dist_for_aux == 0) { aux = au; daux_dau = 1.0; }
      else if (prior_dist_for_aux <= 2) { aux = prior_scale_for_aux * au + prior_mean_for_aux; daux_dau = prior_scale_for_aux; }
      else { aux = prior_scale_for_aux * au; daux_dau = prior_scale_for_aux; }
      dispersion = aux;
      actual_aux = aux;
    }

    // --- beta = transform(z_beta) ---
    beta_.resize(K); dbeta_dz_.resize(K);
    Eigen::VectorXd& beta = beta_;
    Eigen::VectorXd& dbeta_dz = dbeta_dz_;
    for (int k = 0; k < K; ++k) {
      if (prior_dist == 1) { beta[k] = z_beta_p[k] * prior_scale[k] + prior_mean[k]; dbeta_dz[k] = prior_scale[k]; }
      else if (prior_dist == 2) { beta[k] = CFt(z_beta_p[k], prior_df[k]) * prior_scale[k] + prior_mean[k]; dbeta_dz[k] = dCFt(z_beta_p[k], prior_df[k]) * prior_scale[k]; }
      else { beta[k] = z_beta_p[k]; dbeta_dz[k] = 1.0; }
    }

    // --- forward make_theta_L: build lower-tri Cholesky factor T per block,
    //     caching intermediates for the reverse pass (Block type + reused
    //     buffers are members; geometry is call-invariant so resize is a no-op
    //     after the first call) ---
    blocks_.resize(t);
    std::vector<Block>& blocks = blocks_;

    int b_mark = 0, rho_mark = 0, zeta_mark = 0, zT_mark = 0;
    for (int i = 0; i < t; ++i) {
      Block& bl = blocks[i];
      bl.nc = p[i];
      bl.level_count = l[i];
      bl.b_off = b_mark;
      bl.tau_idx = i;
      bl.rho_off = rho_mark;
      bl.zeta_off = zeta_mark;
      const int nc = bl.nc;
      bl.s = tau[i] * re_scale[i] * dispersion;

      if (nc == 1) {
        bl.T.setConstant(1, 1, bl.s);
      } else {
        bl.trace = bl.s * bl.s * nc;
        bl.pi.resize(nc);
        double zsum = 0.0;
        for (int c = 0; c < nc; ++c) { bl.pi[c] = zeta[zeta_mark + c]; zsum += bl.pi[c]; }
        bl.Zsum = zsum;
        bl.pi /= zsum;
        zeta_mark += nc;

        bl.sd.resize(nc);
        for (int c = 0; c < nc; ++c) bl.sd[c] = std::sqrt(bl.pi[c] * bl.trace);

        bl.T.setZero(nc, nc);
        bl.T(0, 0) = bl.sd[0];
        const double T21 = 2.0 * rho[rho_mark] - 1.0;
        bl.T(1, 1) = bl.sd[1] * std::sqrt(1.0 - T21 * T21);
        bl.T(1, 0) = bl.sd[1] * T21;
        int block_rho = rho_mark + 1;  // rho used by onion rows

        bl.sf.assign(nc, 0.0);
        bl.Dm.assign(nc, 0.0);
        bl.zT_row_start.assign(nc, 0);
        for (int m = 2; m < nc; ++m) {          // 0-based row m, T_row length m
          double D = 0.0;
          bl.zT_row_start[m] = zT_mark;
          for (int c = 0; c < m; ++c) { const double v = z_T_p[zT_mark + c]; D += v * v; }
          const double rho_m = rho[block_rho];
          const double sf = std::sqrt(rho_m / D) * bl.sd[m - 1];
          bl.sf[m] = sf;
          bl.Dm[m] = D;
          for (int c = 0; c < m; ++c) bl.T(m, c) = z_T_p[zT_mark + c] * sf;
          bl.T(m, m) = std::sqrt(1.0 - rho_m) * bl.sd[m];
          zT_mark += m;
          ++block_rho;
        }
        rho_mark += nc - 1;
      }
      b_mark += nc * l[i];
    }

    // --- forward make_b: b_level = T_i * z_b_level ---
    b_.resize(q);
    Eigen::VectorXd& b = b_;
    for (const Block& bl : blocks) {
      const int nc = bl.nc;
      for (int lev = 0; lev < bl.level_count; ++lev) {
        const int base = bl.b_off + lev * nc;
        for (int r = 0; r < nc; ++r) {
          double acc = 0.0;
          for (int c = 0; c < nc; ++c) acc += bl.T(r, c) * z_b_p[base + c];
          b[base + r] = acc;
        }
      }
    }

    // --- constrained re-emission (Stan write_array order): only when the
    //     caller wants the stored draw, not on the leapfrog hot path. theta_L
    //     is the column-major vech of each block's Cholesky factor T
    //     (continuous.stan:51-55), so make_b reconstructs the same T. ---
    if (constrained != nullptr) {
      int o = 0;
      if (save_raw) {
        for (int k = 0; k < K; ++k)                  constrained[o++] = z_beta_p[k];
        for (int s = 0; s < q; ++s)                  constrained[o++] = z_b_p[s];
        for (int e = 0; e < len_z_T; ++e)            constrained[o++] = z_T_p[e];
        for (int j = 0; j < len_rho; ++j)            constrained[o++] = rho[j];
        for (int j = 0; j < len_concentration; ++j)  constrained[o++] = zeta[j];
        for (int i = 0; i < t; ++i)                  constrained[o++] = tau[i];
        if (!is_binary)                              constrained[o++] = au;  // aux_unscaled
      }
      if (!is_binary)                                constrained[o++] = aux; // residual sd
      for (int k = 0; k < K; ++k)                    constrained[o++] = beta[k];
      for (int s = 0; s < q; ++s)                    constrained[o++] = b[s];
      for (const Block& bl : blocks) {
        const int nc = bl.nc;
        for (int c = 0; c < nc; ++c)
          for (int r = c; r < nc; ++r) constrained[o++] = bl.T(r, c);
      }
    }

    // --- eta = offset_ + X beta + Z b ---
    // The reused scratch lives in members (reached through `this`), which
    // defeats the no-alias codegen the compiler got for free when these were
    // stack temporaries; the O(N) loops below re-establish it with restrict-
    // qualified local pointers into the (genuinely non-overlapping) buffers.
    // Purely a compiler hint - the iteration order, the values, and the
    // accumulation order are byte-for-byte identical (bitwise-gated).
    eta_ = offset_;
    Eigen::VectorXd& eta = eta_;
    if (K > 0) eta.noalias() += X * beta;
    {
      double* __restrict etap = eta_.data();
      const double* __restrict bp = b_.data();
      const double* __restrict Zwp = Zw.data();
      const int* __restrict Zvp = Zv.data();
      const int* __restrict Zup = Zu.data();
      for (int i = 0; i < N; ++i) {
        double acc = 0.0;
        for (int k = Zup[i]; k < Zup[i + 1]; ++k) acc += Zwp[k] * bp[Zvp[k]];
        etap[i] += acc;
      }
    }

    // --- Gaussian likelihood; g_eta = lambda * (w .* r) ---
    const double lambda = 1.0 / (actual_aux * actual_aux);
    g_eta_.resize(N);
    Eigen::VectorXd& g_eta = g_eta_;
    double S = 0.0;  // sum_i w_i r_i^2
    {
      const double* __restrict yp = y_.data();
      const double* __restrict etap = eta_.data();
      double* __restrict getap = g_eta_.data();
      const double* __restrict wp = has_weights ? weights.data() : nullptr;
      for (int i = 0; i < N; ++i) {
        const double r = yp[i] - etap[i];
        const double w = has_weights ? wp[i] : 1.0;
        S += w * r * r;
        getap[i] = lambda * w * r;
      }
    }
    logp += -N * std::log(actual_aux) - 0.5 * lambda * S;

    // --- backprop eta -> beta, b ---
    // g_beta = X^T g_eta into reused storage (same GEMV, no per-call temporary)
    if (K > 0) g_beta_.noalias() = X.transpose() * g_eta;
    else       g_beta_.resize(0);
    Eigen::VectorXd& g_beta = g_beta_;
    g_b_.setZero(q);
    Eigen::VectorXd& g_b = g_b_;
    {
      double* __restrict gbp = g_b_.data();
      const double* __restrict getap = g_eta_.data();
      const double* __restrict Zwp = Zw.data();
      const int* __restrict Zvp = Zv.data();
      const int* __restrict Zup = Zu.data();
      for (int i = 0; i < N; ++i)
        for (int k = Zup[i]; k < Zup[i + 1]; ++k) gbp[Zvp[k]] += Zwp[k] * getap[i];
    }

    // --- reverse make_b: g_T += g_b_level z_b_level^T (lower tri); g_z_b += T^T g_b_level ---
    // g_T for each block is reused per-block scratch (bl.g_T); g_z_b is reused.
    g_z_b_.setZero(q);
    Eigen::VectorXd& g_z_b = g_z_b_;
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
      Block& bl = blocks[bi];
      const int nc = bl.nc;
      bl.g_T.setZero(nc, nc);
      Eigen::MatrixXd& gT = bl.g_T;
      for (int lev = 0; lev < bl.level_count; ++lev) {
        const int base = bl.b_off + lev * nc;
        for (int r = 0; r < nc; ++r) {
          const double gbr = g_b[base + r];
          for (int c = 0; c <= r; ++c) gT(r, c) += gbr * z_b_p[base + c];     // lower tri only
          for (int c = 0; c < nc; ++c) g_z_b[base + c] += bl.T(r, c) * gbr;   // T^T g_b
        }
      }
    }

    // --- reverse make_theta_L: g_T -> {dispersion, tau, zeta, rho, z_T} ---
    g_rho_.setZero(len_rho);
    g_zeta_.setZero(len_concentration);
    g_tau_chain_.setZero(t);
    g_zT_chain_.setZero(len_z_T);
    Eigen::VectorXd& g_rho = g_rho_;
    Eigen::VectorXd& g_zeta = g_zeta_;
    Eigen::VectorXd& g_tau_chain = g_tau_chain_;
    Eigen::VectorXd& g_zT_chain = g_zT_chain_;
    double g_disp = 0.0;
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
      Block& bl = blocks[bi];
      const int nc = bl.nc;
      const Eigen::MatrixXd& gT = bl.g_T;
      double g_s = 0.0;

      if (nc == 1) {
        g_s = gT(0, 0);
      } else {
        double g_trace = 0.0;
        bl.g_sd.setZero(nc);
        bl.g_pi.setZero(nc);
        Eigen::VectorXd& g_sd = bl.g_sd;
        Eigen::VectorXd& g_pi = bl.g_pi;

        // T(0,0) = sd[0]
        g_sd[0] += gT(0, 0);
        // row 1: T21 and the sqrt(1 - T21^2) diagonal
        const double T21 = 2.0 * rho[bl.rho_off] - 1.0;
        const double c22 = std::sqrt(1.0 - T21 * T21);
        double g_T21 = 0.0;
        g_sd[1] += gT(1, 0) * T21;
        g_T21   += gT(1, 0) * bl.sd[1];
        g_sd[1] += gT(1, 1) * c22;
        g_T21   += gT(1, 1) * bl.sd[1] * (-T21 / c22);
        g_rho[bl.rho_off] += g_T21 * 2.0;

        int block_rho = bl.rho_off + 1;
        for (int m = 2; m < nc; ++m) {
          const double sf = bl.sf[m], D = bl.Dm[m];
          const double rho_m = rho[block_rho];
          const int zt0 = bl.zT_row_start[m];
          double g_sf = 0.0;
          for (int c = 0; c < m; ++c) {
            const double v = z_T_p[zt0 + c];
            g_zT_chain[zt0 + c] += gT(m, c) * sf;
            g_sf += gT(m, c) * v;
          }
          // diagonal T(m,m) = sqrt(1 - rho_m) * sd[m]
          const double dmm = std::sqrt(1.0 - rho_m);
          g_sd[m] += gT(m, m) * dmm;
          g_rho[block_rho] += gT(m, m) * bl.sd[m] * (-0.5 / dmm);
          // sf = q * sd[m-1], q = sqrt(rho_m / D)
          const double qv = sf / bl.sd[m - 1];
          const double g_q = g_sf * bl.sd[m - 1];
          g_sd[m - 1] += g_sf * qv;
          g_rho[block_rho] += g_q * (0.5 / (D * qv));
          const double g_D = g_q * (-0.5 * qv / D);
          for (int c = 0; c < m; ++c) g_zT_chain[zt0 + c] += g_D * 2.0 * z_T_p[zt0 + c];
          ++block_rho;
        }

        // sd[j] = sqrt(pi[j] * trace)
        for (int j = 0; j < nc; ++j) {
          g_pi[j]  += g_sd[j] * 0.5 * bl.trace / bl.sd[j];
          g_trace  += g_sd[j] * 0.5 * bl.pi[j] / bl.sd[j];
        }
        // pi[j] = zeta_block[j] / Zsum
        double dotgp = 0.0;
        for (int j = 0; j < nc; ++j) dotgp += g_pi[j] * bl.pi[j];
        for (int j = 0; j < nc; ++j) g_zeta[bl.zeta_off + j] += (g_pi[j] - dotgp) / bl.Zsum;
        // trace = nc * s^2
        g_s = g_trace * 2.0 * nc * bl.s;
      }

      // s = tau_i * scale_i * dispersion
      g_tau_chain[bl.tau_idx] += g_s * re_scale[bl.tau_idx] * dispersion;
      g_disp += g_s * tau[bl.tau_idx] * re_scale[bl.tau_idx];
    }

    // --- priors + Jacobians + constrain-backprop, writing the free-space grad ---

    // z_beta
    for (int k = 0; k < K; ++k) {
      double g = g_beta[k] * dbeta_dz[k];
      if (prior_dist == 1 || prior_dist == 2) {
        g += -z_beta_p[k];
        logp += -0.5 * z_beta_p[k] * z_beta_p[k];
      }
      grad[off_z_beta + k] = g;
    }

    // z_b ~ N(0,1)
    for (int s = 0; s < q; ++s) {
      grad[off_z_b + s] = g_z_b[s] - z_b_p[s];
      logp += -0.5 * z_b_p[s] * z_b_p[s];
    }

    // z_T ~ N(0,1) over ALL allocated entries (chain nonzero only on consumed)
    for (int e = 0; e < len_z_T; ++e) {
      grad[off_z_T + e] = g_zT_chain[e] - z_T_p[e];
      logp += -0.5 * z_T_p[e] * z_T_p[e];
    }

    // rho ~ Beta(a1, a2) (onion) + logit Jacobian
    for (int j = 0; j < len_rho; ++j) {
      const double rj = rho[j];
      const double prior_d = (rho_a1[j] - 1.0) / rj - (rho_a2[j] - 1.0) / (1.0 - rj);
      grad[off_rho + j] = (g_rho[j] + prior_d) * rj * (1.0 - rj) + (1.0 - 2.0 * rj);
      logp += (rho_a1[j] - 1.0) * std::log(rj) + (rho_a2[j] - 1.0) * std::log(1.0 - rj)
              + std::log(rj) + std::log(1.0 - rj);
    }

    // zeta ~ Gamma(delta, 1) + log Jacobian
    for (int j = 0; j < len_concentration; ++j) {
      const double zj = zeta[j];
      const double prior_d = (delta[j] - 1.0) / zj - 1.0;
      grad[off_zeta + j] = (g_zeta[j] + prior_d) * zj + 1.0;
      logp += (delta[j] - 1.0) * std::log(zj) - zj + std::log(zj);
    }

    // tau ~ Gamma(shape, 1) + log Jacobian
    for (int i = 0; i < t; ++i) {
      const double ti = tau[i];
      const double prior_d = (tau_shape[i] - 1.0) / ti - 1.0;
      grad[off_tau + i] = (g_tau_chain[i] + prior_d) * ti + 1.0;
      logp += (tau_shape[i] - 1.0) * std::log(ti) - ti + std::log(ti);
    }

    // aux_unscaled (continuous): likelihood-direct + dispersion coupling +
    // prior on aux_unscaled + log Jacobian
    if (!is_binary) {
      double g_aux = (-N / aux + S / (aux * aux * aux)) + g_disp;
      double g_au = g_aux * daux_dau;
      if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0.0) {
        if (prior_dist_for_aux == 1) { g_au += -au; logp += -0.5 * au * au; }
        else if (prior_dist_for_aux == 2) {
          const double nu = prior_df_for_aux;
          g_au += -(nu + 1.0) * au / (nu + au * au);
          logp += -0.5 * (nu + 1.0) * std::log(1.0 + au * au / nu);
        } else { g_au += -1.0; logp += -au; }
      }
      grad[off_aux] = g_au * au + 1.0;  // *au chain (au = exp) + log Jacobian
      logp += par[off_aux];             // + log(au) Jacobian
    }

    // Non-finite values (exp overflow in tau/zeta/aux at extreme leapfrog
    // positions) must THROW, matching Stan's log_prob contract: WALNUTS'
    // NoExceptLogpGrad catches and substitutes -inf logp + zero gradient,
    // treating the step as divergent. Silently returning NaN instead poisons
    // the accept probability and, through Adam, the adapted step size.
    if (!std::isfinite(logp) || !grad.allFinite())
      throw std::domain_error("parametric_model: non-finite log density or gradient");
  }
};

}  // namespace stan4bart

#endif  // STAN4BART_PARAMETRIC_MODEL_HPP
