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

#include <cmath>      // std::exp, std::log, std::sqrt
#include <stdexcept>  // std::invalid_argument
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
///        CSR Z, prior constants, block geometry) are held by value so C2 can
///        refresh y_ / offset_ in place between Gibbs transitions. Matches the
///        walnuts::LogpGrad concept: void operator()(const VectorXd&, double&,
///        VectorXd&) const.
struct ParametricModel {
  // --- dimensions ---
  int N = 0, K = 0, q = 0, t = 0;
  bool is_binary = false;
  bool has_weights = false;

  // --- data (y_ and offset_ mutable for the Gibbs refresh in C2) ---
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

  /// \brief Evaluate the log-posterior and its gradient at the unconstrained
  ///        point `par`. Parameter-independent constants are dropped (so the
  ///        value matches Stan's log_prob up to a constant additive offset).
  void operator()(const Eigen::VectorXd& par, double& logp, Eigen::VectorXd& grad) const {
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
    Eigen::VectorXd rho(len_rho), zeta(len_concentration), tau(t);
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
    Eigen::VectorXd beta(K), dbeta_dz(K);
    for (int k = 0; k < K; ++k) {
      if (prior_dist == 1) { beta[k] = z_beta_p[k] * prior_scale[k] + prior_mean[k]; dbeta_dz[k] = prior_scale[k]; }
      else if (prior_dist == 2) { beta[k] = CFt(z_beta_p[k], prior_df[k]) * prior_scale[k] + prior_mean[k]; dbeta_dz[k] = dCFt(z_beta_p[k], prior_df[k]) * prior_scale[k]; }
      else { beta[k] = z_beta_p[k]; dbeta_dz[k] = 1.0; }
    }

    // --- forward make_theta_L: build lower-tri Cholesky factor T per block,
    //     caching intermediates for the reverse pass ---
    struct Block {
      int nc, level_count, b_off, tau_idx, rho_off, zeta_off;
      double s, trace, Zsum;
      Eigen::VectorXd pi, sd;
      Eigen::MatrixXd T;
      std::vector<double> sf, Dm;      // per onion row m (index m in 2..nc-1)
      std::vector<int> zT_row_start;   // per onion row m: start index into z_T
    };
    std::vector<Block> blocks;
    blocks.reserve(t);

    int b_mark = 0, rho_mark = 0, zeta_mark = 0, zT_mark = 0;
    for (int i = 0; i < t; ++i) {
      Block bl;
      bl.nc = p[i];
      bl.level_count = l[i];
      bl.b_off = b_mark;
      bl.tau_idx = i;
      bl.rho_off = rho_mark;
      bl.zeta_off = zeta_mark;
      const int nc = bl.nc;
      bl.s = tau[i] * re_scale[i] * dispersion;

      if (nc == 1) {
        bl.T = Eigen::MatrixXd::Constant(1, 1, bl.s);
      } else {
        bl.trace = bl.s * bl.s * nc;
        bl.pi = Eigen::VectorXd(nc);
        double zsum = 0.0;
        for (int c = 0; c < nc; ++c) { bl.pi[c] = zeta[zeta_mark + c]; zsum += bl.pi[c]; }
        bl.Zsum = zsum;
        bl.pi /= zsum;
        zeta_mark += nc;

        bl.sd = Eigen::VectorXd(nc);
        for (int c = 0; c < nc; ++c) bl.sd[c] = std::sqrt(bl.pi[c] * bl.trace);

        bl.T = Eigen::MatrixXd::Zero(nc, nc);
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
      blocks.push_back(std::move(bl));
    }

    // --- forward make_b: b_level = T_i * z_b_level ---
    Eigen::VectorXd b(q);
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

    // --- eta = offset_ + X beta + Z b ---
    Eigen::VectorXd eta = offset_;
    if (K > 0) eta.noalias() += X * beta;
    for (int i = 0; i < N; ++i) {
      double acc = 0.0;
      for (int k = Zu[i]; k < Zu[i + 1]; ++k) acc += Zw[k] * b[Zv[k]];
      eta[i] += acc;
    }

    // --- Gaussian likelihood; g_eta = lambda * (w .* r) ---
    const double lambda = 1.0 / (actual_aux * actual_aux);
    Eigen::VectorXd g_eta(N);
    double S = 0.0;  // sum_i w_i r_i^2
    for (int i = 0; i < N; ++i) {
      const double r = y_[i] - eta[i];
      const double w = has_weights ? weights[i] : 1.0;
      S += w * r * r;
      g_eta[i] = lambda * w * r;
    }
    logp += -N * std::log(actual_aux) - 0.5 * lambda * S;

    // --- backprop eta -> beta, b ---
    Eigen::VectorXd g_beta = (K > 0) ? Eigen::VectorXd(X.transpose() * g_eta) : Eigen::VectorXd();
    Eigen::VectorXd g_b = Eigen::VectorXd::Zero(q);
    for (int i = 0; i < N; ++i)
      for (int k = Zu[i]; k < Zu[i + 1]; ++k) g_b[Zv[k]] += Zw[k] * g_eta[i];

    // --- reverse make_b: g_T += g_b_level z_b_level^T (lower tri); g_z_b += T^T g_b_level ---
    Eigen::VectorXd g_z_b = Eigen::VectorXd::Zero(q);
    std::vector<Eigen::MatrixXd> g_T(blocks.size());
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
      const Block& bl = blocks[bi];
      const int nc = bl.nc;
      Eigen::MatrixXd gT = Eigen::MatrixXd::Zero(nc, nc);
      for (int lev = 0; lev < bl.level_count; ++lev) {
        const int base = bl.b_off + lev * nc;
        for (int r = 0; r < nc; ++r) {
          const double gbr = g_b[base + r];
          for (int c = 0; c <= r; ++c) gT(r, c) += gbr * z_b_p[base + c];     // lower tri only
          for (int c = 0; c < nc; ++c) g_z_b[base + c] += bl.T(r, c) * gbr;   // T^T g_b
        }
      }
      g_T[bi] = std::move(gT);
    }

    // --- reverse make_theta_L: g_T -> {dispersion, tau, zeta, rho, z_T} ---
    Eigen::VectorXd g_rho = Eigen::VectorXd::Zero(len_rho);
    Eigen::VectorXd g_zeta = Eigen::VectorXd::Zero(len_concentration);
    Eigen::VectorXd g_tau_chain = Eigen::VectorXd::Zero(t);
    Eigen::VectorXd g_zT_chain = Eigen::VectorXd::Zero(len_z_T);
    double g_disp = 0.0;
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
      const Block& bl = blocks[bi];
      const int nc = bl.nc;
      const Eigen::MatrixXd& gT = g_T[bi];
      double g_s = 0.0;

      if (nc == 1) {
        g_s = gT(0, 0);
      } else {
        double g_trace = 0.0;
        Eigen::VectorXd g_sd = Eigen::VectorXd::Zero(nc);
        Eigen::VectorXd g_pi = Eigen::VectorXd::Zero(nc);

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
  }
};

}  // namespace stan4bart

#endif  // STAN4BART_PARAMETRIC_MODEL_HPP
