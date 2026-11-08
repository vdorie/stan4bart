# shrinkage-priors

Status: planned, not built. Branch shrinkage-priors from main 98db852. Ruling: dbarts dec-B645 (VD 2026-11-07,
"Use your recommendation for this item"), building on dec-B644 (the WALNUTS swap). TODO item shrinkage-priors.

agent: opus (an engine slice: five new prior families in the WALNUTS target, src/parametric_model.hpp, with their
  latent scales and a hand-derived reverse-mode gradient; an R slice that lifts the refusal; a gate-only slice that
  records rstanarm reference posteriors first. The BART side and dbarts.h are untouched.)

rng/oracle: draws for normal, student_t, cauchy and no prior must be BITWISE UNCHANGED: the new state segments are
  empty for them and no arithmetic on their path changes, so every existing fit is `identical()` to main's. The five
  new families have no bitwise oracle; their gates are (1) finite differences of the gradient, (2) a value-level
  oracle, an R translation of rstanarm's model block evaluated at the same unconstrained points (the retired Stan
  log_prob cross-check's job), (3) a unit test of the parametric mean and sigma read back from the stored draw under
  both storage modes (where 0.0-13 failed), and (4) a pre-registered comparison against rstanarm's posterior on an
  exactly shared model.

budget: about 250-350 engine lines (parametric_model.hpp: transforms, adjoint, priors, Jacobians, emission and
  names; the builder; the test export), about 25 R lines net (the refusal out, list validation in), about 300 test
  lines (test-17 fixtures and the value oracle, one new test file), about 250 benchmark lines (not shipped) plus one
  recorded reference file, help and NEWS. One engine slice, as dec-B645 priced it, plus a gate-only slice and an R
  slice.

window: none on the dbarts side. User-facing: `stan_args$prior` accepts hs, hs_plus, laplace, lasso and
  product_normal again; NEWS says what 0.0-13 did with each.

## Goal

`stan4bart(..., stan_args = list(prior = <family>()))` takes every family rstanarm's coefficient prior takes: normal,
student_t, cauchy, hs, hs_plus, laplace, lasso and product_normal. The five shrinkage families are drawn by WALNUTS
with rstanarm's latent scales on the unconstrained scale and a hand-derived gradient, and each matches rstanarm's
posterior on a shared model within Monte Carlo error. Fits under the three families already built draw exactly as
on main.

## Context

- The target: `ParametricModel` (src/parametric_model.hpp). Free vector today, in continuous.stan's declaration order:
  `[z_beta(K), z_b(q), z_T, rho_free, zeta_free, tau_free, aux_unscaled_free(!binary)]`. `ParametricModel::eval`
  builds beta from z_beta (affine for normal, the Cornish-Fisher `CFt` for student_t/cauchy), backpropagates
  `g_beta = X^T g_eta`, then adds priors and Jacobians per segment. `finalize()` throws on prior_dist outside
  {0, 1, 2}. The raw prefix lengths feed `rawConstrainedDim`, `betaConstrainedOffset` and `auxConstrainedOffset`,
  which `parametricMean` and `getAux` use to read back the draw the Gibbs bridge hands to BART.
- The builder: `buildParametricModel` (src/walnuts_sampler.cpp) reads the marshaled list; the test export
  `stan4bart_logdensity_grad` (src/logdensity_export.cpp) returns value and gradient at a point.
- The sampler: `WalnutsSampler` draws the initial position uniform(-init_r, init_r) per coordinate (init_r default 2,
  Stan's own default random start, which rstanarm uses) and seeds the mass from the gradient there.
  `WalnutsSampler::ridgeMove` slices each random-effect block's log tau along the curve that holds b fixed; its
  stationarity argument lists the terms of the log density that involve tau_i or z_b(i).
- The R surface: `stan4bart_fit` (R/stan4bart_fit.R) calls `refuse_shrinkage_prior` before `handle_glm_prior`
  (R/rstanarm_functions.R), which already maps all eight families to codes 1-7 and fills `global_prior_df`,
  `global_prior_scale`, `slab_df`, `slab_scale`; `data.stan` already carries them and `num_normals`. Coefficient prior
  scales use `link = "identity"` (stan4bart's latent-variable reading), not the family's link. stan4bart exports no
  prior constructors: a user's `hs()` is rstanarm's.
- Names consumers: `stan4bart` (R/stan4bart.R) splits `upar` by a pattern that already lists global, local, caux,
  mix and one_over_lambda; `as.array.stan4bartFit` reads only gamma, beta, b and aux rows.
- Storage policy (docs/design/walnuts.md, "Sample storage policy"; docs/plans/sample-storage.md): by default only
  the transformed block (aux, beta, b, theta_L) is stored; raw rows under `save_raw_parameters = TRUE`.
- History: docs/plans/walnuts.md Q(b) and docs/design/walnuts.md "The prior-scope restriction (Q(b))" record the
  refusal this item lifts.
- What 0.0-13 did (dec-B645, run 2026-11-07 on 0.0-13): hs and hs_plus ran and gave wrong answers, lasso and laplace
  aborted the R session, product_normal matched rstanarm. Read in the 0.0-13 source (git show 0.0-13), not run: its
  Stan model's target was rstanarm's, but the generated `get_parametric_mean` skipped `2 + K` rows for the hs
  latents where the draw holds `h + h*K + 1` (2K + 3 for hs, 4K + 3 for hs_plus), so BART was handed the mean
  computed from the wrong rows. That mechanism is inferred from the source; NEWS does not state it.
- Sources read for this plan: rstanarm 2.32.2 from CRAN (inst/stan/parameters/parameters_glm.stan,
  tparameters/tparameters_glm.stan, model/priors_glm.stan, functions/common_functions.stan `hs_prior` and
  `hsplus_prior`, tdata/tdata_glm.stan, data/hyperparameters.stan; R/priors.R, R/data_block.R `handle_glm_prior`,
  R/stan_glm.fit.R, R/misc.R `set_prior_scale`); stan4bart 0.0-13's src/stan_files/continuous.stan and
  continuous.hpp. Papers in dbarts scratch/papers: Piironen and Vehtari (2017, arXiv v1 read; EJS 11(2),
  doi:10.1214/17-EJS1337SI), Carvalho, Polson and Scott (2009, AISTATS, PMLR 5), Bhadra, Datta, Polson and Willard
  (2017, arXiv v2 read; Bayesian Analysis 12(4), doi:10.1214/16-BA1028), Park and Casella (2008, JASA 103(482),
  doi:10.1198/016214508000000337), Betancourt and Girolami (arXiv 1312.0906), Bou-Rabee, Carpenter, Kleppe and Liu
  (2025, arXiv 2506.18746, WALNUTS), Vehtari et al. (2021, R-hat, ESS, MCSE), Magnusson et al. (2024, posteriordb).
  Not obtained: Carvalho, Polson and Scott (2010, Biometrika 97(2), doi:10.1093/biomet/asq017), its open host down;
  the 2009 version states the same horseshoe.

## 1. The densities rstanarm states, and the argument map

rstanarm's data codes (data/data_glm.stan): "0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus, 5 = laplace,
6 = lasso, 7 = product_normal". tdata_glm.stan sets `hs = 2` for code 3, `hs = 4` for code 4, else 0. Parameters
(parameters_glm.stan), all positive ones declared `<lower=0>`:

    vector[prior_dist == 7 ? sum(num_normals) : K] z_beta;
    array[hs] real<lower=0> global;
    array[hs] vector<lower=0>[K] local;
    array[hs > 0] real<lower=0> caux;
    array[prior_dist == 5 || prior_dist == 6] vector<lower=0>[K] mix;
    array[prior_dist == 6] real<lower=0> one_over_lambda;

### hs (code 3)

Constructor (R/priors.R): `hs(df = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)`, returning
`dist = "hs", df, location = 0, scale = 1, global_df, global_scale, slab_df, slab_scale`. `handle_glm_prior` copies
`global_scale`, `global_df`, `slab_df`, `slab_scale`; it does not call `set_prior_scale`, and `prior$autoscale` is
absent, so nothing is autoscaled. Transform (tparameters_glm.stan, common_functions.stan):

    real c2 = square(slab_scale) * caux[1];
    if (is_continuous == 1 && family == 1)
      beta = hs_prior(z_beta, global, local, global_prior_scale, aux, c2);
    else beta = hs_prior(z_beta, global, local, global_prior_scale, 1, c2);

    vector[K] lambda = local[1] .* sqrt(local[2]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt( c2 * lambda2 ./ (c2 + square(tau) * lambda2) );
    return z_beta .* lambda_tilde * tau;

Prior (priors_glm.stan):

    target += normal_lpdf(z_beta | 0, 1);
    target += normal_lpdf(local[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(local[2] | 0.5 * prior_df, 0.5 * prior_df);
    target += normal_lpdf(global[1] | 0, 1) - log_half;
    target += inv_gamma_lpdf(global[2] | 0.5 * global_prior_df, 0.5 * global_prior_df);
    target += inv_gamma_lpdf(caux | 0.5 * slab_df, 0.5 * slab_df);

This is Piironen and Vehtari's (2017) regularized horseshoe, their eq. (2.8): "beta_j | lambda_j, tau, c ~ N(0,
tau^2 lambdatilde_j^2), lambdatilde_j^2 = c^2 lambda_j^2 / (c^2 + tau^2 lambda_j^2)", with their eq. (2.11) slab
"c^2 ~ Inv-Gamma(alpha, beta), alpha = nu/2, beta = nu s^2/2, which translates to a Student-t_nu(0, s^2) slab", and
exactly their Appendix C.2 code ("This parameterization was proposed by Peltola et al. (2014)"): `lambda = aux1_local
.* sqrt(aux2_local)`, `tau = aux1_global * sqrt(aux2_global) * scale_global * sigma`, `aux1 ~ normal(0, 1)`, `aux2 ~
inv_gamma(0.5*nu, 0.5*nu)`. With df = 1 the local scales are half-Cauchy, the horseshoe of Carvalho, Polson and
Scott (2009): "(beta_i | lambda_i, tau) ~ N(0, lambda_i^2 tau^2), lambda_i ~ C+(0, 1)". Inferred, not stated by
rstanarm: a half-normal times the square root of an Inv-Gamma(nu/2, nu/2) is a half-t_nu, so lambda_k ~ half-t_df(0,
1) and tau ~ half-t_global_df(0, global_scale * sigma) for a Gaussian response, global_scale alone for binary.

Map to stan4bart's list (`handle_glm_prior` output, already in `data.stan`): df -> `prior_df` (length K),
global_df -> `global_prior_df`, global_scale -> `global_prior_scale`, slab_df -> `slab_df`, slab_scale ->
`slab_scale`; `prior_scale` (1) and `prior_mean` (0) unused. error_scale: stan4bart's aux (the residual sd) for a
Gaussian response, 1 for binary, rstanarm's rule.

### hs_plus (code 4)

Constructor: `hs_plus(df1 = 1, df2 = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)`, returning
`df = df1, location = 0, scale = df2` and the rest as hs ("scale gets used as a second df hyperparameter").
Transform:

    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] eta = local[3] .* sqrt(local[4]);
    real tau = global[1] * sqrt(global[2]) * global_prior_scale * error_scale;
    vector[K] lambda_eta2 = square(lambda .* eta);
    vector[K] lambda_tilde = sqrt( c2 * lambda_eta2 ./ ( c2 + square(tau) * lambda_eta2) );
    return z_beta .* lambda_tilde * tau;

Prior: as hs, plus

    target += normal_lpdf(local[3] | 0, 1) - log_half;
    // unorthodox useage of prior_scale as another df hyperparameter
    target += inv_gamma_lpdf(local[4] | 0.5 * prior_scale, 0.5 * prior_scale);

Bhadra et al. (2017) eq. (3.3): "(theta_i | lambda_i, eta_i, tau) ~ N(0, lambda_i^2), (lambda_i | eta_i, tau) ~
C+(0, tau eta_i), eta_i ~ C+(0, 1)". rstanarm's is that product at df1 = df2 = 1, with the slab of hs added; Piironen
and Vehtari (2017) only "expect our idea to be useful with the horseshoe+", so the regularized horseshoe+ is
rstanarm's construction, not a paper's. Map: df1 -> `prior_df`, df2 -> `prior_scale` (length K), the rest as hs.

### laplace (code 5)

Constructor: `laplace(location = 0, scale = NULL, autoscale = FALSE)`; `handle_glm_prior` sets the scale by
`set_prior_scale` (default 2.5; stan4bart passes `link = "identity"`, so no probit factor), and `autoscale = TRUE`
multiplies by sd(y) (Gaussian) and divides by sd(x) when QR is off. Transform and prior:

    beta = prior_mean + prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);

Help: "The Laplace distribution is also known as the double-exponential distribution ... This distribution can be
motivated as a scale mixture of normal distributions". Inferred: with V ~ Exp(1), sqrt(2V) Z has density
exp(-|x|)/2, so beta_k ~ Laplace(location_k, scale_k), the mixture Park and Casella (2008) write as their eq. (4).
Map: location -> `prior_mean`, scale -> `prior_scale`, autoscale -> `prior_autoscale`.

### lasso (code 6)

Constructor: `lasso(df = 1, location = 0, scale = NULL, autoscale = FALSE)`; scale as laplace. Transform and prior:

    beta = prior_mean + one_over_lambda[1] * prior_scale .* sqrt(2 * mix[1]) .* z_beta;
    target += normal_lpdf(z_beta | 0, 1);
    target += exponential_lpdf(mix[1] | 1);
    target += chi_square_lpdf(one_over_lambda[1] | prior_df[1]);

Help: "a more Bayesian approach would be to place a prior on 'it', or rather its reciprocal in our case (i.e. smaller
values correspond to more shrinkage toward the prior location vector). We use a chi-square prior with degrees of
freedom equal to that specified in the call to lasso or, by default, 1." So beta_k | 1/lambda ~ Laplace(location_k,
scale_k / lambda), 1/lambda ~ chi-square(df[1]), one shared across coefficients, only df's first element read. This
is not Park and Casella's Bayesian lasso, whose prior is conditional on sigma ("pi(beta | sigma^2) = prod lambda /
(2 sqrt(sigma^2)) exp(-lambda |beta_j| / sqrt(sigma^2))", eq. (3); "Conditioning on sigma^2 is important, because it
guarantees a unimodal full posterior") with "gamma priors on lambda^2 (not lambda)", eq. (6). The ruling asks for
rstanarm's density, so stan4bart draws rstanarm's. Map: df -> `prior_df[1]`, location -> `prior_mean`, scale ->
`prior_scale`.

### product_normal (code 7)

Constructor: `product_normal(df = 2, location = 0, scale = 1)` (`stopifnot(all(df >= 1), all(df ==
as.integer(df)))`); `stan_glm.fit` passes `num_normals = as.integer(prior_df)`, and the Stan data declare
`array[prior_dist == 7 ? K : 0] int<lower=2> num_normals`, so df = 1 passes the constructor and is refused by Stan's
data check. Transform and prior:

    beta[k] = z_beta[z_pos]; ... for (n in 2:num_normals[k]) beta[k] *= z_beta[z_pos]; ...
    beta[k] *= prior_scale[k] ^ num_normals[k];
    beta[k] += prior_mean[k];
    target += normal_lpdf(z_beta | 0, 1);

Help: "The product-normal distribution is the product of at least two independent normal variates each with mean
zero, shifted by the location parameter ... the density of a product-normal variate is symmetric and infinite at
location"; "The prior variance of the regression coefficients is equal to the scale raised to the power of 2 times
the corresponding element of df". rstanarm cites no paper for it. Map: df -> `num_normals` (from `prior_df`),
location -> `prior_mean`, scale -> `prior_scale`. No autoscale argument.

## 2. The parameterization drawn

All five are drawn in rstanarm's parameterization: non-centered in the coefficients (beta = location + z * scale,
z ~ N(0, 1)), the local and global half-t scales split as Piironen and Vehtari's C.2 (a half-normal times the square
root of an inverse gamma), and every positive latent on the log scale, the unconstrained scale Stan uses for
`<lower=0>` (x = exp(u), log Jacobian u). The free vector grows by segments placed after z_beta, in rstanarm's
declaration order:

    [ z_beta(Kz), global(h), local(h*K), caux(h > 0), mix(K, codes 5-6), one_over_lambda(1, code 6),
      z_b(q), z_T, rho_free, zeta_free, tau_free, aux_unscaled_free(!binary) ]

with Kz = sum(num_normals) under product_normal and K otherwise, h = 2 (hs), 4 (hs_plus), 0 otherwise; local holds
local[1] (K values), then local[2], and so on. For codes 0-2 every new segment is empty and the offsets of z_b and
after are unchanged in value, so existing fits draw bit for bit as before. Added dimension: hs 2K + 3, hs_plus
4K + 3, laplace K, lasso K + 1, product_normal sum(num_normals) - K.

Why rstanarm's choice. It is the parameterization rstanarm and Piironen and Vehtari's experiments use; keeping it
keeps 0.0-13's raw row names, makes the value oracle a line-by-line translation of rstanarm's model block, and puts
nothing of ours in the construction. Non-centering is the right default when the data are weak relative to the
prior (Betancourt and Girolami: "When the data are weak these non-centered parameterizations yield simpler
posterior geometries"), which is the case a shrinkage prior is for: the null coefficients. Its known cost is a
strongly identified coefficient, where z and the scales form a ridge (z * s fixed). WALNUTS adapts the step within
an orbit, which lets it reach a funnel's neck where NUTS does not, but Bou-Rabee et al. (2025, Figure 11) report on
Neal's funnel that "both samplers mix very slowly" along the scale axis; the within-orbit adaptation does not
remove the ridge, so it does not change the case for rstanarm's choice. Piironen and Vehtari's simpler C.1 (lambda
and tau as single half-t coordinates, K + 1 fewer coordinates) is the alternative; they "have not observed problems
with the code presented in Appendix C.1" but ran their experiments with C.2. Whether the ridge hurts at stan4bart's sizes is what the strong-signal dataset of the comparison (4.2)
measures; a move along that ridge, the counterpart of `ridgeMove` for the coefficient scales, would be our own
construction and needs a proof by independent agents before build (dec-B591). It is not in this plan (For the
maintainer, item 1).

Ours, stated plainly, none of them a sampling construction (each changes the arithmetic, not the target or the
chain's transition rule):
- The order of the new segments inside the state (above), and local's block order. Draw-relevant only through the
  random start's order of uniform draws.
- The regularized scale computed through a logistic weight (3.1) rather than rstanarm's ratio of sums; equal in
  exact arithmetic, stable at the extremes.
- When K = 0 (no fixed-effect columns) the shrinkage segments are left empty; rstanarm would still carry global and
  caux, drawn from their prior and read by nothing.
- product_normal's gradient through prefix and suffix products, so a factor at exactly zero does not divide.

The ridge move needs no change: its conditional lists the terms involving tau_i or z_b(i), and no new term involves
either (the coefficient latents enter only beta; under hs with a Gaussian response beta also depends on aux, which
the move does not touch). Its offsets are read from `off_tau` and `off_z_b`, which `finalize()` recomputes.

## 3. Log-density and gradient

Notation. g_k = dl/dbeta_k from the likelihood, the existing `g_beta = X^T g_eta`. e = aux for a Gaussian response,
1 for binary. Every positive latent x has unconstrained u, x = exp(u), and its log Jacobian u is added. Constants
(normalizers, rstanarm's `- log_half`) are dropped, as elsewhere in `eval`. Two building blocks, on the log scale
with the Jacobian included:

    half-normal(0, 1):          lp(u) = -exp(2u)/2 + u,          dlp/du = 1 - exp(2u)
    Inv-Gamma(a, a):            lp(u) = -a u - a exp(-u),        dlp/du = -a (1 - exp(-u))
    Exp(1):                     lp(u) = -exp(u) + u,             dlp/du = 1 - exp(u)
    chi-square(nu):             lp(u) = (nu/2) u - exp(u)/2,     dlp/du = nu/2 - exp(u)/2

(the inverse gamma's -(a + 1) u from the density plus u from the Jacobian gives -a u; the chi-square's (nu/2 - 1) u
plus u gives (nu/2) u).

### 3.1 hs and hs_plus

Unconstrained: z_k; local u1_k, u2_k (and u3_k, u4_k for hs_plus); global v1, v2; caux uc. Forward:

    log L_k  = u1_k + u2_k/2                          (hs)
    log L_k  = u1_k + u2_k/2 + u3_k + u4_k/2          (hs_plus: lambda_k * eta_k)
    log tau  = v1 + v2/2 + log(global_prior_scale) + log(e)
    gam      = log c2 = 2 log(slab_scale) + uc
    a_k      = 2 log tau + 2 log L_k                  (log of tau^2 L_k^2)
    w_k      = 1 / (1 + exp(a_k - gam))               (= c2 / (c2 + tau^2 L_k^2), in (0, 1))
    s_k      = exp((a_k + log w_k) / 2)               (= tau * lambda_tilde_k)
    beta_k   = z_k s_k

Check: s_k^2 = tau^2 L_k^2 c2 / (c2 + tau^2 L_k^2) = tau^2 lambda_tilde_k^2, rstanarm's `z_beta .* lambda_tilde *
tau`. Log prior: -sum z_k^2/2 + sum_k [hn(u1_k) + ig(u2_k; df_k/2)] (+ sum_k [hn(u3_k) + ig(u4_k; df2_k/2)]) + hn(v1)
+ ig(v2; global_df/2) + ig(uc; slab_df/2), with hn and ig the half-normal and inverse-gamma lp above.

Derivatives of the scale: d log s_k / d a_k = w_k / 2 and d log s_k / d gam = (1 - w_k) / 2, since d log w_k / d a_k
= -(1 - w_k) and d log w_k / d gam = 1 - w_k. With G_k = g_k beta_k = dl / d log s_k and T = sum_k G_k w_k:

    d/dz_k   = g_k s_k - z_k
    d/du1_k  = G_k w_k                 + 1 - exp(2 u1_k)
    d/du2_k  = G_k w_k / 2             - (df_k/2)(1 - exp(-u2_k))
    d/du3_k  = G_k w_k                 + 1 - exp(2 u3_k)                      (hs_plus)
    d/du4_k  = G_k w_k / 2             - (df2_k/2)(1 - exp(-u4_k))            (hs_plus)
    d/dv1    = T                       + 1 - exp(2 v1)
    d/dv2    = T / 2                   - (global_df/2)(1 - exp(-v2))
    d/duc    = sum_k G_k (1 - w_k) / 2 - (slab_df/2)(1 - exp(-uc))
    Gaussian: dl/daux += T / aux       (d log tau / d aux = 1/aux; added to `g_aux` before the existing
                                         daux_dau and aux_unscaled chain, beside the dispersion term g_disp)

Compute w_k as a logistic of gam - a_k and log w_k as -log1p(exp(a_k - gam)) (or its mirror for a_k < gam), so a
point far into either tail gives a finite w_k and s_k.

### 3.2 laplace and lasso

Unconstrained z_k, um_k (mix), and for lasso uo (one_over_lambda). Forward: s_k = scale_k sqrt(2 exp(um_k)) (times
exp(uo) for lasso), beta_k = location_k + z_k s_k. With G_k = g_k (beta_k - location_k):

    lp       = -sum z_k^2/2 + sum_k [-exp(um_k) + um_k]  (+ (df_1/2) uo - exp(uo)/2 for lasso)
    d/dz_k   = g_k s_k - z_k
    d/dum_k  = G_k / 2 + 1 - exp(um_k)
    d/duo    = sum_k G_k + df_1/2 - exp(uo)/2                                  (lasso)

### 3.3 product_normal

z holds n_k = num_normals[k] factors per coefficient, contiguous. beta_k = location_k + scale_k^n_k * prod_j z_kj;
lp = -sum z^2/2 over all Kz factors; d/dz_kj = g_k scale_k^n_k prod_{i != j} z_ki - z_kj, the product formed from
prefix and suffix products. The sign of a pair of factors is not identified: (z_k1, z_k2) and (-z_k1, -z_k2) give the
same beta_k.

The formulas of 3.1-3.3 were checked at planning, 2026-11-07, by an R finite-difference probe written from them
(hs_plus with the aux term, lasso, product_normal with a factor at exactly 0; three random points each, worst
relative error 3.7e-9). That checks the derivation, not the C++ that S1 writes.

### 3.4 What else changes in eval

- `finalize()` accepts codes 0-7 and computes the new offsets; `dim_` grows as in section 2.
- `rawConstrainedDim`, and with it `betaConstrainedOffset` and `auxConstrainedOffset`, counts the new raw rows in
  their emission order (0.0-13's failure was exactly this count).
- `constrainedParamNames` under save_raw emits, after z_beta, `global.j`, `local.j.k` (j fastest, as 0.0-13's
  generated names), `caux.1`, `mix.1.k`, `one_over_lambda.1`, with the constrained (exp'd) values, as Stan's
  write_array did; z_beta has Kz rows under product_normal.
- The non-finite throw stays the guard for exp overflow in the new latents.

## 4. Storage, reporting, initialization, what a user sees

- Stored: by the existing storage policy, only aux, beta, b and theta_L by default; the latents appear only under
  `stan_args = list(save_raw_parameters = TRUE)`, with the names above. rstanarm likewise keeps only transformed
  quantities (alpha, beta, b, aux, theta_L, mean_PPD: `pars` in `stan_glm.fit`). No new storage call.
- Reported: `extract(fit, "fixef")`, `fixef`, `print` and `summary` show beta as for any prior; `extract(fit,
  "stan")` shows the latents under save_raw; `fit$adaptation$inv_mass` grows with the dimension.
- Initialized: uniform(-init_r, init_r) on every new unconstrained coordinate, as Stan's default random start that
  rstanarm uses and as stan4bart does for its existing coordinates; the gradient mass seed covers them.
- A user sees: `stan_args = list(prior = hs())` (rstanarm's constructor) fits instead of stopping. Help
  (man/stan4bart.Rd, the `prior` item) lists the eight families, says the constructors are rstanarm's, that under hs
  and hs_plus with a continuous response the global scale is multiplied by sigma (rstanarm's rule), that
  `adapt_delta` is the knob rstanarm's help points shrinkage users to, and that product_normal's raw factor rows do
  not converge in R-hat by design while beta does.
- Refused, at the R level before any C++ runs: a hand-built hs or hs_plus list missing a positive finite
  global_scale, global_df, slab_df or slab_scale; product_normal with any df below 2 ("product_normal's 'df' must be
  integers of at least 2", where rstanarm stops at Stan's data check). `refuse_shrinkage_prior` and
  `unsupported_shrinkage_dists` go.

## 5. Checks

### 5.1 Finite differences (test-17-gradient.R; the project's existing rule)

New fixtures, on the existing hand-built data lists (`build_stan_data` gains the hs fields and `num_normals`):
hs and hs_plus with a Gaussian response and an nc = 2 block (both the aux -> dispersion and the aux -> tau paths
accumulate into one aux gradient), the same with a binary response (error scale 1), hs weighted, laplace, lasso,
product_normal with num_normals (2, 3) on K = 2 and (2, 4) with a nonzero location, K = 4 for hs to have several
local scales. Points: the existing four moderate random points per fixture, plus two points with every log-scale
latent drawn at +/-5, so w_k sits near 0 and near 1. Rule as today: max |fd - grad| / max(|grad|, 1) < 1e-6 at h =
1e-5. `model_dim` and `gen_point` learn the new segments; the test asserts the C side's gradient length equals
`model_dim`.

### 5.2 The value oracle (test-17, new; replaces the retired Stan cross-check)

An R function `rstanarm_lp(d, par)` translates rstanarm's model block line by line with R's own densities (dnorm,
an inverse-gamma log density written from its formula, dexp, dchisq, all with log = TRUE) plus the log Jacobians, from
section 1 and not from section 3 or the C++. At each fixture's points, value(C) - rstanarm_lp must be constant: its
spread across the points below 1e-8. FD proves the gradient matches the coded value; this proves the coded value is
rstanarm's density. 1e-8 is ours, chosen against the C1 oracle's landed 9.1e-13 spread.

### 5.3 Read-back under both storage modes (new test file)

The export `stan4bart_logdensity_grad` gains an optional third output: the constrained array, `parametricMean` and
`getAux` at the point, for save_raw FALSE and TRUE. For every family, the test checks the mean equals X beta + Z b
with beta and b recomputed in R from the point, sigma equals aux, and the names equal the expected list. This is the
test that would have caught 0.0-13. The same file fits each family end to end on a small Gaussian and a small binary
model (finite draws, the expected rows under save_raw, `fitted` consistent with `extract`), and checks the refusals.

### 5.4 Existing draws unchanged

`identical()` on fixef, ranef, sigma, Sigma and the stored stan block for fits under normal, student_t, cauchy and
the default prior, Gaussian and binary, save_raw FALSE and TRUE, at two seeds, between the slice's library and main
98db852's. Plus the full suite.

### 5.5 The comparison against rstanarm (benchmarks/R/compare-rstanarm.R)

The shared model. `stan4bart(y ~ x1 + ... + xK + (1 | g) + bart(z), offset = f, offset_type = "bart")` replaces
BART's contribution to the parametric target by the fixed vector f, so the parametric posterior is exactly
`rstanarm::stan_glmer(y ~ 0 + x1 + ... + xK + (1 | g), offset = f)` with the same prior lists, the x columns centered
beforehand (stan4bart centers its design; centering a centered column is a no-op), the same default
`prior_aux = exponential(autoscale = TRUE)` and `decov()`. Probed 2026-11-07 on stan4bart 0.0-14 (n = 150, one
chain of 200 draws): the fixed offset holds (beta 0.479, sigma 0.714 against lm on y - f: 0.471, 0.706), and
rstanarm 2.32.2 fits that formula with hs(). BART still runs and is ignored. A binary response cannot share a model
this way (its latents come from the BART side), so binary is covered by 5.1-5.3 only.

Data, generated from fixed seeds in the script: D1, weak signal: n = 200, K = 5, beta = (1, 0.5, 0, 0, 0), 10 groups
of sd 0.5, sigma 1, f = 2 sin(2 pi z). D2, strong signal (the ridge's regime): n = 1000, K = 5, beta = (1, 1, 0.5, 0,
0), 20 groups. Priors per dataset: rstanarm's defaults of hs(), hs_plus(), laplace(), lasso(), product_normal(), and
normal() as a control the current sampler already passes: 12 cells.

Reference (rstanarm, recorded first, S0): 8 chains of 1000 warmup and 5000 draws, adapt_delta 0.99 for hs and
hs_plus (rstanarm's own advice) and 0.95 otherwise; a cell's reference stands only if it meets posteriordb's
criteria (Magnusson et al. 2024, section 3.3: R-hat below 1.01, no divergent transitions, the E-FMI criterion, and
10,000 approximately independent draws, read here as bulk and tail ESS of at least 10,000 for every compared
quantity, in place of thinning, our substitution). A cell whose reference fails after doubling the draws is recorded
and not compared.

stan4bart: 4 chains of 1000 warmup and 2500 draws. A cell's run is read only if R-hat < 1.01 and bulk and tail ESS
> 400 for every compared quantity (Vehtari et al. 2021: "the rank-normalized ESS is greater than 400"); the run is
sized for ESS of at least 2000 and extended to 4 x 10,000 draws when below. ESS below 400 at that length is a mixing
failure, recorded, not a pass (For the maintainer, item 1).

Compared per cell: each coefficient's posterior mean and 5% and 95% quantiles (the 90% interval), and the same three
for sigma and for the group sd: 21 summaries, 252 over 12 cells. Each gives z = (stan4bart - rstanarm) / sqrt(mcse_s^2
+ mcse_r^2), the MCSEs from `posterior::mcse_mean` and `posterior::mcse_quantile` (Vehtari et al. 2021), which apply
directly here: with the offset fixed, the parametric chain has none of the slow BART alternation that made the
walnuts baselines need their own MCSE. Pass: every |z| below c = qnorm(1 - 0.01 / (2 M)), a Bonferroni bound holding
the familywise false-alarm rate at 1% over the M compared summaries; at M = 252, c = 4.11, fixed in S0 from the
registered M before any shrinkage draw exists. At ESS 2000 against 10,000 a mean shift of about 0.10 posterior sd
is detectable (c sqrt(1/2000 + 1/10000)); 0.0-13's hs error, 0.165 against -0.002, would be caught unless that
coefficient's posterior sd exceeded about 1.6. The comparison of two Monte Carlo estimates through their standard errors under the Markov
chain CLT follows posteriordb's reading; the Bonferroni level, the ESS sizing, the quantiles chosen and the one
allowed rerun below are the project's own reading rule, as walnuts.md's k = 4 was.

A failing summary: the cell is rerun once at a new stan4bart seed and both results are recorded; one that fails again
is a failure, localized through 5.1-5.3. The rule is not loosened after the numbers are seen.

product_normal's mixing: R-hat and ESS are read on beta, never on the factors, whose signs are not identified; the
report adds, per chain and coefficient, the R-hat of |z_kj| and the fraction of draws in each sign orthant, as
information. A null coefficient's chains should cross orthants; a strong one's need not, and its beta is correct
either way because the posterior is symmetric under the joint sign flip. Reading beta only is ours.

dec-B645's informal check (a bart() term against stan_lmer, the setting where 0.0-13's hs gave 0.165) is rerun once
on the full model and reported at landing, not gated: the two are not the same posterior.

## 6. Slices, gates, review

Tiers from docs/plans/README.md, "Process by risk" (dbarts). This plan takes a blind critique before S1.

S0. The reference and the harness, before the code. benchmarks/R/compare-rstanarm.R (data, fits, summaries, the
  rule with M and c written into the output), the rstanarm references recorded to
  benchmarks/baselines/rstanarm-reference-<commit>.rds with a MANIFEST row, and the harness proven on the normal()
  control cells against the current sampler (21 of 21 summaries per control cell within c, at M as registered). rstanarm and
  posterior are benchmark requirements, not Suggests. Gate-only slice: the R-surface tier's one sonnet review.

S1. The engine. ParametricModel (section 3), the builder, the export's read-back output, test-17's fixtures and value
  oracle (5.1, 5.2), the read-back unit tests (5.3, the part that needs no fit), and the header comment on the
  shrinkage families. The R refusal stays, so stan4bart() cannot reach the new families yet. Gates: 5.1 and 5.2 on
  every fixture; 5.4 bitwise on every listed fit; the full suite; sanitizers on the R-loaded path (the export reaches
  every new term; stan4bart's sanitizers workflow and a local ASan run of test-17 and the new file); R CMD check
  --as-cran. Engine tier: one opus review that runs mutants, at least: drop the 1/2 on u2_k, drop the aux -> tau term,
  use df + 1 in an inverse gamma, drop one log Jacobian, write the chi-square term as df_1 uo in place of
  (df_1/2) uo, take the product's derivative by division, count the raw rows as 0.0-13 did; each must fail 5.1, 5.2
  or 5.3.

S2. The R surface and the comparison. The refusal and its list go; the list validation (section 4); help; NEWS
  (section 7); the end-to-end part of 5.3. Gate: 5.5 on all 12 cells, run by the implementer and rerun independently
  by the reviewer from the same recorded reference; the full suite; R CMD check --as-cran. It makes the new posteriors
  reachable, so it takes the engine tier: the S1 reviewer, unless S1's fix round changed engine behavior.

Records at landing: this file's Landing note; docs/design/walnuts.md, "The prior-scope restriction (Q(b))" and its
follow-up, rewritten to what shipped; walnuts.md Q(b) pointed here; TODO's item removed; the dbarts ledger entry for
the calls below. The branch merges once, after S2 is green in CI on its tip.

## 7. NEWS

The 0.0-14 item on the refusal goes (it never reached a release) and this item takes its place, under Sampler:

* The shrinkage coefficient priors `hs`, `hs_plus`, `laplace`, `lasso` and `product_normal` (rstanarm's constructors)
  are accepted again for `stan_args$prior`, drawn by the new sampler with their latent scales, and each checked
  against rstanarm's posterior on a shared model. In 0.0-13 and 0.0-13.1:
  - `hs` and `hs_plus` ran and gave wrong answers: a coefficient whose posterior mean under rstanarm was -0.002
    (hs) or 0.001 (hs_plus) came out at 0.165 and 0.252, and under `hs_plus` a coefficient simulated at 0.5 at
    1.07.
  - `laplace` and `lasso` aborted the R session with an uncaught C++ exception.
  - `product_normal` gave rstanarm's answers; its draws now differ, as every fit's do.

0.0-13 was run (dec-B645); 0.0-13.1 carries the same Stan model and accessors and is stated from its source. The
numbers are the dec-B645 run's; S2 rechecks them against its record before NEWS is committed.

## Calls made

1. rstanarm's parameterization (non-centered, half-t scales split as Piironen and Vehtari's C.2, positives on the log
   scale). Alternative: their C.1, single half-t coordinates, K + 1 fewer dimensions for hs; or centered beta for
   strong data. Neither is rstanarm's; C.1 changes the raw rows and the value oracle's translation.
2. New segments after z_beta, in rstanarm's declaration order. Alternative: appended after aux_unscaled; equally
   draw-neutral for existing fits, but the raw emission would then follow a different order from the state.
3. The scale through a logistic weight (3.1): the gradient needs w_k anyway, and it stays finite where tau^2 L^2
   overflows. Alternative: rstanarm's sqrt(c2 L^2 / (c2 + tau^2 L^2)), the same value wherever both are finite.
4. Empty shrinkage segments at K = 0. Alternative: rstanarm's prior-only global and caux.
5. Initialization uniform(-init_r, init_r) on the new coordinates, Stan's and rstanarm's default. Alternative: a
   start at the prior median of each scale; ours, and not rstanarm's.
6. The value oracle in R (5.2). Alternative: FD and the comparison only, leaving no point-level check of the density.
7. The comparison through `offset_type = "bart"`, an exactly shared posterior. Alternative: dec-B645's bart() term
   against stan_lmer, which compares different posteriors.
8. Binary not compared to rstanarm; covered by FD, the value oracle and the read-back tests, the binary-specific term
   being the error scale 1. Alternative: an approximate comparison against stan_glmer(binomial(link = "probit")).
9. The Bonferroni rule at 1% familywise, z through Vehtari et al.'s MCSEs, one rerun allowed. Alternative:
   walnuts.md's fixed k = 4 with its sd-ratio band.
10. product_normal read on beta, the factors as information. Alternative: gate a sign-invariant function of the
    factors too.
11. Refusals for malformed hand-built lists and for product_normal df < 2, in our words. Alternative: let the C++
    builder throw, which aborts the session as 0.0-13's lasso did.
12. NEWS without the inferred mechanism of 0.0-13's hs failure. Alternative: name it ("BART was handed the
    parametric mean read from the wrong rows"), once S1's read-back test reproduces it on 0.0-13's accessor
    arithmetic.

## For the maintainer

Not settled by dec-B645, met by a user; options and costs, not decided here.

1. A shrinkage family that fails to mix under WALNUTS at the stated sizes (ESS below 400, or a comparison failure
   traced to mixing, most likely hs or hs_plus on D2's strong coefficients). (a) Ship it and say in help what
   rstanarm's help says, more iterations and a higher `adapt_delta`; no cost beyond the words, and a user meets slow
   chains with no warning (item 2). (b) A move along the ridge that holds beta fixed while the coefficient scales
   change, the counterpart of `ridgeMove`: our own construction, so a proof by independent agents before build
   (dec-B591), then about 80 engine lines, its invariance and FD tests, and a fresh engine review. (c) Refuse that
   family again, NEWS saying why. Recommendation, if it happens: (b) for hs and hs_plus, the families users are most
   likely to pick, after the proof.
2. Sampler diagnostics. rstanarm warns of divergent transitions, which is how its users learn a horseshoe fit is in
   trouble; WALNUTS reports none and `check_sampler_diagnostics` is a no-op, so a stan4bart user gets no warning.
   (a) As now. (b) After every fit, warn when a coefficient's rank-normalized R-hat exceeds 1.01 or its bulk ESS is
   under 400 (Vehtari et al. 2021), about 40 R lines without a new dependency; a new warning on fits of every family.
   (c) The same check for the shrinkage families only, a rule that reads differently across sibling calls.
3. A binary fit's default laplace and lasso scale. stan4bart sets coefficient prior scales with `link = "identity"`;
   rstanarm uses the family's link and, for probit, multiplies the default 2.5 by dnorm(0)/dlogis(0), about 1.6. A
   user moving a probit model from rstanarm meets a different default scale. This already holds for normal,
   student_t and cauchy. (a) Keep. (b) Use the family's link as rstanarm does: every binary fit's default coefficient
   prior moves, a default-prior change that dec-B559 would want a study for.
4. Where the constructors come from. stan4bart exports none, so `hs()` needs rstanarm attached (a large install), or
   a hand-built list. (a) Keep, help naming rstanarm. (b) Export stan4bart's own copies (R/priors.R and
   R/rstanarm_functions.R both hold unexported ones; the second, collated later, wins): masks rstanarm's when both
   are attached, and the copies carry pre-2.21.1 defaults (`autoscale = TRUE` in R/priors.R, against rstanarm's
   current FALSE), which would then reach users. (c) Document the list form.

## Out of scope

- The intercept and aux priors (already non-shrinkage), the covariance prior (decov only), rstanarm's R2 prior
  (stan4bart's `ok_dists` does not list it, nor rstanarm's `stan_glmer`'s coefficient prior).
- Any change to the default prior, to ridge_move, to the warmup or the step-size seed.
- The duplicate constructor copies in R/priors.R and R/rstanarm_functions.R (item 4).

## Landing note

Not landed.
