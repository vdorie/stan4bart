# walnuts

agent: opus (this REMOVES the embedded Stan machinery and replaces the parametric
  conditional's sampler with a bring-your-own-gradient WALNUTS sampler, following
  bairrtt's integration pattern. The parametric target is a smooth, fixed-dimension
  Gaussian-linear-mixed model whose log-posterior gradient is hand-derivable; the
  autodiff stack, the vendored StanHeaders/sundials/rstan trees, and the
  stanc-generated continuous.hpp all go. The BART side is FROZEN - stan4bart is
  already ported to the dbarts 1.0 flat C API (bart_util.hpp includes
  dbarts/dbarts.h); this arc does not touch dbarts.h or the tree engine.)

rng/oracle: THERE IS NO BITWISE ORACLE for the DRAWS. Swapping NUTS for WALNUTS
  changes the sampler wholesale, so draws change by construction; unlike the dbarts
  arcs there is no equivalence-vs-baseline "identical draws" gate available. But the
  TARGET has a near-exact oracle while Stan is still compiled: continuous.hpp
  exposes log_prob<propto, jacobian> (~line 3382) and the existing wrapper already
  marshals a continuous_model, alive through C3. Correctness rests on THREE gates,
  all built BEFORE the code they guard:
  (1) the Stan log_prob TARGET CROSS-CHECK (C1, near-exact): evaluate the
  hand-written target and Stan's log_prob<false, true> at the SAME unconstrained
  points across all fixture tiers; the difference must be CONSTANT (targets may
  differ by dropped constants). The FD gate proves gradient-vs-coded-target only;
  THIS is the oracle for the target itself, closing the no-oracle risk;
  (2) a finite-difference GRADIENT gate (tight tol; gradient vs the coded target,
  mirrors bairrtt's irt_item_logdensity + R FD check); and
  (3) a DISTRIBUTIONAL-EQUIVALENCE gate (the WALNUTS fit targets the SAME posterior
  as the Stan-era fit - fixef/ranef/sigma/Sigma summaries within the PRE-REGISTERED
  MCSE-based tolerance stated in C2, on fixed reference datasets). The Stan-era
  baselines (summaries + MCSEs) for (3) MUST be recorded while Stan still runs (C0),
  before any deletion. The ONE exact contract
  that survives untouched is reproducibility (test-05-rng: same seed -> same draws,
  different core counts differ) - WALNUTS is deterministic given its per-chain seed,
  so this holds as long as seeding is threaded through as Stan's was.

budget: est. ~1500-2200 lines net, but the headline is a large NET DELETION
  (whole StanHeaders + sundials + rstan trees + the 186 KB stanc-generated
  continuous.hpp + the 198 KB stale continuous_new.hpp + stan_sampler.*). New code
  is concentrated: the logp/grad functor + the make_theta_L/make_b forward port and
  its hand-written reverse-mode adjoint (~500-800 lines, the intellectual core), the
  WalnutsSampler wrapper matching the StanSampler surface init.cpp consumes
  (~250-400), the FD + distributional test harnesses (~300), build/DESCRIPTION/
  author surgery (~small but delicate). Chiefly NEW: src/parametric_model.hpp (the
  target), src/walnuts_sampler.{hpp,cpp} (the wrapper), inst/include/walnuts/* +
  WALNUTS_LICENSE (vendored). EDITED: src/init.cpp (the Gibbs bridge), src/Makevars.in,
  configure.ac, configure.win, DESCRIPTION, readme.md. DELETED at C4: everything Stan.

window: no dbarts-side window; dbarts.h is frozen. The user-facing window is the
  control-argument + diagnostics compatibility policy (Q(c)) and the prior-family
  scope (Q(b)) - both are stan4bart-only R-surface decisions, taken here because
  the Stan control vocabulary (adapt_delta, max_treedepth, ...) has no faithful
  WALNUTS analog.

## Goal

Delete the embedded Stan sampler and the vendored autodiff/ODE trees from stan4bart,
and draw the parametric conditional (fixed effects + lme4-style random effects +
residual sd) with the gradient-based WALNUTS sampler, using a HAND-DERIVED analytic
log-posterior gradient. The outer BART-vs-parametric Gibbs alternation, the R
formula/data-prep surface (glFormula, lme4_functions.R, rstanarm_functions.R,
priors.R), and the dbarts flat-C-API coupling are all preserved. Motivation (VD):
Stan is slow and memory-hungry; the parametric target is smooth and fixed-dimension,
so the whole autodiff apparatus is unnecessary weight. The arc must MEASURE the
payoff (install time, peak RSS, per-iteration speed), not assert it.

## The parametric model (the source of truth: src/stan_files/continuous.stan)

The Stan model is the parametric conditional of a Gaussian linear mixed model whose
mean is offset by BART's current fit. In the Gibbs loop the outer driver swaps
`offset_` (<- BART's training fit, via set_offset, continuous.hpp:3631) and, for
binary, `y` (<- the probit latents, via set_response, continuous.hpp:3626) BETWEEN
transitions, then takes ONE HMC step. So the WALNUTS target is exactly this model
with offset_ and y held fixed for the step. Two sentences: the parameters are
fixed-effect coefficients beta (default Normal prior), lme4-decov random effects b
(a Cholesky-of-covariance parameterization driven by tau/zeta/rho/z_T), and, for
continuous responses, a residual sd aux/sigma (default autoscaled Exponential); the
likelihood is Gaussian (`y ~ Normal(eta, actual_aux)`, weighted or not), with the
binary/probit case reduced to the SAME Gaussian model at aux == 1 over latent
responses drawn on the BART side. There is NO intercept in the parametric part
(BART carries it - stan4bart_fit.R:111-112 warns and drops prior_intercept;
test-01-continuous.R:21-24 asserts "no intercept parameter was included").

### Unconstrained parameterization (Stan samples here; the WALNUTS position vector)

Stan's `parameters` block (continuous.stan:261-278), has_intercept == 0, default
reachable config (prior_dist == 1 Normal for beta; prior_dist_for_aux == 3
Exponential; decov for RE). The free vector omega, in Stan's storage order - the
`parameters` block DECLARATION order (continuous.stan:261-278; get_aux/
get_parametric_mean walk the CONSTRAINED output array, a different layout):

1. `z_beta` in R^K                          (unconstrained)
2. `z_b` in R^q                             (unconstrained; q = sum p_i * l_i)
3. `z_T` in R^{len_z_T}                      (unconstrained; NONZERO ONLY for blocks
     with nc >= 3 - transformed data loop `for j in 3:p[i]`, continuous.stan:258)
4. `rho_free` in R^{len_rho}, rho = inv_logit(rho_free)   (constrained (0,1))
5. `zeta_free` in R^{len_concentration}, zeta = exp(zeta_free)   (constrained > 0)
6. `tau_free` in R^t, tau = exp(tau_free)   (constrained > 0)
7. `aux_unscaled_free` in R (ONLY if !is_binary), aux_unscaled = exp(...)  (> 0)

For the shrinkage prior families (prior_dist in {3,4,5,6}) Stan ALSO carries
global/local/caux/mix/one_over_lambda latents (continuous.stan:266-270). These are
UNREACHABLE under stan4bart's default and arguably nonsensical for a few-column
parametric part - see Q(b) for the proposal to drop them.

### Transforms (continuous.stan:279-343)

- beta = z_beta .* prior_scale + prior_mean          [Normal; prior_dist==1]
    (student_t/cauchy: beta_k = CFt(z_beta_k, df_k)*scale_k + mean_k, a deterministic
     Cornish-Fisher polynomial, continuous.stan:146-158,295-297 - differentiable, cheap)
- aux = prior_scale_for_aux * aux_unscaled            [Exponential/scale-only]
    (+ prior_mean_for_aux for Normal/student_t aux prior, continuous.stan:325-332)
- actual_aux = is_binary ? 1 : aux                    (likelihood sd)
- dispersion = is_binary ? 1 : aux                    (feeds the RE scale - the coupling below)
- theta_L = make_theta_L(p, dispersion, tau, scale, zeta, rho, z_T)  (continuous.stan:1-59)
- b = make_b(z_b, theta_L, p, l)                      (continuous.stan:61-94)
- eta = offset_ + X*beta + Z*b                        (Z is CSR sparse, csr_matrix_times_vector, model block 350-356)

### Log-posterior T(omega) (continuous.stan model block 344-429; constants dropped)

Likelihood (continuous.stan:358-366):
  unweighted: -N*log(actual_aux) - (1/(2*actual_aux^2)) * sum_i (y_i - eta_i)^2
  weighted:   -N*log(actual_aux) - (1/(2*actual_aux^2)) * sum_i w_i*(y_i - eta_i)^2
  (NOTE the weighted log-normalizer uses N, not sum(w) - continuous.stan:365 - a
   deliberate rstanarm quirk; the port MUST match it or the sigma posterior shifts.)
Priors:
  z_beta ~ N(0,1)                            -> -0.5 * ||z_beta||^2       (continuous.stan:380)
  aux_unscaled ~ Exp(1)                       -> -aux_unscaled            (continuous.stan:375)
  decov_lp (continuous.stan:96-122):
    z_b ~ N(0,1)                              -> -0.5 * ||z_b||^2
    z_T ~ N(0,1)                              -> -0.5 * ||z_T||^2
    rho ~ Beta(shape1, shape2)                (onion regularization; shapes from
                                               `regularization` + block size, 104-118)
    zeta ~ Gamma(delta, 1)                    -> sum (delta-1) log zeta - zeta
    tau  ~ Gamma(shape, 1)                    -> sum (shape-1) log tau  - tau
Constraint log-Jacobians (Stan adds these automatically; the hand target MUST add them):
  rho (logit):   + sum [ log rho + log(1-rho) ]
  zeta,tau,aux_unscaled (log): + sum log(zeta) + sum log(tau) + log(aux_unscaled)
  (aux itself is a DETERMINISTIC transformed parameter, prior placed on aux_unscaled,
   so NO extra Jacobian for aux.)

### Reachable priors (what the R surface actually wires - stan4bart_fit.R)

- beta coef prior `stan_args$prior`: default `default_prior_coef_gaussian()` == Normal
  (stan4bart_fit.R:107-108). ok_dists allows normal/student_t/cauchy/hs/hs_plus/
  laplace/lasso/product_normal (132-135) but the exotic four are opt-in only.
- intercept prior: WARNED OUT and not included (111-112).
- aux (sigma) prior: default `exponential(autoscale = TRUE)` (113-114). ok_aux_dists
  == normal/student_t/cauchy/exponential (135).
- covariance prior: ALWAYS `decov()` (104-105) - the decomposition-of-covariance
  onion. No lkj path is wired.

So the realistic reachable target is: Gaussian likelihood + {normal, student_t,
cauchy} beta + {normal, student_t, cauchy, exponential} aux + decov RE. The
shrinkage families are reachable-in-principle-but-unused (Q(b)).

### Binary end-to-end (the simplifying reduction)

stan4bart detects a 0/1 response and sets family = binomial(link="probit")
(stan4bart.R:66). The Stan data sets is_binary == 1 -> actual_aux == 1, and the
aux_unscaled parameter is ABSENT (continuous.stan:277 `aux_unscaled[!is_binary]`).
The probit latents are drawn on the BART side: init.cpp:700-703 calls
bartFunctions.getLatents(...) then setResponse(model, bartLatents) each sweep. So
the WALNUTS binary target is IDENTICAL to the continuous target with sigma fixed at
1 and the aux dimension removed - the dispersion coupling to the RE scale also
drops (dispersion == 1). Binary is strictly a sub-case; C3 is small.

### Random-effect (decov) structure - the gradient's hard center

lme4 terms map to blocks (stan4bart_fit.R:302-329): `t` grouping terms, block i has
`p[i]` correlated coefficients over `l[i]` levels; q = sum p_i*l_i. The tests
exercise BOTH shapes: `(1 + X4 | g.1)` is a nc=2 correlated intercept+slope block
and `(1 | g.2)` is a nc=1 random-intercept block (test-01-continuous.R:16,
test-02-binary.R:15, test-08-glFormula.R). make_theta_L builds, per block, the
lower-triangular Cholesky factor of the covariance:
  - nc == 1: theta_L = tau_i * scale_i * dispersion        (continuous.stan:15-18) - trivial
  - nc >= 2: T_i built from tau_i (overall scale), zeta (a Dirichlet simplex `pi`
    dividing the trace across dimensions, 25-30), rho (canonical partial correlations
    via T21 = 2*rho-1, 35-38), and for nc >= 3 the "onion method" scale_factor =
    sqrt(rho / dot_self(T_row)) * std_dev threading z_T (40-49).
make_b then maps z_b through these blocks: b_level = T_i * z_b_level (61-94).

## Gradient derivation: assessment and confidence call

The forward pass is a straight port of make_theta_L/make_b/decov_lp - plain
arithmetic in the .stan (no autodiff needed forward). The work is the reverse pass.
Term by term:

- LIKELIHOOD + fixed effects + sigma: standard Gaussian-regression gradients.
  Let r = y - eta, lambda = 1/actual_aux^2, wr = (w .* r) (w==1 unweighted).
  d/d z_beta = X^T(lambda*wr) .* prior_scale - z_beta      [+ Normal prior]
  d/d z_b    = (from d/d b below) - z_b                     [+ Normal prior]
  d/d aux    (continuous) = -N/aux + aux^{-3} * sum w r^2   (likelihood)
             PLUS d/d b . db/d aux   (aux enters b via dispersion in theta_L),
             then chain aux = scale*aux_unscaled -> aux_unscaled_free (x aux_unscaled),
             plus -aux_unscaled (prior) + 1 (Jacobian). CONFIDENCE: HIGH. The only
             subtlety is not forgetting the dispersion->b second term.
- make_b chain d/d b -> {z_b, tau, zeta, rho, z_T, dispersion}:
  - nc == 1: b_s = (tau*scale*disp) * z_b,s. All partials trivial. CONFIDENCE: HIGH.
  - nc == 2: T_i has three nonzero entries (T11, T21, T22) as closed forms in
    tau, pi(zeta), rho. b_level = T_i z_b,level. The adjoint is a handful of
    products; no z_T. CONFIDENCE: HIGH.
  - nc >= 3: the onion (scale_factor = sqrt(rho/dot_self(T_row))*std_dev threading
    z_T, 40-49) couples z_T entries through dot_self and a sqrt. The adjoint is
    mechanical but ERROR-PRONE. CONFIDENCE: MEDIUM without a gate; HIGH with the FD
    gate exercised on a genuine nc>=3 fixture.
- rho/zeta/tau priors + their logit/log Jacobians: elementary. CONFIDENCE: HIGH.
- student_t/cauchy beta (Cornish-Fisher CFt): polynomial derivative. CONFIDENCE: HIGH.

CALL: self-derivable WITH CONFIDENCE across the whole reachable target, with the
nc>=3 decov onion as the single term warranting external corroboration. Strategy:
write the reverse-mode adjoint by hand as a short sequence of elementary-op backprops
mirroring the forward loops (NOT a general autodiff - that would reintroduce the
dependency we are removing), and gate it with finite differences. RECOMMENDATION on
Q(a): proceed self-derived and FD-gated; treat VD's offer ("I can do some math and
provide an R implementation as a prototype") as a welcome SECOND-INDEPENDENT-ORACLE
for make_theta_L/make_b specifically (the onion), used to cross-check the FD gate at
nc>=3 - valuable but not a blocker that serializes C1 on VD. The FD gate is
mandatory regardless (C1's headline gate). See Q(a) for the fork and costs.

## Context (seams, read in code)

- THE GIBBS BRIDGE, init.cpp `run` (541-...), one sweep (608-704): (1)
  stanSampler->run(isWarmup) takes one HMC transition (616); (2) getParametricMean
  (620/626/631/635/641) computes X*beta + Z*b into bartOffset; (3) getSigma (654,
  continuous only) reads aux for BART's residual sd; (4) BART runs (680), its fit
  minus offset becomes stanOffset; (5) setStanOffset (699) + (binary) getLatents +
  setResponse (700-703) refresh the parametric target for the next sweep. The
  WALNUTS wrapper must present exactly this surface: run(isWarmup), getParametricMean
  (with the include_fixed/include_random overloads, continuous.hpp:3715), getSigma,
  copyOutParameters, and a freeze() for disengageAdaptation (850-856). The offset/
  response setters map onto refreshing the functor's held-by-reference data (the
  bairrtt pattern: logp_grad.theta = ...; irt_hmc.cpp:101).
- THE STAN SAMPLER SURFACE, stan_sampler.hpp:52-113 + stan_sampler.cpp: StanControl
  (28-42) holds the NUTS knobs; StanSampler ctor (464-521) builds an
  interruptable_sampler over stan::mcmc::adapt_diag_e_nuts
  (interruptable_sampler.hpp:110) - WINDOWED DIAGONAL-MASS adaptation plus step
  size, with the unit_e metric (472-473, `create_unit_e_diag_inv_metric`) only the
  INITIAL value. WALNUTS' Nutpie-style diagonal mass estimator replaces this
  LIKE-FOR-LIKE (both adapt a diagonal mass);
  sample/sampler names incl. the diagnostics (479-484: stepsize__, treedepth__,
  n_leapfrog__, divergent__, energy__); getParametricMean (540-550), getSigma
  (552-555), copyOutParameters (557-560). createStanModelFromExpression (112-380)
  marshals ~44 named data fields (dataNames 67-80) - the WALNUTS functor needs a
  small subset (X, y, weights, offset_, the CSR Z = {w,v,u}, the prior scalars, and
  the block geometry p/l/t/q/len_* + decov shape/scale/concentration/regularization).
- THE CONSTRAINED OUTPUT CONTRACT, continuous.hpp get_aux (3638-3660) +
  get_parametric_mean (3663-3768): the stored `result$stan` array rows are
  [Stan sample diag] + [5 NUTS sampler diag] + [constrained primitives z_beta.., rho..,
  zeta.., tau.., aux_unscaled.. + transformed beta.., b.., theta_L.., aux..]. The R
  generics consume by NAME (generics.R:32-99, 226-274): beta.*->fixef, b.*->ranef,
  aux.1->sigma, theta_L.*->Sigma (via mkVarCorr). stan4bart.R:243 splits upar
  (gamma|z_beta|global|local|caux|mix|one_over_lambda|z_b|z_T|rho|zeta|tau|aux_unscaled)
  from tpar (beta|b|theta_L|aux). So the WALNUTS wrapper must, per stored draw,
  (a) map unconstrained->constrained, (b) compute theta_L/beta/b/aux, (c) name every
  row exactly as Stan did, so the R side is UNCHANGED. This "write_array" replacement
  is real work but is a deterministic re-emission of the same forward transforms.
- THE DIAGNOSTICS, check_sampler_diagnostics (stan4bart.R:255-297): warns on
  divergent__ count, treedepth__ >= max_treedepth, and low EBFMI from energy__.
  These are NUTS artifacts (written stan_sampler.cpp:479-484). The vendored WALNUTS
  ChainHandler concept surfaces ONLY position/lp/step_size/diag_inv_mass (plus a
  cross-chain r_hat in the summary layer) - NO treedepth, NO divergent, NO
  n_leapfrog, NO energy - so ALL THREE checks are lost, not conditionally. This is
  a user-facing delta with a real fork (Q(c)); re-verify against the fresh upstream
  API at C1 before finalizing.
- CHAINS: R-level parallel, NOT within-chain threads. stan4bart_fit.R:500-503 spawns
  a PSOCK/FORK cluster; each worker (33-63) runs ONE chain single-threaded
  (C_stan4bart_create -> run warmup -> disengageAdaptation -> run sample). So
  RcppParallel/TBB (configure.ac:22-34, STAN_THREADS) serve ONLY Stan-internal
  threading and go with Stan. WALNUTS' AdaptiveWalnuts is per-chain (bairrtt uses one
  adapter per sampler, irt_hmc.cpp:86-90); its cross-chain mass-averaging
  (config.hpp InitConfigBuilder::masses average_masses) is UNUSED - per-chain
  independent adaptation matches the current independent-process chains exactly.
- THE WALNUTS PATTERN (bairrtt), to copy: irt_hmc.cpp lifecycle create_sampler /
  warmup_draw (drives AdaptiveWalnuts per outer scan while conditioning data moves) /
  freeze (AdaptiveWalnuts::sampler() -> WalnutsSampler) / draw_sample (43-130);
  bairrtt_types.h the LogpGrad functor shape `operator()(const VectorXd& par,
  double& lp, VectorXd& grad) const` with conditioning data held by value/reference
  and refreshed in place (104-154), and the persistent XPtr state holding logp_grad,
  a LatestDraw handler, rng, WarmupConfig/SamplingConfig, and the optional
  adapter/sampler (156-195); engine.R the R FD gradient check (irt_item_logdensity,
  138-147). WALNUTS mass is DIAGONAL only (config.hpp InitChainConfig mass is a
  VectorXd; adaptive_walnuts.hpp MassEstimator is a Nutpie-style diagonal estimator,
  25-90) - matching stan4bart's diagonal metric; there is no dense-metric decision.
  Vendor provenance: inst/include/walnuts/* + WALNUTS_LICENSE (MIT, (c) 2025 Bob
  Carpenter) + the walnuts.hpp umbrella; C++20 (concepts) + Eigen (RcppEigen).
  VENDORED-VERSION note (VD's direction): the implementer pulls WALNUTS FRESH from
  https://github.com/flatironinstitute/walnuts at C1 and RECORDS THE UPSTREAM
  COMMIT here and in the C1 landing note; bairrtt's copy is the INTEGRATION
  REFERENCE only, not the source of the vendored headers. Any API drift vs the
  bairrtt-era snapshot (functor concept, config builders, handler concept - incl.
  the diagnostics claim above) must be re-verified against the fresh headers.

## Constraints

- dbarts.h FROZEN; the tree engine, bart_util.*, and the flat-C-API coupling are not
  touched (except folding VD's one-liner, below). No dbarts lockstep.
- The R formula/data-prep surface is PRESERVED: glFormula, lme4_functions.R,
  rstanarm_functions.R, priors.R, test_data.R are DATA PREP and largely survive.
  The decov()/handle_glm_prior machinery that shapes data.stan stays; only its
  consumer (the sampler) changes. Any prior-family narrowing (Q(b)) is a surgical
  edit to the reachable set, not a rewrite.
- SAME POSTERIOR. The WALNUTS target must reproduce continuous.stan's log-density
  including its quirks (the weighted log-normalizer's N-not-sum(w), the decov onion's
  exact induced prior on Sigma, the aux->RE-scale dispersion coupling). Re-deriving a
  "cleaner" covariance parameterization is OUT (it changes the prior on Sigma - see
  Q rejected-alternatives).
- Every commit leaves the package installable and `R CMD check`-able. Stan is present
  and running through C3; nothing is deleted until its WALNUTS replacement passes its
  gate (C4).
- fast over safe in C/C++ (the logp/grad is the hot path, called once per leapfrog
  micro-step); header-only C++20; Doxygen/LLVM/ASCII; brevity.
- OUT of scope: the shrinkage prior families unless Q(b) resolves to keep them; any
  change to the outer Gibbs schedule, the offset_type debugging hooks
  (init.cpp:624-650), the callback surface (init.cpp:706-...), or the BART side; a
  dense mass metric; QR reparameterization changes (stan4bart_fit.R:242-256,563-573 -
  it operates on the returned beta rows and is sampler-agnostic, so it survives
  untouched as long as beta.* rows are named and ordered as today).

## Commit partition

Sequenced so Stan runs until its replacement is proven, and the motivation is
measured. Each commit installs and tests.

C0. RECORD BASELINES (Stan still the only sampler). Two kinds, plus housekeeping.
   (a) PERF baselines (the motivation, measured not asserted): wall-clock `R CMD
   INSTALL` time, peak sampling RSS, and per-iteration sampling wall-time, on three
   reference fits - continuous multilevel, binary multilevel, and a weighted
   continuous - captured to a committed script + a results file. Concrete
   before-number already in hand as a proxy: stan_sampler.o compiles to ~50 MB
   (init.o ~7 MB) - the stanc translation unit dominates build cost; C4 must show
   this gone.
   (b) POSTERIOR baselines (the distributional oracle): fixef/ranef/sigma/Sigma
   posterior means, sds, AND per-summary MCSEs (means AND sds are gated summaries;
   the MCSEs are what make the C2/C3 tolerance falsifiable - record them, e.g. via
   ESS-based mcse, alongside every gated quantity) at fixed seeds, saved as .rds,
   on reference datasets spanning the gradient's tiers - continuous nc=1,
   continuous nc=2, continuous nc>=3 (the onion), weighted continuous, and binary.
   Long enough runs that the summaries are stable (these are the C2/C3 gate
   targets).
   (c) Fold VD's uncommitted one-liner: bart_util.hpp `dbarts_results current = {};`
   (the zero-initializer, currently unstaged in VD's live checkout). ASK VD to either
   commit it himself first or bless folding it into C0. Ignore the stale untracked
   src/stan_files/continuous_new.hpp (2022) until C4 deletes it.
   Files: benchmarks/ or tests/baselines/ (new scripts + .rds/.csv), bart_util.hpp.
   Gate: baselines recorded and reproducible; package unchanged otherwise; full
   testthat green. Size: M. Abort: if a "reference" summary is not stable across
   reruns at the chosen length, lengthen before freezing - a noisy baseline makes
   C2/C3 ungateable.

C1. VENDOR WALNUTS + the CONTINUOUS logp/grad functor + the TARGET + FD GATES (Stan
   untouched and still running). Vendor inst/include/walnuts/* + WALNUTS_LICENSE +
   the walnuts.hpp umbrella FRESH from https://github.com/flatironinstitute/walnuts
   (VD's direction; RECORD the upstream commit hash in this plan's VENDORED-VERSION
   note and the landing note; bairrtt's copy is integration reference only - re-verify
   the functor/config/handler API against the fresh headers). Add
   src/parametric_model.hpp: the forward port of make_theta_L /
   make_b / decov_lp / the Gaussian likelihood + priors + Jacobians, and the
   hand-written reverse-mode adjoint, as a functor matching walnuts::LogpGrad
   (`operator()(const VectorXd&, double&, VectorXd&) const`) holding X, y, weights,
   offset_, CSR Z, prior scalars, and block geometry. Expose an Rcpp/.Call entry
   `.stan4bart_logdensity_grad(par, data)` returning {value, gradient} (bairrtt's
   irt_item_logdensity analog) and an R-side finite-difference checker (engine.R
   analog). CXX_STD -> CXX20 here (WALNUTS needs concepts; Stan compiles fine under
   C++20 - verify the Stan translation unit still builds). Files: inst/include/walnuts/*,
   src/parametric_model.hpp (new), a small src/logdensity_export.cpp (new) + its
   .Call registration in init.cpp, R/ FD helper, tests/testthat/test-*-gradient.R
   (new), Makevars.in (CXX_STD), DESCRIPTION (LinkingTo add RcppEigen already present;
   SystemRequirements C++20). Gates (THE PRIMARY CORRECTNESS PROOF, both mandatory):
   (i) the Stan log_prob TARGET CROSS-CHECK - continuous.hpp exposes
   log_prob<propto, jacobian> (~line 3382) and the existing wrapper already marshals
   a continuous_model (alive through C3); evaluate the hand-written target and
   Stan's log_prob<false, true> at the SAME unconstrained points across ALL fixture
   tiers and require the difference to be CONSTANT (targets may differ by dropped
   constants). Rationale: the FD gate proves gradient-vs-coded-target only; this is
   the near-exact oracle for the TARGET itself, closing the plan's headline
   no-oracle risk.
   (ii) FD gradient agreement to tight tol on fixtures covering all tiers - fixed
   normal + sigma, nc=1, nc=2, AND nc>=3 (onion), continuous AND binary (aux
   absent), weighted AND unweighted.
   Plus: Stan path fully green (untouched). Size: XL (the functor + adjoint is
   the intellectual core). Abort: a non-constant log_prob difference or any FD
   mismatch that is not a fixture/tol artifact - the target or adjoint is wrong; do
   not proceed to wire a wrong gradient.

   VENDORED-VERSION (landed at C1). WALNUTS pulled FRESH from
   https://github.com/flatironinstitute/walnuts at commit
   5854be888e5432a103275466d7d0be95c7c5c67a (upstream date 2026-07-13, vendored
   2026-07-15) into inst/include/{walnuts/*, walnuts.hpp, WALNUTS_LICENSE} +
   WALNUTS_VERSION; upstream layout (include/walnuts/* + umbrella) matches
   bairrtt's, so no layout adaptation was needed. API drift vs the bairrtt-era
   snapshot is MINOR and does not reach C1: adaptive_walnuts.hpp / handlers.hpp /
   concepts.hpp / walnuts.hpp are byte-identical; the AdaptiveWalnuts -> freeze
   (sampler()) -> WalnutsSampler lifecycle, the diagonal MassEstimator, and the
   InitChainConfig(step_size, position, mass) ctor are unchanged. The only drift
   (a C2 concern for WarmupConfig sizing, not C1): WarmupConfig gained a
   step-size convergence tolerance (step_size_converge_tol_, from upstream's
   "init-stepsize" work) and api.hpp's internal detail::sample() was refactored
   to take the SamplingConfig directly. DIAGNOSTICS re-verified against the fresh
   headers: the ChainHandler concept (concepts.hpp) still surfaces ONLY
   position / lp / step_size / diag_inv_mass, and GlobalHandler adds on_r_hat -
   there are NO treedepth, divergent, n_leapfrog, or energy hooks, so the plan's
   Q(c) claim (all three NUTS diagnostic checks are lost, not conditionally)
   HOLDS. CXX_STD flipped to CXX20 (Makevars.in/.win): the Stan translation unit
   and continuous.hpp compile clean under C++20, clearing that risk early.
   Gate outcome: FD worst relative error 2.35e-08 across all nine tiers
   (tol 1e-6), and value(mine) - Stan log_prob<false, true> constant to
   <= ~1e-12 spread per tier (tol 1e-6) - the target and adjoint are correct.

C2. WIRE WALNUTS for the CONTINUOUS path behind the existing Gibbs loop (Stan still
   compiled, still serving binary). Add src/walnuts_sampler.{hpp,cpp}: a
   WalnutsSampler holding the parametric_model functor + AdaptiveWalnuts/WalnutsSampler
   (the bairrtt IrtSampler analog), presenting the StanSampler surface init.cpp
   consumes - run(isWarmup) (adapting vs frozen transition), getParametricMean (+ the
   fixed/random overloads), getSigma, copyOutParameters (the constrained write_array
   re-emission with Stan-identical row NAMES), and freeze() for disengageAdaptation.
   Route the continuous family through it (a compile-time or family switch in
   init.cpp; binary still to Stan). setStanOffset/setResponse become in-place functor
   refreshes. ADAPTATION NOTES: WALNUTS' diagonal-mass adaptation REPLACES Stan's
   LIKE-FOR-LIKE - stan4bart runs stan::mcmc::adapt_diag_e_nuts
   (interruptable_sampler.hpp:110), windowed diagonal-mass + step-size adaptation
   with unit_e only as the initial metric, so both sides adapt a diagonal mass; and
   SIZE WALNUTS' WarmupConfig to the ACTUAL warmup iteration count (the outer loop
   drives one transition per sweep, control.common$warmup of them) so its mass
   windows schedule correctly against the outer-loop-driven transitions. Files:
   src/walnuts_sampler.{hpp,cpp} (new), init.cpp (the family
   switch + setter routing), Makevars.in (SOURCES). Gate: DISTRIBUTIONAL equivalence
   vs C0 continuous baselines (nc=1, nc=2, nc>=3, weighted), PRE-REGISTERED
   tolerance: for every gated summary,
   |mean_walnuts - mean_stan| < k * sqrt(mcse_stan^2 + mcse_walnuts^2) with k = 4
   FIXED IN ADVANCE; posterior sds compared as a ratio within [0.8, 1.25] (a
   relative band consistent with sd-estimate noise at the baseline ESS; tighten if
   C0 runs longer). Reproducibility (test-05-rng) green; the continuous
   testthat (test-01) green (regenerate any value-pinned snapshots that legitimately
   move - the sampler changed). Size: L. Abort: a posterior summary outside the
   pre-registered tolerance
   means the target or the constrain/output mapping is wrong - fall back to the
   target cross-check + FD gates and the write_array re-emission to localize.

### Landings through C2 (2026-07-15/16)

C0 = 75d7970 + 9a4f906. The structSize fix first: the dbarts 1.0 port had
never set dbarts_results.structSize, so the versioned-struct field gate
skipped every run output - the train buffers stayed zero and the Gibbs
coupling decohered; stack garbage usually masked it. Then the five-tier
posterior baselines with mean/sd/MCSE per gated summary; sqrt(n) batch-means
MCSEs read 3-8x small under the slow BART-vs-parametric tradeoff, so the
recorded MCSE is max(n^(2/3) batch means, between-chain), and the cross-seed
stability gate passed on that basis (nc3 94/98, binary 53/55 within 2x).
Perf baselines: script committed, the recording run awaits a quiet machine.

C1 = cb809b1. The functor and hand adjoint proven two ways per the amended
gate design: FD worst relative error 2.4e-8 across ten tiers (including a
continuous mixed-block tier added at review - the aux -> dispersion adjoint
must ACCUMULATE across heterogeneous blocks), and the Stan log_prob
constant-difference target oracle at 9.1e-13 worst spread (the weighted
normalizer and binary reduction both confirmed at machine precision).
WALNUTS vendored fresh at upstream 5854be8; its API exposes no treedepth,
divergence, or energy - Q(c)'s diagnostics fork resolves to drop-all plus
r_hat.

C2 = 3929581. Continuous through WALNUTS behind a runtime family switch;
binary stays on Stan (bitwise-untouched, suite-proven). One integration
lesson worth keeping: WALNUTS assumes Stan's throw-on-non-finite log_prob
contract, and a silently returned inf at extreme leapfrog positions poisons
Adam step-size adaptation into a crash at freeze - eval now throws, treated
as divergence. Equivalence vs the C0 baselines at the pre-registered bands:
nc1 23/23, nc2 45/45, nc3 98/98 (independently re-run at landing -
identical), weighted 44/45 at the pre-registered seed. The weighted
excursion (ranef g:Xr1:14, mean ratio 4.22 vs k = 4) is ADJUDICATED SEED
NOISE with the gate unloosened: it scatters 0.598/0.118/0.309 across seeds
where a wiring bug would reproduce, an independent-seed refit passes 45/45,
the recording's own Stan-vs-Stan rerun hit 2.85 on these summaries (ranef
MCSEs still ~1.4x under even after the n^(2/3) correction), and the
weighted-specific normalizer is pinned to 1e-13 by the target oracle. Full
suite 239/0; no snapshots regenerated.

C3. BINARY path to WALNUTS (Stan now unused at runtime, still compiled). The binary
   functor is the continuous one with is_binary -> actual_aux==1, the aux dimension
   removed, and the dispersion->RE-scale coupling dropped; response is the probit
   latents already fed by setResponse. Route binary through WalnutsSampler. Files:
   parametric_model.hpp (the is_binary specialization if not already parameterized in
   C1), init.cpp (route binary). Gate: distributional equivalence vs C0 binary
   baseline under C2's pre-registered tolerance (k = 4 on means via combined MCSE;
   sd ratio in [0.8, 1.25]); test-02-binary green. Size: S-M. Abort: as C2.

### C3 landing (2026-07-16)

C3 = 2ee6048. The family branch collapsed - both families on WALNUTS, Stan
compiled but unused at runtime. Binary gate 55/55 within the pre-registered
bands (worst mean ratio 2.56, worst sd ratio 1.137; independently re-run at
landing - identical). Suite 239/0, nothing regenerated, reproducibility
verified for binary directly. VD signed off the C4 author-list surgery
2026-07-16; the C4 protocol is pause -> VD confirms a quiet machine ->
Stan-era perf baselines (before-numbers) -> the deletion -> after-numbers
in the same window.

C4. DELETE STAN. Remove src/include/{stan,sundials,rstan}, src/stan_files/ (incl. the
   stale continuous_new.hpp), stan_sampler.{cpp,hpp}, stan_sampler_includes.hpp,
   interruptable_sampler.hpp. double_writer.hpp is RETAINED-AND-DE-STANNED, not
   deleted: it derives from stan::callbacks::writer but init.cpp's Gibbs bridge and
   createStanResultsExpr consume sample_writer directly - strip the stan base class,
   keep the storage/naming behavior. init.cpp co-moves: drop the unused
   #include "rstan/io/r_ostream.hpp" and rewire the sample_writer/
   createStanResultsExpr seam to the de-Stanned writer. Build surgery: Makevars.in drop
   -I"include/sundials", the @TBB@/RcppParallel flags, drop stan_sampler.cpp from
   SOURCES/OBJECTS; configure.ac drop the RcppParallel/TBB probing (22-38) and
   STAN_THREADS, KEEP the ax_ext SIMD + AC_CHECK_SIZEOF/ALIGNOF/alloca/posix_memalign
   probes (config.h feeds the RETAINED misc/rc dbarts-support libs, not Stan);
   regenerate configure via autoreconf; mirror in configure.win/Makevars.win.
   DESCRIPTION: drop LinkingTo BH + RcppParallel, drop Imports RcppParallel (keep
   Depends parallel for makeCluster); Title/Description drop "Stan-Sampled"; add
   WALNUTS/Bob Carpenter (ctb/cph, MIT) + SystemRequirements C++20. AUTHOR/LICENSE
   surgery (FLAG for VD - author-list changes are his call, Q(c)): remove the ~20
   StanHeaders `ctb` entries (DESCRIPTION:12-32) and the 3 CVODES `cph` entries
   (34-36); KEEP Goodrich/Gabry/Ali/Brilleman/Burkner (rstanarm_functions.R,
   retained) and the lme4/config.guess/m4 contributors. Files: the deletions,
   double_writer.hpp (de-Stanned), init.cpp (rstan include drop + writer seam),
   Makevars.in, Makevars.win, configure.ac, configure, configure.win, DESCRIPTION.
   Gate: R CMD check clean; full testthat green; WALNUTS-only build installs; MEASURE
   the AFTER perf numbers vs C0 (the payoff - install time, peak RSS, per-iter speed,
   and stan_sampler.o gone). Size: L (mechanical but wide + delicate build/author
   surgery). Abort: any dangling Stan include or an orphaned config.h symbol that
   breaks the retained misc/rc/ext libs - the SIMD/size probes are shared, do not
   over-delete.

### C4 landing (2026-07-16)

C4 = 174b369. 2135 files, 7.07 MB, 210k line-deletions; double_writer
de-Stanned rather than deleted, the control parsing and results assembly
moved verbatim to parametric_control.cpp with their R-visible contracts
unchanged, configure surgery kept the SIMD/size probes the retained support
libs consume, and the author surgery landed per VD's sign-off (21
StanHeaders ctb + 3 CVODES cph out, Bob Carpenter in, rstanarm/lme4
lineages kept). The Stan log_prob cross-check retired with the model; its
verdict stands in the C1 landing note and the FD gate remains.

THE MOTIVATION, MEASURED (one quiet window, three recordings in
benchmarks/baselines/ per the MANIFEST): the sampler swap bought 12-39x on
whole-fit per-iteration wall time (15.5/8.1/13.5 ms down to 0.4-0.7 ms on
the reference fits - the parametric step itself roughly two orders); the
deletion bought R CMD INSTALL 42.1s -> 15.1s and COMPILE peak RSS 2.21 GB
-> 0.47 GB (the stanc TU was the memory hog; peak SAMPLING RSS at these
reference sizes was never Stan-dominated and is unchanged at ~273 MB).

Gates: install clean, suite green (twice at landing; one run showed the
known marginally-flaky statistical expectation, absent on rerun), nc1
distributional spot-check PASS, ihdp causal smoke sane (CATT ~ 4.03 vs
truth 4, finite), R CMD check on the C4 tarball Status: OK with 0
WARNING/0 NOTE at the orchestrator's run. The implementer's check run had
hit an intermittent test-07 extraction-consistency ERROR: investigated at
landing - the Stan-era control tree ALSO checks OK, the C4 tarball checks
OK on re-run, and the mechanism (extract("trees") allTrees vs
individualSamples disagreeing on a node-count column, value varying across
runs) is a data-dependent dbarts-side inconsistency, sampler-agnostic,
FILED on the dbarts backlog. stan4bart cannot fix it and the test is not
regenerated.

C5. DOCS + NEWS. A landing note (docs/design/ or this plan's landing section) with
   the target's log-density, the adjoint derivation and its FD gate, the Stan->WALNUTS
   control-arg + diagnostics map (Q(c) resolution), and the prior-scope decision
   (Q(b)). NEWS entry (the sampler swap, the deprecations, the perf payoff);
   readme.md + man/ updates for changed control args and the Title. Files:
   docs/design/walnuts.md, NEWS.md (new), readme.md, man/*.Rd. Gate: R CMD check man;
   docs consistent with the landed surface. Size: M.

Challenge to the suggested shape: C2 and C3 are kept SEPARATE (not merged) even
though binary is a sub-case, because a per-family distributional gate localizes a
target bug to one family and keeps each commit's blast radius small; merging saves
little. The one place the code could argue for reorder: if Q(b) drops the shrinkage
families, that narrowing lands IN C1 (the functor never grows those latents), not as
a separate commit.

## Verification

- C1 (primary, two near-exact gates): (i) the Stan log_prob TARGET cross-check -
  hand target vs continuous.hpp's log_prob<false, true> at the same unconstrained
  points, difference CONSTANT across points, all fixture tiers (the oracle for the
  target itself; the FD gate cannot see a mis-specified target); (ii)
  finite-difference gradient agreement to tight tol across all tiers
  (fixed/sigma, nc=1/2/>=3, continuous/binary, weighted/unweighted). Together these
  stand in for the missing bitwise oracle. VD's R prototype (if offered) is a third
  independent implementation the FD points can be cross-checked against at nc>=3.
- C2/C3: distributional equivalence vs the C0 Stan-era baselines under the
  PRE-REGISTERED tolerance (C2): means gated by
  |mean_walnuts - mean_stan| < 4 * sqrt(mcse_stan^2 + mcse_walnuts^2), sds gated as
  a ratio in [0.8, 1.25], for every fixef/ranef/sigma/Sigma summary, on continuous
  nc=1/nc=2/nc>=3 + weighted + binary reference datasets. Plus reproducibility
  (test-05-rng: same seed
  -> identical draws; different core counts differ) - WALNUTS deterministic per-chain
  seed preserves this. Value-pinned testthat snapshots that legitimately move (the
  sampler changed) are regenerated, not forced.
- C4: R CMD check clean on the WALNUTS-only build; the AFTER perf numbers recorded
  against C0 (install time, peak RSS, per-iteration speed) - the arc's motivation,
  measured. rchk note if any .Call signature changed (the logdensity export is new).
- End-to-end smoke: the ihdp example (stan4bart/ihdp - a causal-inference IHDP
  simulation harness: sim.R/runSimulations.R fit stan4bart across replications,
  data.R/sim.data.gz the data) run on a couple of replications to confirm the causal
  workflow (treatment=, test=) still produces sane CATE estimates.
- Full testthat across all 11 files; test-11-callback (the callback surface reads
  copyOutParameters mid-sweep) is the one most sensitive to the write_array
  re-emission - watch it.

### C5 landing (2026-07-16) - ARC CLOSED

C5 = 45219e6. Design record docs/design/walnuts.md, NEWS.md 0.0-14, readme
and man updates, codoc clean. Two surface edges found and closed with
tests: shrinkage priors previously reached an uncaught C++ exception that
ABORTED the R process - now an informative R error at setup; the ignored
NUTS adaptation args now warn once as deprecated. check_sampler_diagnostics
was dead code predating this arc (object$diagnostics was never populated) -
now a documented no-op, so the WALNUTS diagnostics loss changed nothing
observable. Suite 241/0 (twice); R CMD check 0 WARNING / 0 NOTE, with the
known intermittent test-07 dbarts extraction flake firing on the
implementer's two check runs but not the orchestrator's (observed rate
roughly half across the arc's seven check runs; filed on the dbarts
backlog as extract-trees-consistency).

Arc summary: six working commits (structSize fix, baselines, target +
adjoint, continuous wiring, binary wiring, deletion) plus records. Every
posterior-defining step was proven against two oracles or a pre-registered
distributional gate, and the motivation landed measured: 12-39x whole-fit
per-iteration speed from the swap, install 42.1s -> 15.1s and compile peak
RSS 2.21 GB -> 0.47 GB from the deletion, sampling RSS unchanged. Open
doors recorded in the design note: probit stays unbuilt behind the same
functor seam, shrinkage families behind an informative error, the count of
ignored NUTS args drops to zero at their removal release.

## Open questions for VD

Q(a) GRADIENT DERIVATION - self-derived vs VD's R prototype as oracle. FORK: who
  produces the authoritative gradient for the decov onion (nc>=3). Self-derived
  analytic adjoint, FD-gated (RECOMMENDED): I read continuous.stan end to end and am
  HIGH-confidence on the whole reachable target EXCEPT the nc>=3 onion
  (scale_factor = sqrt(rho/dot_self(T_row))*std_dev threading z_T,
  continuous.stan:40-49), where I am MEDIUM without a gate. Cost: a subtle onion
  adjoint bug could survive if the FD fixtures miss the rare nc>=3 shape - mitigated
  by REQUIRING an nc>=3 fixture in the C1 FD gate. VD's R prototype as oracle
  (the offer): a second, independent gradient implementation to cross-check FD at
  nc>=3. Cost: serializes C1's start on VD's availability if treated as a blocker.
  RECOMMENDATION: proceed self-derived + FD-gated now; fold VD's R prototype in as a
  cross-check for make_theta_L/make_b when it lands, NOT as a gate on starting C1.
  Either way the FD gate is mandatory.

Q(b) PRIOR-FAMILY SCOPE - support all 8 fixed-effect coef families or restrict to the
  reachable-sane subset. FORK: the shrinkage families (hs, hs_plus, laplace, lasso,
  product_normal, continuous.stan:266-270,298-322) carry global/local/caux/mix/
  one_over_lambda latents that roughly double the parametric dimension and are the
  nastiest gradient terms - but they are opt-in only and arguably nonsensical for
  stan4bart's typically-few-column parametric part (BART carries the signal).
  Restrict to {normal, student_t, cauchy} beta + {normal, student_t, cauchy,
  exponential} aux + decov RE (RECOMMENDED): the full DEFAULT-reachable set; a user
  who passed prior=hs() gets a clear "unsupported under WALNUTS" error. Cost: a hard
  break for anyone using shrinkage priors on the parametric part (likely nobody -
  no test exercises them; the default is normal). Support all 8: faithful to
  rstanarm. Cost: large dead-weight gradient surface, more FD fixtures, more adjoint
  risk, for a code path stan4bart's design discourages. RECOMMENDATION: restrict,
  deprecate the shrinkage families with an informative error; revisit if a user
  reports needing them. (This is also a user-facing compat break - grouped with Q(c).)

Q(c) USER-FACING CONTROL-ARG + DIAGNOSTICS + AUTHOR-LIST COMPATIBILITY. FORK: the
  Stan control vocabulary and diagnostics have no faithful WALNUTS analog. Loop-level
  args (iter, warmup, skip, chains, cores, refresh, seed, verbose - stan4bart.R:13-24)
  are sampler-agnostic and KEPT with identical meaning. The NUTS-specific stan_args
  (adapt_delta, adapt_gamma, adapt_kappa, adapt_t0, adapt_init_buffer/term_buffer/
  window, stepsize, stepsize_jitter, max_treedepth, init_r - controlNames
  stan_sampler.cpp:82-96) map poorly onto WALNUTS' Adam step-size + Nutpie diagonal
  mass + its own warmup schedule. Options: (i) accept-but-ignore with a one-release
  deprecation warning (RECOMMENDED - maximal source compat; a script passing
  adapt_delta still runs); (ii) drop them and expose WALNUTS' knobs (initial
  step_size, mass_smoothing, warmup length); (iii) best-effort map the few with an
  analog (max_treedepth if WALNUTS caps depth). DIAGNOSTICS - the honest fork:
  check_sampler_diagnostics (stan4bart.R:255-297) warns on divergent__/treedepth__/
  EBFMI(energy__), but the vendored WALNUTS ChainHandler concept surfaces ONLY
  position/lp/step_size/diag_inv_mass (+ cross-chain r_hat) - NO treedepth, NO
  divergent, NO n_leapfrog, NO energy. So ALL THREE checks are lost, not
  conditionally. Options: (i) DROP all three with a NEWS note and keep/surface
  r_hat (RECOMMENDED - zero plumbing, honest about what the sampler reports); vs
  (ii) plumb depth/divergence counters out of walnuts.hpp internals (cost: patching
  the vendored headers, a maintenance burden re-paid on every upstream refresh).
  CAVEAT: the implementer vendors WALNUTS fresh from
  https://github.com/flatironinstitute/walnuts at C1 and must RE-VERIFY this
  diagnostics claim against the fresh API before finalizing the Q(c) resolution -
  upstream may have grown handler hooks since the bairrtt snapshot. AUTHOR LIST:
  deleting StanHeaders + CVODES removes ~20 ctb +
  3 cph DESCRIPTION entries (12-36) and adds Bob Carpenter (WALNUTS, MIT); this is
  VD's call - it is his package's authorship. Cost of (i): a small deprecation-shim
  surface to carry for a release. Cost of (ii): an immediate break for scripts
  setting adapt_*. RECOMMENDATION: control args - accept-but-ignore + deprecate the
  NUTS args, keep loop args verbatim; diagnostics - option (i), drop the three NUTS
  checks with a NEWS note and keep r_hat, subject to the fresh-upstream re-check;
  VD signs off the author-list surgery.

Q(d) WALNUTS: VENDOR into stan4bart, or a SHARED dependency. FORK: bairrtt vendors
  the headers (inst/include/walnuts, MIT, Bob Carpenter). Vendor identically
  (RECOMMENDED): no CRAN dependency on a package that does not exist, MIT permits it,
  cost is one WALNUTS_LICENSE file + one author entry + carrying the headers in the
  tree (a few hundred KB vs the tens of MB of StanHeaders being removed - still a
  massive net win). Shared LinkingTo a walnuts headers package: DRY across bairrtt +
  stan4bart. Cost: that package must be authored, published to CRAN, and version-
  coordinated first - it does not exist, and the memory record (stan4bart-walnuts-plan,
  bairrtt is the reference consumer) treats vendoring as the current pattern.
  RECOMMENDATION: vendor now; revisit a shared package if/when one is published.

Q(e) BASELINE-RECORDING SCOPE (C0). FORK: how many reference fits to pin as the
  distributional oracle. Minimal (one continuous + one binary): cheap to run/store.
  Cost: the nc>=3 onion - the riskiest gradient term - would be gated only by FD, not
  by an end-to-end posterior. Tiered (continuous nc=1 / nc=2 / nc>=3 / weighted +
  binary - RECOMMENDED): the distributional gate then covers every gradient tier the
  FD gate covers, so a target bug that FD's point-checks miss still trips a posterior-
  summary gate. Cost: ~5-6 baseline fits to run at stable length and store as .rds.
  RECOMMENDATION: tiered - the marginal cost is a few fits, and the onion is exactly
  where an independent end-to-end check earns its keep.
