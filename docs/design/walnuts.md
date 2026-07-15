# walnuts (design record)

This is the retrospective design record for the WALNUTS arc (docs/plans/walnuts.md
holds the drafting/execution record - commit partition, open questions, landing
notes with the measured numbers this file cites). Read that file for the "why now"
and the commit-by-commit history; this file is the "what shipped and why it is
correct" record, written at C5 once the arc closed.

## The parametric target (src/parametric_model.hpp)

The WALNUTS target is `ParametricModel`, a forward port of continuous.stan's
`parameters`/`model` blocks plus a hand-written reverse-mode adjoint. It reproduces
continuous.stan's log-density exactly (up to a parameter-independent constant),
including its quirks. The free (unconstrained) position vector, in continuous.stan's
`parameters` declaration order:

```
[ z_beta(K), z_b(q), z_T(len_z_T), rho_free(len_rho),
  zeta_free(len_concentration), tau_free(t), aux_unscaled_free(!is_binary) ]
```

Transforms (parametric_model.hpp:246-266): `rho = inv_logit(rho_free)`,
`zeta = exp(zeta_free)`, `tau = exp(tau_free)` (all per-block/per-element), and for
continuous responses `aux_unscaled = exp(aux_unscaled_free)` then
`aux = prior_scale_for_aux * aux_unscaled [+ prior_mean_for_aux]` depending on the
aux prior family; `beta` is either the affine Normal transform
(`z_beta * prior_scale + prior_mean`) or, for student_t/cauchy, the Cornish-Fisher
polynomial `CFt(z_beta, df)` (parametric_model.hpp:36-55, ported from
continuous.stan:146-158) composed with the same affine map. `make_theta_L`/`make_b`
(parametric_model.hpp:269-349) are a direct port of continuous.stan's decov "onion"
construction: each grouping-factor block gets a lower-triangular Cholesky factor `T`
built from `tau` (overall scale), a Dirichlet simplex over `zeta` (splitting the
trace across the block's `nc` correlated coefficients), `rho` (canonical partial
correlations, `T21 = 2*rho-1` for `nc==2`), and, for `nc>=3`, the "onion method"
threading `z_T` through `scale_factor = sqrt(rho/dot_self(T_row)) * sd`; `b_level =
T_i * z_b_level`.

The log-density (parametric_model.hpp:382-556; parameter-independent constants
dropped so the value matches Stan's `log_prob` up to a constant, not exactly):
Gaussian likelihood `-N*log(actual_aux) - 0.5*lambda*sum_i w_i*(y_i-eta_i)^2` with
`lambda = 1/actual_aux^2`, `eta = offset_ + X*beta + Z*b` (Z sparse, CSR); Normal(0,1)
priors on `z_beta`/`z_b`/`z_T`; Beta(a1,a2) onion-regularization priors on `rho`;
Gamma priors on `zeta`/`tau`; the matching logit/log constraint Jacobians added by
hand exactly where Stan's `transformed parameters` block would add them
automatically. **The weighted-normalizer quirk**: the likelihood's normalizing term
is `-N*log(actual_aux)`, using `N` (the observation count), NOT `sum(w)` (the
summed weight) even in the weighted case - continuous.stan:365's deliberate choice,
reproduced verbatim at parametric_model.hpp:392. Getting this wrong shifts the
posterior on sigma under weighting; the C1 Stan `log_prob` cross-check (below)
pinned this term to the two implementations agreeing to ~1e-13 on the weighted
tier specifically. Binary responses reduce to the same Gaussian model with
`actual_aux == 1` (the `aux` dimension absent entirely, not merely fixed) over the
BART-side probit latents that `setResponse` refreshes each Gibbs sweep; the
`aux -> dispersion -> RE-scale` coupling (random-effect blocks are scaled by
`tau*scale*dispersion`, and `dispersion == aux` for continuous responses) drops
along with it.

## The hand adjoint and its two-oracle validation history

The forward pass above is plain arithmetic - no autodiff needed. The work is the
reverse pass: a hand-written backprop mirroring the forward loops term by term
(parametric_model.hpp:394-565), NOT a general autodiff (reintroducing that
dependency is exactly what this arc removes). The riskiest term is the `nc>=3`
onion, where `z_T` couples through a `dot_self` and a square root
(parametric_model.hpp:450-471); every other term (the Gaussian likelihood, the
affine/Cornish-Fisher beta transforms, the `nc==1`/`nc==2` closed-form theta_L, the
rho/zeta/tau priors and Jacobians) is mechanically simple.

Two independent oracles validated this code, built BEFORE it was trusted (C1):

1. **Finite-difference gradient gate** (tests/testthat/test-12-gradient.R, still
   live): central differences of the hand-written value vs the hand-written
   gradient at fixed-seed random unconstrained points, across ten fixtures
   spanning every tier - fixed effects only, nc=1/nc=2/nc>=3 (onion),
   weighted, a mixed multi-block fit (added at review: the aux -> dispersion
   adjoint must ACCUMULATE across heterogeneous blocks, not just apply once),
   binary, and student_t betas/aux priors. Landed result: worst relative error
   2.4e-8 across all ten tiers (tol 1e-6).
2. **Stan `log_prob` target cross-check** (retired at C4, verdict frozen in the
   C1 landing note of docs/plans/walnuts.md): while Stan's `continuous.hpp` was
   still compiled, the hand-written value was evaluated against Stan's own
   `log_prob<false, true>` at the SAME unconstrained points across all fixture
   tiers, requiring the difference to be CONSTANT across points (targets may
   differ by dropped constants, but not by a function of the parameters - that
   would mean a mis-specified term). This is the oracle the FD gate cannot
   provide: FD only proves gradient-vs-coded-target agreement, not that the
   coded target itself is continuous.stan's target. Landed result: worst
   spread 9.1e-13 across all tiers (tol 1e-6), including the weighted
   normalizer and the binary reduction confirmed independently. This gate was
   REMOVED at C4 along with the rest of the Stan machinery it depended on
   (there is no `continuous.hpp` left to call) - it could retire precisely
   because its job was to validate the TARGET once, not to run forever; the
   FD gate remains as the standing regression gate on gradient-vs-target
   agreement; a future target bug (e.g. from a refactor) would need a new
   value-level oracle since this one is gone for good.

## The WALNUTS integration

WALNUTS is vendored (not a CRAN/shared dependency - Q(d) of docs/plans/walnuts.md)
into inst/include/{walnuts/*, walnuts.hpp, WALNUTS_LICENSE, WALNUTS_VERSION}, pulled
FRESH from https://github.com/flatironinstitute/walnuts at commit
5854be888e5432a103275466d7d0be95c7c5c67a (upstream date 2026-07-13, vendored
2026-07-15; MIT, (c) 2025 Bob Carpenter) - bairrtt's copy was used only as an
integration reference, re-verified against the fresh headers at C1 (see
inst/include/WALNUTS_VERSION for the drift notes; the only drift, a `WarmupConfig`
step-size convergence tolerance and an internal `detail::sample()` signature
change, does not touch the per-chain path this package uses).

`src/walnuts_sampler.{hpp,cpp}` presents exactly the surface init.cpp's Gibbs
bridge consumes (the `ParametricSampler` base: `run(isWarmup)`,
`getParametricMean`, `getSigma`, `setOffset`/`setResponse`, `freeze()`), mirroring
bairrtt's `IrtSampler` lifecycle: construct in the adapting phase (one
`AdaptiveWalnuts` per chain, matching the current one-process-per-chain
parallelism - WALNUTS' cross-chain mass averaging is unused by design), drive one
ADAPTING transition per warmup sweep (`run(true)` -> `(*adapter)()`), `freeze()` at
`disengageAdaptation` (`AdaptiveWalnuts::sampler()` hands off to a
`WalnutsSampler`, `adapter.reset()`), then one FIXED-TUNING transition per sampling
sweep (`run(false)` -> `(*sampler)()`). `WarmupConfig` is sized to the outer loop's
actual warmup sweep count so WALNUTS' internal mass-adaptation windows schedule
correctly against the one-transition-per-sweep cadence the Gibbs loop drives
(walnuts_sampler.cpp:224-247). `setOffset`/`setResponse` are in-place `memcpy`s
into the functor's held `offset_`/`y_` vectors (the bairrtt
`logp_grad.theta = ...` pattern) - no functor reconstruction between sweeps.

**The non-finite throw-contract lesson** (parametric_model.hpp:558-564): WALNUTS
assumes Stan's `log_prob` contract that a non-finite result is signaled by a
THROW, not by silently returning `inf`/`nan`. At extreme leapfrog positions,
`exp()` overflow in `tau`/`zeta`/`aux` can produce a non-finite log density or
gradient; returning it silently, as a naive port might, poisons the Adam
step-size adaptation (which integrates the bad value) into a crash much later,
at `freeze()`, far from the actual cause. `ParametricModel::eval` now explicitly
checks `isfinite`/`allFinite` and throws `std::domain_error`, which WALNUTS'
`NoExceptLogpGrad` wrapper catches and turns into `-inf` logp / zero gradient,
correctly treated as a divergence at that step. This was a real integration bug
caught during C2, not a hypothetical.

## The distributional-equivalence gate design

There is no bitwise oracle for the DRAWS (swapping NUTS for WALNUTS changes the
sampler wholesale). Correctness instead rests on the two target/gradient oracles
above plus a THIRD gate: distributional equivalence between the WALNUTS fit and
Stan-era baselines on the same reference data, PRE-REGISTERED before any WALNUTS
draws existed (the tolerance was fixed at C0/C1, before C2 could see the numbers
it would be judged against). For every gated summary (fixef/ranef/sigma/Sigma,
mean AND sd):

```
means: |mean_walnuts - mean_stan| < 4 * sqrt(mcse_stan^2 + mcse_walnuts^2)
sds:   ratio in [0.8, 1.25]
```

`k = 4` and the sd band were fixed in advance (docs/plans/walnuts.md C2), not
tuned post hoc. Five reference tiers cover every gradient-adjoint tier the FD gate
covers, so a target bug FD's point-checks miss would still trip a posterior-summary
gate: continuous nc=1, nc=2, nc>=3 (onion), weighted continuous, and binary
multilevel. MCSEs (both Stan-era and WALNUTS-era) were computed as
`max(n^(2/3)-batch-means, between-chain)` since plain `sqrt(n)` batch-means MCSEs
read 3-8x too small under this workload's slow BART-vs-parametric alternation
(docs/plans/walnuts.md C0 landing note).

Landed results (independently re-run at each landing - identical both times):
nc1 23/23 within band, nc2 45/45, nc3 98/98, binary 55/55 (worst mean ratio 2.56,
worst sd ratio 1.137). The weighted tier landed 44/45: one ranef summary
(`g:Xr1:14`, mean ratio 4.22 vs `k=4`) sat outside the band.

**Weighted-tier seed-noise adjudication.** This excursion was investigated, not
waived: the gate itself was not loosened. Evidence it is seed noise, not a wiring
bug: (1) the same summary's ratio scatters 0.598/0.118/0.309 across independent
seeds, where a real bug would reproduce; (2) an independent-seed refit passes
45/45 outright; (3) a Stan-vs-Stan rerun (same sampler both sides, pure seed
variation) hit 2.85 on these exact summaries - the baseline's own ranef MCSEs run
about 1.4x under even after the `n^(2/3)` correction, so a 4.2 excursion on one
summary out of 45 is within the tail this baseline noise implies; (4) the
weighted-specific term (the `N`-not-`sum(w)` normalizer) is independently pinned
to 1e-13 agreement by the C1 Stan `log_prob` cross-check, ruling out a target bug
specific to weighting. Full testthat suite was 239/0 at that landing; nothing was
regenerated.

## Diagnostics: dropped, not conditionally

Stan's `check_sampler_diagnostics` (R/stan4bart.R) used to warn on divergent
transitions, max-treedepth transitions, and low E-BFMI (from `energy__`). The
vendored WALNUTS `ChainHandler` concept (inst/include/walnuts/concepts.hpp)
surfaces only `position`/`lp`/`step_size`/`diag_inv_mass` (`GlobalHandler` adds a
cross-chain `on_r_hat`, but this package never wires a `GlobalHandler` in: chains
run as independent R-level processes/threads, not a single joint WALNUTS run, so
there is no cross-chain handle to attach it to). There is NO treedepth, divergent,
n_leapfrog, or energy hook at all - re-verified against the freshly-vendored
headers at C1, not assumed from the bairrtt-era snapshot. So all three checks are
lost categorically, not case-by-case. `WalnutsSampler::run` still writes
`stepsize__`/`treedepth__`/`n_leapfrog__`/`divergent__`/`energy__` rows
(walnuts_sampler.cpp:290-293), but as constant zero - a Stan-row-layout
placeholder so `result$stan`'s dimnames are unchanged, not real diagnostics.
`check_sampler_diagnostics` (R/stan4bart.R) is now a documented no-op: inspecting
the actual pre-WALNUTS code turned up that `object$diagnostics` - the field every
one of its checks read - was never populated in the first place (the real data
lives in `result$stan`, keyed by a `"parameters"` dimname, not a `"diagnostic"`
one), so every check had already been unreachable dead code before this arc
touched it. The function is kept (rather than deleted) as a single documented
call site, in case WALNUTS or a future handler grows a real diagnostic to plug in.

## The prior-scope restriction (Q(b))

continuous.stan's shrinkage coefficient priors (`hs`, `hs_plus`, `laplace`,
`lasso`, `product_normal`) carry `global`/`local`/`caux`/`mix`/`one_over_lambda`
latents that roughly double the parametric dimension and are the nastiest
gradient terms in the whole model - and they are opt-in only (the default is
`normal`; no test exercised them) and arguably nonsensical for stan4bart's
typically-few-column parametric part, since BART carries the signal. `
ParametricModel` implements only `{none, normal, student_t/cauchy}` for beta and
`{none, normal, student_t, exponential}` for aux (parametric_model.hpp:12-15);
`ParametricModel::finalize()` throws `std::invalid_argument` if a `prior_dist`
code outside that set reaches it. Before C5, nothing stopped a `prior_dist` code
for the shrinkage families from reaching that throw: the R surface's
`handle_glm_prior` (R/rstanarm_functions.R) still accepts `hs`/`hs_plus`/`lasso`/
`laplace`/`product_normal` (`ok_dists`, stan4bart_fit.R:132-133, unchanged - this
list is shared rstanarm-derived code, not narrowed), so a user passing
`stan_args = list(prior = hs())` reached `ParametricModel::finalize()`'s throw
uncaught across the `.Call` boundary - `libc++abi: terminating due to uncaught
exception ...`, a hard process abort, not a catchable R error (verified directly
before landing this fix). C5 adds `refuse_shrinkage_prior()` (R/stan4bart_fit.R),
called before `handle_glm_prior`, raising an informative `stop()` naming the
unsupported families and the supported alternative, entirely at the R level,
before any C++ runs. `prior_intercept` and `prior_aux` were already scoped to
non-shrinkage families (`ok_intercept_dists`/`ok_aux_dists`) and needed no change.

## The control-arg compatibility resolution (Q(c))

The NUTS control vocabulary (`adapt_delta`, `adapt_gamma`, `adapt_kappa`,
`adapt_t0`, `adapt_init_buffer`, `adapt_term_buffer`, `adapt_window`, `stepsize`,
`stepsize_jitter`, `max_treedepth`) has no faithful WALNUTS analog (WALNUTS uses
Adam step-size adaptation and a Nutpie-style diagonal mass estimator with its own
warmup schedule). These are parsed into `StanControl`
(src/parametric_control.cpp - the `controlNames` array, ~line 31) purely for
accept-but-ignore source compatibility; `WalnutsSampler` never reads them. Before
C5 this was silent - a script setting `adapt_delta` ran and nothing told the user
their setting had no effect. C5 adds `warn_ignored_nuts_args()`
(R/stan4bart_fit.R), called as soon as `stan_args` is available (before any of the
prior/QR/control processing that follows), issuing a single warning per fit
naming every deprecated arg the caller supplied: `"'<arg>' is ignored by the
gradient-based sampler and is deprecated; it will be removed in a future
release"`. `init_r` is deliberately excluded - it still has an effect: it sizes
the uniform half-width WALNUTS draws its initial unconstrained position from
(`WalnutsSampler`'s ctor, walnuts_sampler.cpp:236-240), mirroring Stan's `init_r`
random-start radius. `seed`/`skip` are loop-level (threaded to the WALNUTS rng
seed and the outer sweep-skip cadence respectively) and also excluded, alongside
the top-level `iter`/`warmup`/`chains`/`cores`/`refresh`/`verbose` arguments,
which were never part of `stan_args` in the first place.

## The measured payoff (benchmarks/baselines/MANIFEST)

Three recordings, one quiet window, arm64 macOS, reference fits (continuous nc=2 /
binary multilevel / weighted continuous):

| stage                          | install | compile peak RSS | per-iter (cont/bin/wtd)      | sampling RSS   |
|---------------------------------|---------|-------------------|-------------------------------|----------------|
| before (cb809b1, Stan, all runtime) | 42.1s | 2.21 GB | 15.5 / 8.1 / 13.5 ms | 273/274/271 MB |
| pre-deletion (2ee6048, WALNUTS runtime, Stan still compiled) | 43.5s | (unchanged) | 0.4 / 0.7 / 0.6 ms | 273-283 MB |
| after (174b369, WALNUTS only, Stan deleted) | 15.1s | 0.47 GB | 0.5 / 0.5 / 0.6 ms | 271-274 MB |

The sampler swap (before -> pre-deletion) bought 12-39x on whole-fit
per-iteration wall time - the parametric step itself roughly two orders of
magnitude, dwarfed at the whole-fit level by BART's own per-iteration cost.
Sampling RSS was NEVER Stan-dominated at these reference sizes and stays flat
throughout - the win is elsewhere. The deletion (pre-deletion -> after) bought
`R CMD INSTALL` 42.1s -> 15.1s (2.8x) and COMPILE peak RSS 2.21 GB -> 0.47 GB
(4.7x) - the `stanc`-generated `continuous.hpp` translation unit was the memory
and compile-time hog, not the runtime. Both numbers were MEASURED, per the
plan's charter, not asserted.

## Recorded follow-ups

- **The probit-path door stays closed.** stan4bart requires `family$link ==
  "probit"` for binary responses (R/stan4bart_fit.R) and hard-codes it in
  `stan4bart.R`'s family detection; a generic-link check is present in
  commented-out form (`#if (!length(link)) ...`) but was never enabled, before
  or after this arc. The binary reduction (`actual_aux == 1`, probit latents
  drawn on the BART side) is specific to probit; nothing in this arc opens the
  door to other links, and nothing in this arc closes it further than it
  already was.
- **Shrinkage coefficient priors, if ever demanded.** Q(b) restricts to
  `{normal, student_t, cauchy}` because no test and no known user exercises the
  shrinkage families, and BART already carries the signal they would compete
  for. If a user reports needing `hs()`/`lasso()`/etc. on the parametric part,
  the path back in is `ParametricModel`'s prior-scope check
  (parametric_model.hpp:110-117) plus a matching forward/adjoint extension for
  each family's extra latents - not a redesign.
- **The dbarts `extract("trees")` flake.** The C4 landing hit an intermittent
  `test-07-extractedTrees.R` extraction-consistency `ERROR` under `R CMD check`:
  `allTrees` vs `individualSamples` disagreeing on a node-count column, with the
  disagreement's presence varying across runs. Investigated at landing: both the
  Stan-era control tree and the post-deletion WALNUTS tarball check OK on rerun,
  and the mechanism is data-dependent and sampler-agnostic (reproduces
  independent of whether Stan or WALNUTS drew the parametric part) - it is a
  dbarts-side inconsistency in tree extraction, filed on the dbarts backlog, not
  something stan4bart can fix locally. The test is not regenerated or weakened;
  a rerun of `R CMD check` is the documented workaround if it fires.
