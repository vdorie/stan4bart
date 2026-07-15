# stan4bart 0.0-14

## Sampler

* The parametric conditional (fixed effects, lme4-style random effects, and
  the residual sd for continuous responses) is now drawn with a hand-derived
  analytic log-posterior gradient and the WALNUTS sampler (vendored from
  <https://github.com/flatironinstitute/walnuts>, commit 5854be8; MIT
  licensed, (c) 2025 Bob Carpenter), replacing the embedded Stan/NUTS
  sampler for both continuous and binary response families. This is a
  sampler swap only: the outer BART-vs-parametric Gibbs alternation, the R
  formula/data-prep surface (`glFormula`, lme4 grouping syntax, priors), and
  the posterior targeted are unchanged, verified by a pre-registered
  distributional-equivalence gate against the Stan-era posterior on
  reference fits spanning every gradient tier (nc=1/nc=2/nc>=3 random-effect
  structures, weighted, and binary). See `docs/design/walnuts.md` for the
  full design record, including the gate design and results.

* Measured payoff (one quiet window, arm64 macOS, reference fits spanning
  continuous/binary/weighted-continuous): whole-fit per-iteration wall time
  dropped 12-39x (15.5/8.1/13.5 ms down to 0.4-0.7 ms). Deleting the Stan/
  StanHeaders/sundials/rstan machinery dropped `R CMD INSTALL` time from
  42.1s to 15.1s (2.8x) and compile-time peak RSS from 2.21 GB to 0.47 GB
  (4.7x); peak *sampling* RSS was never Stan-dominated at these reference
  sizes and is unchanged (~273 MB throughout).

* Dependencies: StanHeaders, BH, and RcppParallel are no longer linked to or
  imported. C++20 is now required (`SystemRequirements: C++20`). WALNUTS
  (MIT license, Bob Carpenter) is vendored under `inst/include/walnuts`.

## Storage

* Warmup draws are no longer stored by default (`save_warmup = FALSE`), and the
  returned parametric ("stan") store keeps only the transformed rows every
  consumer reads (`beta`, `b`, `theta_L`, `aux`) plus the two live diagnostics
  (`lp__`, `stepsize__`). The raw unconstrained rows (`z_beta`, `z_b`, `z_T`,
  `rho`, `zeta`, `tau`, `aux_unscaled`) and the five constant-zero placeholder
  diagnostic rows (`accept_stat__`, `treedepth__`, `n_leapfrog__`,
  `divergent__`, `energy__`) - none of which is read by any computed surface -
  are dropped from default storage. Not storing warmup roughly halves the fit
  object at scale (`bart_train` is ~50% of the object and its warmup copy is the
  same again); every quantitative convergence diagnostic in the 2026 toolchain
  is defined on post-warmup draws only, matching cmdstanr/rstanarm/brms
  defaults. In place of full warmup the fit gains an `adaptation` component: per
  chain the frozen step size and diagonal inverse mass (which WALNUTS tuned and
  the previous build discarded), a warmup-end position snapshot, and a thinned
  warmup trace of the monitored scalars. Recomputing `bart_train` from stored
  trees was considered and left out of scope; `bart_train` and the derivable
  `f` functionals are unchanged. See `docs/plans/sample-storage.md` for the
  measurements and rationale.

* Two opt-ins restore the old behavior: `save_warmup = TRUE` stores the full
  per-draw warmup under `warmup` (and re-enables the `include_warmup = TRUE`
  accessors; on a default fit they now error informatively), and
  `stan_args = list(save_raw_parameters = TRUE)` restores the raw unconstrained
  rows (funnel forensics). The stored draws themselves are unchanged - the
  sampler path is untouched - so posterior summaries are identical to prior
  releases regardless of the storage flags.

## Deprecated

* The NUTS-specific `stan_args` controls - `adapt_delta`, `adapt_gamma`,
  `adapt_kappa`, `adapt_t0`, `adapt_init_buffer`, `adapt_term_buffer`,
  `adapt_window`, `stepsize`, `stepsize_jitter`, and `max_treedepth` - have
  no analog under WALNUTS. They are still accepted (a script that sets one
  does not break), but are now ignored, and a warning naming every supplied
  deprecated argument is issued at fit time. They will be removed in a
  future release. `init_r` (the initial-position radius) and the loop-level
  arguments (`iter`, `warmup`, `skip`, `chains`, `cores`, `refresh`, `seed`,
  `verbose`) keep their exact prior meaning and do not warn.

## Removed

* Sampler diagnostics tied to NUTS - divergent transitions, max-treedepth
  transitions, and low E-BFMI, previously warned on by
  `check_sampler_diagnostics` - are gone. The vendored WALNUTS sampler
  reports no analog of any of the three (it exposes only position, log
  density, step size, and the diagonal mass estimate), so the corresponding
  warnings no longer fire.

## Breaking changes

* The shrinkage coefficient prior families (`hs`, `hs_plus`, `lasso`,
  `laplace`, `product_normal`) are no longer supported for
  `stan_args$prior`. Supplying one now raises an informative error at fit
  setup ("prior families hs, hs_plus, lasso, laplace, and product_normal are
  not supported by the gradient-based sampler; use normal, student_t, or
  cauchy"); previously this reached the sampler and crashed the R session
  with an uncaught C++ exception. Use `normal`, `student_t`, or `cauchy`
  instead; the residual-sd (`prior_aux`) and covariance (`prior_covariance`)
  priors are unaffected - they were already restricted to non-shrinkage
  families.

## Bug fixes

* Fixed `dbarts_results.structSize` never being set by the dbarts 1.0
  flat-C-API port, which caused the versioned-struct field gate to skip
  populating every run's output buffers (usually masked by the buffers
  otherwise reading as zero). Found and fixed while recording this arc's
  Stan-era baselines.
