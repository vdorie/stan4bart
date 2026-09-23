# stan4bart 0.0-14

## Sampler

* The parametric conditional (fixed effects, lme4-style random effects, and
  the residual sd for continuous responses) is now drawn with a hand-derived
  analytic log-posterior gradient and the WALNUTS sampler (vendored from
  <https://github.com/flatironinstitute/walnuts>, commit f3c1833; MIT
  licensed, (c) 2025--2026 by the Walnutpie Developers), replacing the
  embedded Stan/NUTS sampler for both continuous and binary response
  families. This is a sampler swap only: the outer BART-vs-parametric Gibbs
  alternation, the R formula/data-prep surface (`glFormula`, lme4 grouping
  syntax, priors), and the posterior targeted are unchanged, verified by a
  pre-registered distributional-equivalence gate against the Stan-era
  posterior on reference fits spanning every gradient tier (nc=1/nc=2/nc>=3
  random-effect structures, weighted, and binary). See
  `docs/design/walnuts.md` for the full design record, including the gate
  design and results.

* Measured payoff (one quiet window, arm64 macOS, reference fits spanning
  continuous/binary/weighted-continuous, with dbarts held at one version so
  the comparison isolates the sampler): whole-fit per-iteration wall time
  dropped 53-95x (15.4/7.5/11.9 ms down to 0.16/0.14/0.16 ms). Deleting
  the Stan/StanHeaders/sundials/rstan machinery dropped `R CMD INSTALL`
  time from 42.1s to 15.1s (2.8x) and compile-time peak RSS from 2.21 GB
  to 0.47 GB (4.7x); peak *sampling* RSS was never Stan-dominated at these
  reference sizes and is unchanged (~273 MB throughout).

* Dependencies: StanHeaders, BH, and RcppParallel are no longer linked to or
  imported. C++20 is now required (`SystemRequirements: C++20`). WALNUTS
  (MIT license, Bob Carpenter) is vendored under `inst/include/walnuts`.

* The warmup-phase initial mass is now seeded from the log-posterior gradient
  at each chain's initial unconstrained position (`(1 - s) * |grad| + s`, a
  Nutpie-style heuristic; `s` is the existing mass-smoothing constant),
  replacing the previous identity (unit-mass) start; the initial step size is
  unchanged. Sampled values for a given seed change as a result - this is a
  draw-moving change, verified to stay within the pre-registered
  distributional-equivalence band against the Stan-era baseline on all
  reference tiers. It improves tuning at short warmup, especially for
  large-`n` fits with many random-effect levels.

* `fit$adaptation` gains `mean_leapfrog` and `mean_leapfrog_warmup`: the mean
  number of leapfrog (gradient evaluation) steps per transition, per chain,
  during the sampling and warmup phases respectively. This is the
  runtime-relevant tuning diagnostic - a high `mean_leapfrog` relative to a
  well-tuned fit signals under-warming - and is exact and draw-neutral (an
  eval counter, not a new sampling path).

* New `print` and `summary` methods for `stan4bartFit` fits. Both warn when
  the sampling-phase `mean_leapfrog` looks poor ("parametric sampler tuning
  looks poor ...; consider increasing warmup"), and `?stan4bart` documents a
  warmup floor: if `iter` is shortened below the default, keep `warmup` at or
  above roughly 100 for large-`n` fits with many random-effect levels.

* `stan_args$adapt_delta` is a live control again, mapped to WALNUTS'
  step-size acceptance-rate target (the analog of Stan's `adapt_delta`),
  validated to `(0, 1)`. Its default (0.8) reproduces prior behavior exactly
  - an unset or explicit-0.8 fit is bit-identical to earlier builds - so this
  is draw-neutral unless set. A higher value targets a smaller step size
  (more gradient evaluations, finer geometry tracking); a lower value takes
  larger, cheaper steps. Setting it no longer warns.

* `fitted(type = "ppd")` notes, once per session, that the posterior-predictive
  mean equals the expected-value mean, so `type = "ev"` computes it exactly
  and faster. The computation itself is unchanged: under a fixed seed,
  `fitted(type = "ppd")` still reproduces an average of
  `extract(type = "ppd")` bitwise.

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

* New `store` argument (`c("fits", "trees")`, default `"fits"`). `store =
  "trees"` keeps only the sampled trees and recomputes the `bart_train` /
  `bart_test` blocks on demand through the `dbarts` predict path, instead of
  retaining the `n x draws x chains` blocks in the fit. It implies `keepTrees`
  (and errors on a contradictory `bart_args = list(keepTrees = FALSE)`), keeps
  the parametric `stan` and `bart_varcount` blocks, and leaves `bart_train` /
  `bart_test` absent - read the BART fits through `extract`, `fitted`, or
  `predict`, all of which route through the recompute seam. This is an opt-in
  memory/time trade: the `bart_train` block is 94-98% of a large fit object
  (`n >= 10000`), so dropping it is a large saving above a crossover near
  `n = 550`; in exchange `extract` materializes the whole block on each call
  (~15 s at `n = 10000`, 4 chains x 1000 draws) while `fitted` streams the
  posterior mean in row blocks to stay memory-bounded. The recomputed values
  match the stored path to a tight tolerance (~1e-13; the stored block carries
  an extra offset round-trip and is the noisier quantity), and sampling is
  untouched, so parametric draws are bit-identical to a `store = "fits"` fit
  with the same seed. `store = "fits"` remains the default (no silent change);
  a `store = "fits"` fit retaining more than a gigabyte of per-draw BART blocks
  notes the alternative once per session. A size-adaptive default was
  considered and rejected: storage semantics should not depend on data size.
  See `docs/plans/bart-train-recompute.md`.

## Deprecated

* The NUTS-specific `stan_args` controls - `adapt_gamma`,
  `adapt_kappa`, `adapt_t0`, `adapt_init_buffer`, `adapt_term_buffer`,
  `adapt_window`, `stepsize`, `stepsize_jitter`, and `max_treedepth` - have
  no analog under WALNUTS. They are still accepted (a script that sets one
  does not break), but are now ignored, and a warning naming every supplied
  deprecated argument is issued at fit time. They will be removed in a
  future release. `init_r` (the initial-position radius), `adapt_delta` (see
  the Sampler section), and the loop-level arguments (`iter`, `warmup`,
  `skip`, `chains`, `cores`, `refresh`, `seed`, `verbose`) keep a live
  meaning and do not warn.

## Removed

* Sampler diagnostics tied to NUTS - divergent transitions, max-treedepth
  transitions, and low E-BFMI, previously warned on by
  `check_sampler_diagnostics` - are gone. The vendored WALNUTS sampler
  reports no analog of any of the three (it exposes only position, log
  density, step size, and the diagonal mass estimate), so the corresponding
  warnings no longer fire.

## New features

* `bart_args` now reaches the whole model-level half of a `dbarts`
  specification, not just `dbartsControl`'s formals plus a hand-wired
  `k`/`power`/`base`/`split.probs`. The BART component's control/model/data
  triple is built by `dbarts::dbartsSpec()`, exported by dbarts 1.0-0 for
  exactly this purpose, in place of a hand-assembled `dbartsModel` and a
  `dbarts:::parsePriors` call through a `:::` shim, so `tree.prior`,
  `node.prior`, `proposal.probs`, `monotone`, `interactions()`, `blocks()`,
  and `seed` all pass through and are validated by dbarts itself. Priors are
  resolved in dbarts's own vocabulary, so `node.prior = normal(k = chi(1.25,
  Inf))` works whether or not dbarts is attached, and `k`/`power`/`base`/
  `split.probs` remain as shorthand for the priors they write into (giving
  both a shorthand and its prior is now an error). Draws for every previously
  expressible `bart_args` are bit-identical.

  `sigma`, `resid.prior`, `resid.dist`, and `variance` are refused: the
  parametric component draws the residual standard deviation and the BART
  component conditions on it each sweep, so the forest has no residual model
  of its own to configure. `family` is refused unless it names the one the
  response already implies.

* New `mvbart()` jointly fits BART to two or more continuous outcomes that
  share a predictor set but have correlated residuals - the seemingly-
  unrelated-regressions (SUR) analogue of BART. Each outcome gets its own
  sum-of-trees mean surface; the surfaces are coupled through a residual
  covariance matrix with a conjugate inverse-Wishart prior, sampled jointly,
  which calibrates the joint predictive distribution better than one
  independent BART per outcome. Composed from existing machinery - no engine
  support was required. Results are `mvbartFit` objects with a `print`
  method; see `?mvbart`.

* Tracking dbarts 1.0-0 retiring `resid.prior` onto the family object: the
  BART component's fixed unit residual variance is now set with `family =
  gaussian(sigma = fixed(1))` instead of the now-tombstoned `resid.prior =
  fixed(1)`, in both `stan4bart_fit` and `mvbart`; no observable change.

## Breaking changes

* Factor variables in the `bart()` part of the formula are now encoded with
  `dbarts`'s categorical splits (`factors = "categorical"`) instead of one
  indicator column per level (`factors = "indicators"`, the previous, port-era
  default). A `bart()` factor is now a single design column whose tree prior
  chooses level subsets to split on, rather than treating each level as an
  independent 0/1 predictor - a different prior over factor structure, so
  draws move for any fit with a factor in the `bart()` part; factor-free fits
  are bit-identical. `extract(type = "varcount")` now has one row per factor
  (named by the factor's variable name) instead of one row per level. An
  unseen `bart()` factor level in `newdata` or `test` is still an error, as
  before. Fixed-effect and random-effect factor handling - `model.matrix`
  contrasts and the `lme4` grouping-factor machinery, including its new-level
  semantics under `sample_new_levels` - are unchanged by this change.

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

* Fixed `bart_args` silently dropping any name that did not match a
  `dbarts::dbartsControl`/`dbarts::dbartsSpec` formal, rather than erroring:
  most consequentially, `dbarts::dbartsControl`'s `rngSeed` argument was
  renamed to `seed`, so `bart_args = list(rngSeed = ...)` reached no formal
  and the seed was silently ignored. `rngSeed` is now mapped to `seed` (with
  a once-per-session warning naming the rename); any other unmatched name
  now errors, naming the offending argument, instead of being dropped.

* `mvbart`'s `bart_args` dropped unmatched names in the same silent way, and
  now refuses them by name too. The rng seed is `mvbart`'s own argument
  there, since the seeds are derived per chain and per equation, so neither
  `seed` nor `rngSeed` is accepted in `bart_args`.

* Fixed `bart_args = list(n.cuts = ...)` corrupting the `dbartsControl`
  object when given a non-integer numeric: it was written into the `n.cuts`
  slot with `attr<-`, which bypasses the coercion `dbartsControl()` itself
  applies, so a plain numeric `n.cuts` landed in the slot uncoerced and
  `validObject` later failed ("should be or extend class integer"). `n.cuts`
  is already picked up, coerced, and validated by the existing
  `dbartsControl()` construction call; the separate `attr<-` assignment is
  removed.

* All-zero `weights` are now refused up front with an error naming the
  cause. Previously they failed inside the `lme4`/`glm` fits that seed the
  sampler, with messages such as "object 'fit' not found" (binary) or
  "sigma_init is NaN" (continuous).

* Fixed the last-resort initialization fit building the formula `y` instead
  of `y ~ 1` when the model has only `bart()` and random-effect terms, which
  failed with "object 'y' not found".

* The random-effect standard deviation is now drawn once per sweep by an exact
  slice move along the curve that leaves the linear predictor fixed, and this is
  on by default. Rescaling a grouping factor's standard deviation and dividing
  that factor's standardized effects by the same amount leaves the random
  effects, the linear predictor and the likelihood untouched, so the conditional
  along that curve is one-dimensional, log-concave and free of any data term - a
  slice sampler draws it exactly, at a cost that does not register against a BART
  sweep. That curve is the one direction the gradient-based sampler could not
  travel, and it was the whole of the random-effect scale's mixing problem. On a
  20-group, 100-observations-per-group Gaussian random intercept with the default
  `skip`, the scale's lag-1 autocorrelation drops from 0.96 to 0.09 and its ESS
  rises from 7 to 526 per 1000 kept draws, worst of five seeds, at a wall-time
  ratio of 1.01 on a four-chain fit with a counterfactual test surface. Across
  four grouped designs at the package defaults the group standard deviation now
  clears the dbarts hand-off bar - lag-1 below 0.8 and ESS at least 100 per 1000
  - in eight of twelve fits, where before it cleared it in none. Posterior means
  are unmoved within Monte Carlo error. This is a draw-moving change for every
  model with a random-effect block; `stan_args = list(ridge_move = FALSE)`
  restores the previous behavior, which is how the comparison above was
  measured. `skip` remains the escape hatch for what the move does not reach -
  notably the residual standard deviation, which pays the BART-versus-parametric
  alternation whether or not there is a random-effect block.

* Fixed the `"stan"` element of `skip` reaching the sampler and being consumed
  by nothing: the parametric block took one transition per BART sweep whatever
  it was set to, while the `"bart"` element had been honored throughout. It now
  takes that many parametric transitions per sweep during sampling and keeps the
  last, at no change to the number of draws returned. This is the available
  remedy for a poorly mixing random-effect standard deviation: the non-centered
  block ties the log scale to the standardized effects on a product ridge, which
  a diagonal metric cannot rescale, so that one direction of the parametric
  block mixes an order of magnitude slower than the rest. On a 20-group,
  100-observations-per-group Gaussian random intercept fit,
  `skip = c(bart = 1, stan = 8)` takes the scale's lag-1 autocorrelation from
  about 0.96 to about 0.77 and its ESS from about 14 to about 115 per 1000 kept
  draws at 1.4x wall time, and `stan = 16` to about 0.63 and 250 at 2.0x. Warmup
  keeps one parametric transition per sweep at any `skip`, so that the
  parametric block cannot outrun the still-growing forest and park the response
  mean in the random intercepts. Draws at the default `skip = 1` are unchanged,
  so the remedy has to be asked for.

* Fixed `dbarts_results.structSize` never being set by the dbarts 1.0
  flat-C-API port, which caused the versioned-struct field gate to skip
  populating every run's output buffers (usually masked by the buffers
  otherwise reading as zero). Found and fixed while recording this arc's
  Stan-era baselines.

* Fixed `extract`/`fitted` returning a dimensionless result for `type = "ev"`
  and `"ppd"` when the fit was made with `offset_type = "bart"`. Composing
  those types reads the result's shape from the BART component, but that
  `offset_type` substitutes the length-n offset vector for it, leaving the
  composed draws without dimensions.

* Fixed `extract(type = "ev")` and `type = "ppd"` silently ignoring
  `offset_type = "fixef"` and `"ranef"`: the stored values were compared
  against the internally-checked spellings `"fixed"` and `"random"`, so
  neither ever matched and the offset was dropped from the composition.

* Fixed three defects in drawing random effects for brand-new grouping
  levels under `sample_new_levels = TRUE`, all of which required a grouping
  factor with more than one predictor (a random slope); intercept-only
  factors were unaffected. The per-draw covariance was sliced with
  `drop = FALSE`, leaving a 3-D array where `chol()` requires a square
  matrix, so the path errored outright; the standard normal draws were
  allocated with `n_predictors^2` rows instead of `n_predictors`; and the
  Cholesky factor was both transposed at construction and applied with
  `crossprod`, so the draws carried covariance `R R'` rather than the
  target `Sigma` (equal only when `Sigma` is diagonal).

* Fixed `predict` capturing a positionally-supplied `type` as `offset`:
  `predict(fit, newdata, "ev")` either errored on a non-numeric addition or,
  for the `indiv.*` types, quietly returned the default `"ev"` component.

* Fixed `as.matrix` collapsing the whole draw array into a single column
  instead of reshaping it to (iterations x chains) by parameters.

* Fixed `as.array`'s `include_warmup` silently returning the post-warmup
  draws: the internal warmup-expression rewriter recognized only an object
  literally named `object`, and `as.array`'s argument is named `x`.

* Fixed a spurious length-recycling warning whenever a grouping factor
  gained new levels, from comparing random-effect row names with `==`
  rather than `identical()`.

* Fixed a binary response's 0/1 `weights` being silently dropped instead of
  installed as the BART component's active-row mask. `dbarts` 1.0-0 changed
  what 0/1 weights mean for a probit fit: they no longer weight the
  likelihood (a weighted probit has no tractable latent-variable form) and
  instead name the rows in and out of it, resolved by `dbarts::dbartsSpec`
  into a mask the caller must install on the sampler it builds. `stan4bart`
  already called `dbartsSpec` but discarded that mask, so every row re-entered
  the likelihood regardless of `weights`. The mask is now installed on every
  sampler a binary fit builds, including the one rebuilt after a
  `saveRDS`/`readRDS` round trip; a non-0/1 weight vector on a binary response
  is refused with `dbarts`'s own message, as before.
