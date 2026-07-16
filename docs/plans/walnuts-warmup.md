# walnuts-warmup

agent: opus (a SAMPLER-QUALITY + runtime arc, not a correctness rewrite. It
  tunes the WALNUTS warmup path - the initial step-size/mass SEED and the warmup
  ALLOCATION - to cut the leapfrog-step count per transition, which perf-review.md
  identified as the single largest stan4bart-owned runtime lever. The parametric
  target (src/parametric_model.hpp), its hand adjoint, and the FD gate are FROZEN;
  this arc changes only how WALNUTS is initialized and driven, plus docs and a
  diagnostic. The BART side and dbarts.h are untouched.)

rng/oracle: THERE IS NO BITWISE ORACLE. Two draw-classes per commit:
  (1) DRAW-NEUTRAL - the measurement harness, an eval-call counter, docs, and the
  surfaced diagnostic touch nothing on the sampling/warmup path; gate is the full
  testthat suite green + `identical()` draws vs the pre-commit build.
  (2) DRAW-MOVING - anything that changes the warmup trajectory (the seed) or the
  warmup length (allocation) moves EVERY subsequent draw, because warmup consumes
  RNG and fixes the frozen tuning. These are NOT bitwise-gated; they are gated by
  the pre-registered DISTRIBUTIONAL-EQUIVALENCE tiers vs the Stan-era baseline
  (benchmarks/baselines/posterior-baselines-75d7970.rds): means within
  4*sqrt(mcse_stan^2 + mcse_walnuts^2), sd ratio in [0.8, 1.25], every
  fixef/ranef/sigma/Sigma summary, all five tiers (compare-posterior.R). The
  Stan-era baseline is the FIXED equivalence oracle and is NOT re-blessed: a seed
  or warmup change keeps targeting the SAME posterior, so it must still clear the
  SAME band. What legitimately regenerates is the value-pinned testthat snapshots
  (test-01 etc.), which move because the sampler tuning path changed. Plus ESS/sec
  as a SECOND floor (below): a change that speeds iterations but hurts mixing is a
  net loss and is rejected even if the distributional band passes.

budget: est. ~300-600 lines net, concentrated in benchmarks/R (the measurement
  harness, the largest piece) and a ~20-line change to WalnutsSampler's ctor (the
  seed), plus small R-surface docs/defaults and a mutable eval counter in
  ParametricModel. No new .Call signature is strictly required (getAdaptationInfo
  already exists; the leapfrog diagnostic reuses that seam). Chiefly NEW:
  benchmarks/R/bench-warmup.R (per-iter wall + ESS/sec + warmup-grid over a
  large-n-RE tier). EDITED: src/walnuts_sampler.cpp (the seed), src/parametric_model.hpp
  (the eval counter), R/stan4bart.R / R/stan4bart_fit.R (a documented floor +
  the diagnostic), man/ + NEWS.

window: the user-facing surface is (a) whether the warmup DEFAULT or a floor
  changes (a draw- and runtime-affecting default flip - every user), and (b) a new
  cheap post-fit diagnostic. Both are stan4bart-only R-surface decisions; no dbarts
  lockstep, no dbarts.h change.

## Goal

Cut the WALNUTS parametric step's per-iteration cost at large n with random
effects by making the FROZEN tuning good with LESS warmup and from a better start,
and protect users who under-warm. perf-review.md #1 measured the mechanism: at
n=1e5 continuous with q=1000 random-intercept levels, per-iteration sampling wall
is 69 ms/iter after a 20-iteration warmup but 28 ms/iter after warmup >= 100 (2.2x),
and the driver is the leapfrog-step count per transition, set by the frozen
step-size/mass tuning - monotone-improving in warmup, non-monotone in q. The
parametric step is 27-65% of per-iter wall depending on tuning quality. The arc
must MEASURE each candidate's payoff (per-iter wall AND ESS/sec) on the reference
tiers before shipping it, and must clear the distributional-equivalence floor -
this is a sampler-quality change, not a numerically-inert one.

## Verdict up front (the thesis, to be confirmed by C0 measurement)

The SEED and the warmup ALLOCATION are ONE coupled decision, not two. At the
production default (iter=2000, warmup=iter/2=1000, stan4bart.R:13-14) the tuning
is already well past the warmup>=100 knee, so a better seed barely speeds the
SAMPLING phase there. The seed's real payoff is twofold: (i) it lets the frozen
tuning converge in fewer warmup transitions, which is what would make a SMALLER or
dimension-aware warmup default SAFE; and (ii) it protects the user who explicitly
under-warms (short iter, or warmup set low), moving them off the 69 ms end of the
curve. So "improve the seed" and "should warmup allocation change" are the same
lever seen from two ends. Both payoffs are currently UNMEASURED (perf-review
measured the warmup-length SENSITIVITY, not whether a step-probe + Nutpie seed
reaches the knee faster). Therefore the first commit is the MEASUREMENT, and the
seed/allocation changes are conditional on what it shows. The honest floor is
ESS/sec: fewer leapfrog steps per iteration is only a win if it does not cost
proportionally more iterations to the same effective sample size.

## Context (file:line ground truth, read in code)

THE CURRENT SEED (src/walnuts_sampler.cpp). Confirmed: the sampler starts unit
mass, step size 0.1, and uses NEITHER upstream seed heuristic.
  - `impl_->step_size = 0.1` (walnuts_sampler.cpp:213), fed to the Adam optimizer
    as log(0.1) - a hardcoded constant, not probed.
  - initial mass = `Eigen::VectorXd::Ones(dim)` (walnuts_sampler.cpp:277-278),
    i.e. an identity diagonal metric.
  - the ctor builds `InitChainConfig` DIRECTLY (single chain); it never calls
    `InitConfigBuilder`, so the two builder-only seed heuristics below are unused.
  - `init_r` (Stan's random-start radius) defaults to 2.0 (stan4bart_fit.R:547);
    the initial position is uniform(-2, 2) per coordinate (walnuts_sampler.cpp:240),
    a wide random start. Any seed heuristic evaluated AT this point sees a random,
    possibly poorly-conditioned location (esp. the nc>=3 onion z_T coords).

THE WARMUP DRIVE (R/stan4bart_fit.R, src/init.cpp). One WALNUTS transition per
outer Gibbs sweep, driven manually.
  - defaults: `iter = 2000L`, `warmup = iter %/% 2L` (= 1000) (stan4bart.R:13-14).
  - the worker runs `C_stan4bart_run(sampler, warmup, TRUE, ...)` then
    `disengageAdaptation` (-> `freeze()`) then `run(iter - warmup, FALSE)`
    (stan4bart_fit.R:51-75). `run()` loops numIter times calling
    `paramSampler->run(isWarmup)` (init.cpp:583-591); `run(true)` calls the
    AdaptiveWalnuts functor, `run(false)` the frozen WalnutsSampler.

CRITICAL: `WarmupConfig::min_max_iter` IS INERT IN THIS DRIVE. The wrapper sizes
`WarmupConfig().min_max_iter(num_warmup, num_warmup)` (walnuts_sampler.cpp:244-247),
but `AdaptiveWalnuts::operator()` (adaptive_walnuts.hpp:234-251) reads NEITHER
min_iter nor max_iter nor the convergence tolerances - those live only in the
THREADED `detail::adapt`/`controller_loop`/`AdaptWorker` path (adapt.hpp:110-259,
the high-level api.hpp driver), which stan4bart does not use. The manual per-sweep
drive's warmup LENGTH is purely how many times `run(true)` is called (= warmup).
So the `num_warmup` argument threaded into the ctor currently affects nothing
observable, and the design note's claim that WarmupConfig "is sized so the mass
windows schedule correctly" (docs/design/walnuts.md:118-122) is imprecise: the
mass estimator's discount uses `mass_init_count + iteration` (adaptive_walnuts.hpp:76),
not max_iter, and Adam decays on its own `t` counter (adam.hpp:80), not a horizon.

WHAT UPSTREAM EXPOSES WITHOUT PATCHING (vendored at 5854be8; the seed levers).
  - `walnuts::detail::adapt_step(rng, logp_grad, theta, M, step, D)` (util.hpp:287-304):
    a NUTS-style step-size probe - doubles the step while the one-leapfrog
    acceptance error > log(0.9), then multiplies by sqrt(0.5) while < log(0.6).
    A handful of `leapfrog_error` calls, each 2 gradient evals. It is a free
    function in `walnuts::detail`, directly callable; also wrapped by
    `InitConfigBuilder::adapt_step_build` (config.hpp:465-472).
  - `InitConfigBuilder::masses(logp_grad, mass_smoothing, average_masses=false)`
    (config.hpp:356-378): the Nutpie gradient mass seed,
    `mass = (1 - mass_smoothing)*|grad| + mass_smoothing`, ONE gradient eval at
    the init position. (The wrapper would call the builder, or replicate this
    one-liner directly against its `ParametricModel`.)
  - the `MassEstimator` seeds BOTH its draw-variance and score-variance online
    estimators from the init mass, regularized by pseudo-count
    `mass_init_count` (default 4.0, adaptive_walnuts.hpp:54-62). So the init mass
    biases the whole warmup mass estimate, most strongly under SHORT warmup - a
    better seed helps beyond t=0.

WARMUPCONFIG KNOBS THAT ARE LIVE IN THE MANUAL DRIVE (builder-settable, no patch):
  - `mass_init_count` (4.0): the regularization pseudo-count toward the init mass.
  - `max_macro_steps_target` (15.0, config.hpp:628): the MinMicroStepsAdaptHandler
    target; it raises the minimum micro-steps-per-macro-step to hold expected
    macro steps near this target, directly shaping leapfrog count vs tree-doubling.
  - `step_accept_rate_target` (0.8): the Adam target acceptance rate - the WALNUTS
    ANALOG of Stan's adapt_delta. A lower target yields a larger step and fewer
    leapfrogs (a direct speed/mixing knob).
  - `step_learning_rate` (0.05), `step_learn_rate_decay` (0.5): Adam convergence
    speed; faster convergence -> good tuning in fewer warmup transitions.
SamplingConfig (frozen phase, sampler.hpp): `max_trajectory_doublings` (5, caps
macro steps at 32), `max_step_halvings` (5), `max_hamiltonian_error` (0.5),
`min_micro_steps` (adaptation-overridden). All at defaults today.

OBSERVABILITY (what the fit already carries, sample-storage.md).
  - `getAdaptationInfo` (init.cpp:838-859) surfaces the frozen step size and
    diagonal inverse mass per chain, landed as `fit$adaptation$step_size` /
    `fit$adaptation$inv_mass` (stan4bart.R:443-476); plus `fit$adaptation$snapshot`
    (warmup-end position) and `fit$adaptation$trace` (a thinned warmup trace of
    the monitored scalars incl. lp__ and stepsize__, stan4bart_fit.R:56-66).
  - the MEAN LEAPFROG / EVAL COUNT per transition - the quantity that IS the
    runtime - is NOT surfaced: `n_leapfrog__` is a constant-zero placeholder row
    (walnuts_sampler.cpp:301), and the `LatestDraw` handler captures only
    position/lp/step_size/inv_mass (the ChainHandler `on_warmup`/`on_sample`
    signatures carry no trajectory depth). AdaptiveWalnuts observes `1 << depth`
    internally (adaptive_walnuts.hpp:248) but never hands it out.
  - CHEAP DRAW-NEUTRAL ROUTE to it: `ParametricModel::eval` is called exactly once
    per leapfrog micro-step (the hot path). A `mutable` call counter on the model
    (eval is const, reached through `std::cref`, parametric_model.hpp:134-136) gives
    exact evals; mean leapfrog per iter = evals-during-sampling / sampling-transitions.
    No vendoring, no draw change.

THE HARNESS THAT EXISTS (benchmarks/R). Reused, not reinvented.
  - `bench-perf.R`: per-iteration sampling wall by DIFFERENCING two runs at matched
    warmup, different iter (slope = per-iter cost); plus install time and peak RSS.
    Reference fits: continuous_nc2, binary_multilevel, weighted_continuous. It has
    NO large-n-RE tier and NO ESS/sec - both must be added.
  - `compare-posterior.R`: the distributional gate vs posterior-baselines-75d7970.rds,
    k=4 means, sd ratio [0.8,1.25], tiers continuous_nc1/nc2/nc3/weighted/binary;
    reuses fit_tier/summarize_fit/chain_health/batch_means_summary from
    record-posterior-baselines.R. This IS the correctness floor for every
    draw-moving commit, run as-is.

## Design

### The measurement, first (Q4 - because the payoff is unmeasured)

Add `benchmarks/R/bench-warmup.R`, the instrument that judges every candidate. It
must report, per candidate x tier, BOTH:
  - per-iteration SAMPLING wall (bench-perf.R's matched-warmup differencing), AND
  - ESS/sec: bulk and tail ESS (posterior::ess_bulk / ess_tail) of the monitored
    scalars - the fixef, sigma, and the log-scale group SDs derived from Sigma
    (sample-storage.md's monitored set, the honest slow-mixing canaries) - divided
    by the sampling wall. ESS/sec is the HONEST metric: a tuning change that halves
    the leapfrog count but halves ESS is not a win.
  - the mean leapfrog / eval count per iteration (from the new eval counter), which
    is the mechanistic quantity the wall tracks.
Tiers: the five distributional tiers PLUS a new large-n-RE tier reproducing
perf-review's case (n=1e5 continuous, q=1000 random-intercept levels) - this is
where the 2.2x lives and where the seed must prove out. Over that tier, sweep
warmup in {20, 50, 100, 200, 400, 1000} to trace the 69->28 ms knee and locate
where each candidate seed lands on it (does a better seed shift the knee LEFT?).
This tier is expensive; run it on the quiet bench box (memory: dbarts x86 bench
box, `ssh dbarts-bench`) and never concurrently with other load. The correctness
floor for any change is compare-posterior.R unchanged.

### The seed (Q1)

Options (all draw-moving; each moves the whole warmup trajectory):
  - A1 STEP PROBE ONLY: call `detail::adapt_step` at the init position before
    constructing the AdaptiveWalnuts, replacing the constant 0.1 with a probed
    step. ~10-30 gradient evals once per chain (negligible vs a full fit). Gives
    Adam a start near the target so it converges in fewer warmup transitions.
  - A2 NUTPIE MASS ONLY: seed the init mass with `(1-s)*|grad| + s` (one gradient
    eval) instead of ones. Improves conditioning immediately and, via
    mass_init_count=4, biases the whole warmup mass estimate - the larger help
    under short warmup.
  - A3 BOTH (the natural pairing): seed the mass first, then probe the step USING
    that mass (adapt_step takes M). This is what upstream's `adapt_step_build`
    after `masses(...)` does, and matches Nutpie/Stan practice. RECOMMENDED shape
    IF C0 shows payoff.
  - A4 STRUCTURAL MASS from the parametric block instead of the gradient: seed the
    diagonal from the prior scales (prior_scale for beta, the decov scale for the
    RE block) and the response scale, i.e. a curvature guess that needs NO gradient
    eval and cannot be poisoned by a bad random init point. Weaker than the true
    gradient curvature but robust. Best role: the FALLBACK when the gradient at the
    random init is non-finite or ill-conditioned (the onion dot_self-near-zero
    risk the design note already flags).
  - orthogonal refinement: evaluate the seed at a CENTERED point (zero / prior
    mean) rather than the uniform(-2,2) init draw, to cut seed noise. Cheap, and it
    de-randomizes the seed (same seed still reproduces via the position draw that
    follows).
Non-finite guard: the seed evals must use the same isfinite/allFinite contract
eval already enforces (parametric_model.hpp:558-564); on a non-finite gradient at
the init point, fall back to A4 (structural) or unit mass + step 0.1, never seed
from a poisoned value.
RECOMMENDATION: A3 (Nutpie mass + step probe) with the A4 structural fallback and
the centered-eval refinement - but LAND IT ONLY IF C0 shows it reaches the
warmup>=100 tuning quality materially sooner (knee shifts left) without an ESS/sec
regression. If C0 shows the seed barely moves the knee, stop at docs (Q3) and do
not ship a draw-moving change for no measured gain.

### Warmup allocation (Q2)

Options (B2/B3 draw-moving and runtime-changing for every user):
  - B1 KEEP iter/2, change nothing here; rely on the seed (Q1) + docs/diagnostic
    (Q3). Conservative: no user's draws or runtime shift from allocation.
  - B2 A DIMENSION-AWARE FLOOR: `warmup = max(user_warmup, floor(dim))` where dim
    scales with K+q (the parametric dimension) or with n/q. perf-review shows
    warmup>=100 suffices even at n=1e5/q=1000, so the floor is modest. Protects
    under-warmers but RAISES warmup for anyone who set it below the floor - each
    added warmup iter is a full BART sweep (~28 ms at n=1e5), and it moves their
    draws. The floor's FORM and VALUE must come from the C0 warmup-grid, not a
    guess.
  - B3 LOWER THE DEFAULT below iter/2 (with a good seed making it safe) - e.g. a
    convergence-sized default. Largest potential runtime win (warmup is half the
    iterations by default; at n=1e5 that is ~14 s of the fit), but the riskiest:
    under-warming silently degrades every fit unless the seed genuinely closes the
    gap, and it changes every user's draws.
RECOMMENDATION: B1 as the default resolution (no allocation change), PLUS a
documented floor and the diagnostic warning (Q3) rather than an enforced floor -
unless C0's warmup grid shows a specific dimension-aware floor buys ESS/sec
broadly across tiers, in which case B2 with the measured form. Reject B3 unless the
seed measurement is strong; a silent default reduction that under-warms is the
worst failure mode (it degrades correctness invisibly).

### Docs + the diagnostic (Q3, draw-neutral)

  - a documented WARMUP FLOOR / guidance for large-n-RE fits in ?stan4bart and the
    vignette: "for large n with many random-effect levels, warmup below ~100 leaves
    WALNUTS mis-tuned and can run 2x slower per iteration; the default iter/2 is
    ample, but if you shorten iter keep warmup >= <measured floor>." Pure docs, no
    behavior change.
  - the LEAPFROG DIAGNOSTIC: surface mean leapfrog / eval per sampling iteration
    (from the eval counter, via the getAdaptationInfo seam) as
    `fit$adaptation$mean_leapfrog` (per chain), and have print/summary emit a
    one-line warning when it is high relative to dimension (a heuristic threshold
    from C0, e.g. mean leapfrog per iter >> the well-warmed tier value) -
    "parametric sampler tuning looks poor; consider more warmup." This is the
    honest, cheap, sampler-agnostic canary for the geometry-mismatch risk
    sample-storage.md flags (tuning freezes against warmup-era f, then f keeps
    moving, and WALNUTS exposes no divergence/E-BFMI hook). It reuses the existing
    adaptation slot, needs no new .Call, and is draw-neutral.

### RNG / compatibility classes per commit (Q5)

  - DRAW-NEUTRAL (suite green + `identical()` draws; no tier gate):
    the bench-warmup.R harness, the eval counter + its getter, the surfaced
    diagnostic, all docs/NEWS. None touches the sampling or warmup code path.
  - DRAW-MOVING (distributional tiers + ESS/sec floor; value-pinned snapshots
    regenerated; Stan baseline NOT re-blessed):
    the seed change (C1) and any allocation change (C2). Each re-runs
    compare-posterior.R on all five tiers (must stay in band vs the FIXED
    posterior-baselines-75d7970.rds - same posterior, so it must still pass) and
    bench-warmup.R for ESS/sec (must not regress). test-01/test-02 etc. value pins
    are regenerated because the tuning path moved, exactly as the C2/C3 sampler
    landings did (docs/plans/walnuts.md).

## Commits (each installs and tests; per-commit gates)

C0. MEASUREMENT HARNESS + THE INSTRUMENT (draw-neutral). Add bench-warmup.R
   (per-iter wall + ESS/sec + mean-leapfrog, over the five tiers plus a large-n-RE
   tier with the warmup grid {20,50,100,200,400,1000}). Add the `mutable` eval
   counter to ParametricModel and a `getMeanLeapfrog`-style accessor on the
   sampler surface (reuse the getAdaptationInfo path). RECORD the baseline curve:
   reproduce perf-review's 69->28 ms knee at n=1e5/q=1000, and record ESS/sec and
   mean-leapfrog per tier for the CURRENT unit-mass/step-0.1 seed. Files:
   benchmarks/R/bench-warmup.R (new), src/parametric_model.hpp (counter),
   src/walnuts_sampler.{hpp,cpp} + src/init.cpp (accessor, if the diagnostic is
   wired now). Gate: full testthat green; draws `identical()` to the pre-commit
   build (the counter is inert); the harness reproduces the measured 2.2x. Size:
   M. Abort: if the harness cannot reproduce the knee, the premise is wrong -
   stop and reassess before changing the seed.

C1. THE SEED (draw-moving), CONDITIONAL ON C0. Replace the constant unit
   mass/step 0.1 with the recommended A3 seed (Nutpie mass + probed step) evaluated
   at a centered point, with the A4/unit fallback on non-finite gradient. Files:
   src/walnuts_sampler.cpp (the ctor). Gate: compare-posterior.R in-band on all
   five tiers vs posterior-baselines-75d7970.rds; bench-warmup.R shows the knee
   shift LEFT and no ESS/sec regression on any tier; test-05-rng reproducibility
   green (same seed -> identical draws); value-pinned snapshots regenerated.
   Size: S-M code, L gate. Abort/DO-NOT-LAND: if the knee does not move or ESS/sec
   regresses - a draw-moving change with no measured runtime gain is not worth the
   snapshot churn; fall back to docs-only.

C2. WARMUP ALLOCATION (draw-moving), CONDITIONAL ON C0's grid. Default B1 (no
   change); land B2's dimension-aware floor ONLY if the grid justifies a specific
   form. Files: R/stan4bart.R (default/floor), R/stan4bart_fit.R. Gate: as C1 plus
   a runtime accounting (the floor's added warmup cost is less than the sampling
   speedup it buys). Size: S code, L gate. Abort: reject any silent default
   reduction (B3) that under-warms a tier below its measured floor.

C3. DOCS + DIAGNOSTIC SURFACE + NEWS (draw-neutral). The documented warmup-floor
   guidance, the `fit$adaptation$mean_leapfrog` field, the print/summary warning
   at the C0-measured threshold, NEWS. Files: man/*.Rd, R/generics.R (print/summary),
   R/stan4bart.R (assembly), NEWS.md. Gate: R CMD check man; suite green; draws
   `identical()`. Size: S-M.

## Open questions for VD

Q(a) SHIP THE SEED, OR STOP AT DOCS. Fork: whether to land a draw-moving seed
  change at all. The upstream heuristics (adapt_step + Nutpie mass) are free to
  call and are the standard warm-start, but their payoff HERE is unmeasured, and at
  the production default (warmup=1000) the tuning is already past the knee, so the
  seed only helps the warmup PHASE and under-warmers. Options: (i) land A3 (Nutpie
  mass + probed step) IF C0 shows the knee shifts left without an ESS/sec
  regression (RECOMMENDED - it is cheap, principled, and unlocks a smaller safe
  warmup); (ii) mass-only A2 (one gradient eval, less code, smaller draw move);
  (iii) docs/diagnostic only, no seed change (zero draw churn, zero snapshot
  regeneration, but leaves the under-warmer on the 69 ms end). Cost of (i)/(ii):
  regenerate the value-pinned testthat snapshots and re-run all five tiers; a
  one-time gate. RECOMMENDATION: (i), gated strictly on C0 - if the measurement is
  flat, drop to (iii). The plan is built so C0 decides this, not assertion.

Q(b) WARMUP DEFAULT - LEAVE iter/2, ADD A FLOOR, OR LOWER IT. Fork: changing
  allocation changes EVERY user's draws AND runtime. Options: (i) keep iter/2 and
  rely on the seed + a DOCUMENTED (not enforced) floor + the diagnostic warning
  (RECOMMENDED - conservative; no silent behavior change; the diagnostic tells a
  user when they under-warmed); (ii) an ENFORCED dimension-aware floor
  `warmup = max(user, f(dim))` (protects under-warmers automatically, but silently
  raises warmup and cost for anyone below it, and moves their draws - the form must
  come from C0); (iii) LOWER the default below iter/2 with the seed making it safe
  (biggest runtime win - warmup is half the iterations - but under-warming degrades
  correctness invisibly if the seed does not close the gap). Cost of (i): a user
  who ignores the docs and warning still pays 2x. Cost of (ii): a surprising
  runtime/draw change on fits that set warmup low deliberately. Cost of (iii):
  correctness risk. RECOMMENDATION: (i); escalate to (ii) only if C0's grid shows a
  clean dimension-aware floor that buys ESS/sec broadly; never (iii) without a
  strong seed result. The perf-review's own knee (warmup>=100 at n=1e5/q=1000)
  suggests any floor is modest.

Q(c) THE LEAPFROG DIAGNOSTIC - EXACT COUNTER, DERIVED PROXY, OR PATCH WALNUTS.
  Fork: how to surface "tuning looks bad." Options: (i) a `mutable` eval-call
  counter in ParametricModel -> mean leapfrog per iter, exact, draw-neutral, no
  vendoring, ~10 lines (RECOMMENDED); (ii) derive a weaker proxy from the existing
  stepsize trace (a still-drifting stepsize at warmup end implies under-convergence)
  - free, uses only stored data, but indirect and easy to misread; (iii) patch the
  vendored walnuts ChainHandler/transition to surface true per-transition trajectory
  depth / n_leapfrog - the most faithful, but re-paid on every upstream refresh
  (we track 5854be8) and against the arc's no-patch stance. Cost of (i): a tiny
  always-on increment on the hot eval (a single `++` - measured-negligible, but
  worth confirming it does not perturb the __restrict codegen that perf-review #2
  found sensitive). RECOMMENDATION: (i); it is exact, local, and reuses the
  adaptation slot. Confirm in C0 that the counter increment is codegen-inert
  (draws `identical()`, per-iter wall unchanged) - if it is not, fall back to (ii).

Q(d) EXPOSE WALNUTS TUNING KNOBS TO USERS. Fork: WALNUTS' live knobs
  (step_accept_rate_target - the adapt_delta analog; the initial step size;
  max_macro_steps_target) are builder-settable and could be surfaced as stan_args
  for power users who want a manual speed/mixing lever. Options: (i) keep them all
  internal, fixed at upstream defaults tuned only by this arc (RECOMMENDED - keeps
  the surface small; the C4/C5 arc just DEPRECATED the old NUTS knobs incl.
  adapt_delta, and resurrecting an accept-rate knob under a new name reopens that
  surface and its support/deprecation debt); (ii) expose ONLY an initial step size
  as a stan_arg (harmless override, useful if C0 finds a case the auto-probe
  mistunes); (iii) wire step_accept_rate_target to the already-ignored adapt_delta
  for source-compat (rejected - it un-deprecates a knob we just retired and blurs
  the WALNUTS/NUTS semantics). Cost of (i): a user with a pathological target has no
  escape hatch but more warmup. RECOMMENDATION: (i), with (ii) held in reserve only
  if C0 surfaces a concrete mistuning the auto-seed cannot fix. Do NOT do (iii).
