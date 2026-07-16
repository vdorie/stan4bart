# bart-train-recompute

scope: IMPLEMENTATION PLAN. Replace, or make optional, the stored
  `bart_train` block (n x draws x chains doubles; 94-98% of a large fit
  object once warmup is default-off) with on-demand recomputation from the
  kept trees through dbarts's predict path. The BART engine's flat C API
  (dbarts.h) is FROZEN and append-only; the predict / tree-storage entries
  this arc needs (predictBART, getTrees, setTreeStorage, numSavedSamples,
  exportBARTState/createStoredBARTSampler) ALREADY EXIST and are bound in
  src/init.cpp. Approved in principle 2026-07-16 "as long as it is reasonable
  work that fits within our general speed/flexibility objective"; whether it
  lands as a default flip, an opt-in, or a size-adaptive default is decided
  here by measurement.

budget: est. M-L, mostly R-side, plus ONE dbarts ENGINE bugfix (not a dbarts.h
  ABI change - see C1/Q2). C0 harness S; C1 dbarts setState-at-scale fix M and
  the risky item (engine investigation); C2 serialization retention S-M (R +
  reuse of existing C entries); C3 recompute seam + opt-in M (R in generics.R,
  stan4bart.R, stan4bart_fit.R); C4 default-policy decision S (policy + docs).
  No new .Call entry is strictly required; predictBART already recomputes.

window: the only non-stan4bart change is the C1 dbarts engine fix (setState
  consistency at scale); it touches src/R_interface_bartcore.cpp / the bartcore
  setState validation, NOT the dbarts.h signatures, so the frozen-ABI /
  LinkingTo contract holds. Everything else is stan4bart-local R.

## Goal

Stop paying the always-on n x draws x chains memory cost of `bart_train` for
users who do not need per-draw BART fits materialized, by recomputing them from
the kept trees on access (dbarts predict over the training design matrix). Keep
the `fitted`/`extract`/`predict` surface intact and its returned values equal to
the stored path within a tight numeric tolerance. Decide the default policy from
the measured trade, not assumption: small-n users must not pay recompute latency
for a memory saving that does not exist at their scale, and large-n users (the
ones hurting) must be able to get the saving without a silent value change.

## Verdict up front

Recompute is a large, clean win for large n and a net loss for small n, with a
crossover near n ~= 550. The memory saving reaches 94-98% of the object at
n >= 10k. BUT two blockers sit in front of it, both PRE-EXISTING and both
independent of whether recompute becomes the default:

1. SERIALIZATION (the hardest constraint). The fit carries only a LIVE external
   pointer (`$sampler.bart`); the serializable tree state (`state.bart`) is
   built, consumed to make the pointer, then DISCARDED. After saveRDS/readRDS the
   pointer is a dead externalptr and `predict`/`extract(type="trees")` already
   error today ("predictBART called on NULL external pointer"). A recompute path
   MUST retain `state.bart` and lazily rebuild the pointer, or a reloaded fit can
   produce nothing.

2. setState-AT-SCALE. Building the stored-tree sampler
   (createStoredBARTSampler -> dbarts setState) REJECTS large fits today:
   reproduced at n=50000 x 200 draws ("state is not consistent with this
   sampler"); n=45000 x 200 and n=50000 x 40 both pass. keepTrees is therefore
   already broken at exactly the scale where recompute pays off. This is a dbarts
   engine bug, fixable without a dbarts.h change.

Exactness is NOT bitwise (~1e-13), and that is fine - the divergence is in the
STORED value, not in tree replay (below). Recommendation: land opt-in first
(store trees INSTEAD of fits), keep store-fits the default, revisit a
size-adaptive default once the two blockers are fixed and the contracts are
lived-in.

## Context read in code (file:line)

Fit assembly and the discarded state (the serialization gap):
  - stan4bart_fit.R:79-80: `if (control.bart@keepTrees) results$state.bart <-
    .Call(C_stan4bart_exportBARTState, sampler)` - per chain, serializable.
  - stan4bart_fit.R:645-656: combine step assembles `all_state` and calls
    `C_stan4bart_createStoredBARTSampler(control.bart, data.bart, model.bart,
    all_state)`, storing the resulting EXTERNAL POINTER as
    `attr(chainResults, "sampler.bart")`. `state.bart` is NOT propagated past
    this point.
  - stan4bart.R:456-457: `result$sampler.bart <- attr(chain_results,
    "sampler.bart")` - the fit keeps ONLY the pointer. package_samples
    (stan4bart.R:273-357) copies bart_train/test/varcount but never state.bart.
    => the fit object cannot survive a save/load for any tree-based op.

The recompute mechanism already exists:
  - init.cpp:332-383 predictBART: takes the StoredBARTSampler pointer + a test X,
    calls dbarts predict, returns n_test x numSavedSamples x numChains on the
    ORIGINAL response scale ("the restored state carries the fit's transform").
    Used for test data at generics.R:685.
  - Feeding it the TRAINING X reproduces `bart_train`. Training X is
    reconstructible from retained data: getTestDataFrames ->
    dbarts::makeTestModelMatrix(object$bartData, mf.bart) (test_data.R:135); the
    fit already retains `$frame` and `$bartData`, so recompute needs NO extra
    stored design matrix.
  - init.cpp:88, 351, 430, 489-492 numSavedSamples; init.cpp:478-542 getTrees;
    init.cpp:71-72,231,608-610 setTreeStorage (trees stored only for the
    non-warmup run: `setTreeStorage(sampler, isWarmup?0:1, isWarmup?0:numIter)`).

How bart_train is produced (source of the ~1e-13 exactness gap):
  - init.cpp:688-694: `run(...,&bartSamples->current)` then `trainingSample[j]
    -= sampler.bartOffset[j]`. The stored value is fl((offset + treeSum) -
    offset); predictBART returns treeSum directly. The two differ by one
    subtractive round-trip, not by tree-replay error.

Every reader of `$bart_train` (the compat surface):
  - generics.R extract (296-446): type in {ev, ppd, indiv.bart} reads
    `get_samples(object$bart_train, ...)` (424) or `$bart_test` (426); dims read
    at 315-319. Sizing at 656-661 (predict) and 656.
  - generics.R fitted (492-524): calls extract() then averages -> materializes
    the full n x draws x chains array to return a length-n mean.
  - generics.R predict (630-737): for a NEW newdata it ALREADY recomputes via
    predictBART (685) and requires keepTrees (681-682); with no newdata it falls
    back to extract() of the stored block (646-647).
  - stan4bart.R:315-357 packaging; tests/testthat/test-05-rng.R:26,44,60,62
    (expect_equal on `$bart_train`, stored-vs-stored reproducibility);
    test-11-callback.R:91-98 (expects `$bart_train` NULL when keepFits off).
  - EXTERNAL: bartCause consumes stan4bart fits; it should read through
    extract()/fitted(), not `$bart_train` directly - to be confirmed (Q3/e).

The consistency rejection (C1 blocker) is engine-side:
  - dbarts src/R_interface_bartcore.cpp:3944-3952: `restored = ... &&
    sampler.setState(state, currentPredictors); if (!restored) Rf_error("state is
    not consistent with this sampler")`. The false comes from the bartcore
    setState flat-tree validation, not from stan4bart.

Queued framing: TODO "bart-train-recompute" ("keepTrees semantics, a lazy seam
  behind fitted/extract, and an eager-storage option for fast-repeated-access
  users; no dbarts.h change"); NEWS.md:47-48 recorded it as left out of scope at
  the warmup-off landing.

## Measurements (fixed seed=99 friedman causal+ranef; installed stan4bart
0.0.14 / dbarts 1.0.0; n.trees=75 default; warmup NOT stored by default)

All projected to 4 chains x 1000 post-warmup draws. Storage is exactly linear
in draws x chains (bart_train object.size matched 8*n*S*C to 3 digits at every
scale); recompute latency is linear in n x total-draws (3.68e-7 s per obs-draw,
stable across scales); tree state is linear in stored draws and only WEAKLY
n-dependent (deeper trees at larger n).

1. SIZE and CROSSOVER. Kept-tree state measured directly via the dbarts sampler
   `$state` (== exportBARTState; the true serialized footprint):
   per-draw-per-chain ~= 4.15 KB (n=1k), 5.1 KB (n=10k), 6.65 KB (n=50k).

     scale        | bart_train stored | kept-tree state | saving | %
     -------------|-------------------|-----------------|--------|----
     n=1k  4x1000 |     32.0 MB       |    16.6 MB      |  15 MB | 48%
     n=10k 4x1000 |    320.0 MB       |    20.4 MB      | 300 MB | 94%
     n=50k 4x1000 |   1600.0 MB       |    26.6 MB      |1573 MB | 98%

   Tree state is ~n-independent (0.83 -> 1.02 -> 1.33 MB at 200 draws/1 chain for
   n=1k/10k/50k); bart_train is strictly linear in n. CROSSOVER (trees >= fits):
   8*n ~= per-draw-per-chain bytes ~= 4200-4400 -> n ~= 520-550. Below n ~= 550,
   recompute SAVES NOTHING (trees cost as much or more than the fits) AND adds
   latency; above it the saving grows with n toward 98%.

2. RECOMPUTE LATENCY (full n x draws x chains via predictBART; measured 0.13 s at
   n=1k/500 draws and 3.68 s at n=10k/1000 draws), projected to 4x1000 = 4000
   total draws:
     n=1k  ~1.0 s ; n=10k ~14.7 s ; n=50k ~74 s (over a minute).
   - extract(type in {ev,ppd,indiv.bart}) = the full matrix: unavoidable, and at
     n=50k it is ~74 s per call, plainly. This is the cost users trade the memory
     for.
   - fitted() = posterior mean: today reads stored then averages (0.08 s at
     n=10k). It CAN stream without holding the full matrix: recompute in
     row-blocks and accumulate the mean (measured blocked-stream matches the
     stored mean to 5e-14 and runs at ~the full-recompute latency but with peak
     memory = block_n x draws x chains, not n x draws x chains). So the memory
     saving is preserved for fitted(); the cost is latency (0.08 s -> ~15 s at
     n=10k x 4000 draws). No C change needed for streaming - predictBART already
     accepts an arbitrary test-row block; sample-subsetting is NOT available
     (predictBART uses all saved samples), so streaming is over ROWS only.

3. keepTrees SAMPLING OVERHEAD. Wall-time OFF vs ON at n=10k: 3.38 s vs 3.18 s,
   i.e. negligible (within run-to-run noise). Extra memory held during/after
   sampling = the kept-tree state itself (row "kept-tree state" above; ~20 MB at
   n=10k x 4000 draws, ~n-independent). Storing trees is essentially free in time
   and cheap in memory; the cost is entirely on the recompute side.

4. EXACTNESS. predictBART over the training X vs stored `bart_train`: max abs
   diff 7.2e-14 (n=1k), 1.1e-13 (n=10k). NOT bitwise identical(). The divergence
   is the stored value's offset round-trip (fl((offset+treeSum)-offset) vs
   treeSum, init.cpp:688-694), so recompute is arguably the CLEANER quantity;
   tree replay itself is exact. Verdict: the exactness gate must be a tight
   tolerance (max abs diff < 1e-10, or relative < 1e-12), never identical().
   Reproducibility (test-05-rng) is unaffected: same seed -> same trees -> same
   recompute; and stored-vs-stored comparisons there do not change if the stored
   path is retained.

5. SERIALIZATION (round trip). After saveRDS/readRDS: `$sampler.bart` is a
   non-null externalptr whose address is dead; `predict(fit2, ...)` errors
   "predictBART called on NULL external pointer"; `state.bart` is not present in
   the fit. Recompute-after-reload is IMPOSSIBLE today. This is the design's
   hardest constraint and the reason C2 is a prerequisite, not a nicety.

6. setState-AT-SCALE (blocker). createStoredBARTSampler fails for large fits:
   n=50000 x 200 draws -> "state is not consistent with this sampler";
   n=45000 x 200 draws OK; n=50000 x 40 draws OK; n=40000 x 350-iter OK. The
   failure is data/scale dependent (dbarts setState flat-tree validation,
   R_interface_bartcore.cpp:3944-3952) and already breaks keepTrees in the
   large-n regime recompute targets. Must be root-caused before the arc delivers
   at scale.

## Options (with costs)

(A) Keep store-fits as today; do nothing. Cost: the 94-98% memory cost stays for
    every large-n user. Baseline.

(B) OPT-IN "store trees instead of fits" (recommended first step). A new arg
    (e.g. bart_args store = c("fits","trees")) that turns keepTrees on and DROPS
    the accumulation of `bart_train`/`bart_test` into the fit; fitted/extract/
    predict route through the lazy recompute seam (C3). Saves the full block at
    the cost of recompute latency on access (M2). Predictable API, no silent
    value change, `$bart_train` simply absent under the opt-in. Cost: "an option
    nobody flips helps nobody" - large-n users must know to set it.
    Impl: M, R-side; depends on C1+C2.

(C) SIZE-ADAPTIVE default: recompute when n (or n x draws) exceeds a threshold
    (~a few thousand, where saving > ~90% and full-matrix latency < ~15 s),
    store-fits below it. Helps the hurting users automatically. Cost: fitted()
    latency and object size become data-dependent (surprising), and returned
    bart-component values shift by ~1e-13 at the threshold vs a store-fits fit.
    Impl: +S over (B) (a size check at assembly). Needs C1+C2 solid first.

(D) LAZY-ALWAYS + mean cache: never store bart_train; recompute on every access;
    memoize only the posterior-mean vector (n doubles, negligible) for repeated
    fitted() calls. Simplest single contract. Cost: small-n users pay recompute
    latency for a saving that does not exist below n ~= 550. Impl: M.

Caching policy (question c) for (B)-(D): do NOT cache the full recomputed matrix
    - that reintroduces exactly the memory we removed. Recompute-every-call for
    extract(full matrix); streamed row-blocked recompute for fitted(mean);
    optionally memoize the length-n mean vector only. Per-op, not global.

## Commit partition (each commit installable + testable; the exactness gate
lands BEFORE any behavior change)

C0 - Gate infrastructure (no behavior change). Add the exactness + latency
  harness: fixtures across tiers (n small/med, gaussian + binary) that fit with
  keepTrees, compare predictBART-over-training-X against stored `bart_train` to
  the tolerance in M4, and record baseline sizes + latencies. Sizes: S. Gate:
  harness green on the CURRENT build (proves the tolerance is real before any
  code depends on it). Abort: if the divergence is ever > 1e-10, stop and
  re-derive the source before proceeding.

C1 - dbarts engine: setState-at-scale fix (PREREQUISITE, dbarts-side, no .h
  change). Root-cause the n=50000 x 200-draw rejection at
  R_interface_bartcore.cpp:3944 / bartcore setState (candidate: an int-width or
  data-dependent flat-tree validation edge). Gate: n=50000 x 200-draw keepTrees
  fit assembles and round-trips predict; add a tests/cpp state round-trip at
  scale. Size: M (investigation-heavy). Abort/branch: if it cannot be fixed
  without a dbarts.h ABI change, escalate (Q2) - the arc's large-n value is
  gated on this.

C2 - Serialization retention (stan4bart, prerequisite; ALSO fixes a standing
  bug). Retain `state.bart` in the fit object (package_samples / stan4bart.R
  assembly) and make `$sampler.bart` lazy: an accessor that returns the live
  pointer or, when it is null/dead (fresh reload), rebuilds it from `state.bart`
  via createStoredBARTSampler. Sizes: the retained state is the "kept-tree state"
  row (16-27 MB at scale), replacing nothing yet. Gate: saveRDS/readRDS then
  predict/extract(type="trees") equal the in-session result to tolerance (today
  they error); no change to any store-fits value. Size: S-M, R + existing C
  entries. Depends on C1 for large n.

C3 - Recompute seam behind readers + the opt-in (Option B). Route
  fitted()/extract(type in {ev,ppd,indiv.bart})/predict(no newdata) through the
  lazy recompute when the fit was built trees-only; stream the mean in row-blocks
  (M2). Add the store=c("fits","trees") arg (default "fits"). Gate: a trees-only
  fit's fitted/extract equal a store-fits fit's to tolerance across tiers;
  record extract(full matrix) latency and the memory saving; `$bart_train` absent
  under the opt-in with extract() documented as the supported surface. Keep
  store-fits the DEFAULT (no silent flip). Size: M, R-side. Update
  test-11-callback expectations if the opt-in shares the keepFits path.

C4 - Default-policy decision (Q1). With C0-C3 contracts in hand, the author picks
  among: keep opt-in (B), size-adaptive default (C), or lazy-always+cache (D).
  Gate: the recorded exactness tolerance + latency table; a DOCUMENTED threshold
  if adaptive. Size: S (policy + docs + NEWS). No code lands here beyond wiring
  the chosen default.

## Landings through C2 (2026-07-16)

C1 = dbarts 5b98ea5. The setState-at-scale failure was not a restore bug:
the empty-leaf veto's finite -1e7 sentinel was out-penalized by legitimate
branch scores at scale, empty leaves entered live trees, and restore
correctly rejected the degenerate state. The sentinel is now -HUGE_VAL -
draw-neutral at every gated scale (both sentinels underflow exp to zero;
all three bitwise anchors identical) - with a fast regression test and the
empty-leaf-veto design note corrected.

C0 + C2 = d6777d6. The exactness harness: recompute-from-trees matches
stored bart_train to 6.4e-14 worst-case across gaussian/binary tiers
including past the old restore threshold (gate 1e-10; the stored value
carries the offset round-trip and is the noisier quantity). C2 retains
state.bart on keepTrees fits and rebuilds the sampler pointer lazily
(liveness-probed, cached per session, in-session path untouched,
keepTrees = FALSE fits byte-identical), fixing the standing
predict-after-reload failure; test-14-serialization proves the round trip
in a fresh R process. Suite 275/0; continuous_nc1 tier PASS (sampling
untouched). The C0+C2 implementer was lost to a usage limit mid-task; the
orchestrator finished its remainder (one test-harness fix: deparse of a
long .libPaths wraps - collapse before scripting).

### C3 landing and C4 resolution (2026-07-16) - ARC CLOSED

C3 = 2bd83a1. R-side seam only (the C-side keepFits gate bundles too much;
the accumulation point in package_samples is the smaller correct cut).
store = "trees" forces keepTrees, drops bart_train/bart_test, and routes
extract/fitted/predict(no newdata) through recompute - means stream in
~16 MB row blocks, full matrices materialize with the latency documented
in the Rd. Equality gates: parametric draws identical between modes, bart
surfaces within 1e-10 (measured 1e-12 to 1e-14), gaussian + binary,
in-session + reloaded. Suite 307/0; harness ALL PASS; nc1 tier PASS;
check Status OK. Known fallbacks: fitted streaming defers to full
materialization for component-replacing offsets, na.exclude gaps, or a
missing parametric design (guarded in code).

C4 RESOLVED to Option B held: store = "fits" stays the default, the
opt-in is the shipped policy, and a size-adaptive default is revisited
only after the contracts are lived-in - the conservative branch changes
nothing for existing users, honors the measured latency of full-matrix
recompute at scale (~74 s at n=50k), and leaves the 94-98% memory win one
argument away for the users who need it.

## Verification

- Exactness (C0, re-run every commit): max abs diff stored-vs-recompute < 1e-10
  across gaussian + binary, small + medium n; NEVER identical()-based.
- Serialization (C2): saveRDS/readRDS round trip -> predict/extract(trees) equal
  in-session to tolerance; the standing "NULL external pointer" error is gone.
- Scale (C1): n=50000 x 200-draw keepTrees assembles + predicts; dbarts tests/cpp
  state round-trip at scale passes.
- Reproducibility: test-05-rng still green (same seed -> same draws); if a value
  path moves to recompute, refresh only the snapshots that were bitwise, since
  they shift by ~1e-13.
- Latency + memory: record extract(full matrix) seconds and peak-RSS delta
  (store-fits vs trees-only) at each tier; confirm the projected saving.
- Compat: fitted/extract/predict unchanged for store-fits fits; under the opt-in,
  a documented, informative error/absence for `$bart_train` direct reads.

## Open questions (for the author)

Q1 (default policy - the genuine fork). Opt-in (B) vs size-adaptive default (C)
  vs lazy-always+cache (D). Measured trade: (B) is API-predictable and never
  silently changes a returned value, but needs the large-n user to know the flag;
  (C) auto-helps the hurting users (recompute above ~n=2-5k where saving > 90%
  and a full extract is < ~15 s) at the price of data-dependent fitted() latency
  and a ~1e-13 value shift at the threshold; (D) is one simple contract but
  charges small-n users (< n ~= 550) latency for zero saving. Recommendation:
  ship (B) at C3, revisit (C) at C4 after the contracts are lived-in.

Q2 (the C1 constraint). The setState-at-scale rejection is a dbarts ENGINE bug
  (not a dbarts.h ABI change). Is fixing it in dbarts acceptable within this
  arc's "no dbarts.h change" framing, or must recompute be scoped to n below the
  failure threshold until dbarts is fixed on its own track? The arc's large-n
  payoff depends on the answer.

Q3 (exactness + downstream contract). Accept the ~1e-13 non-bitwise divergence in
  returned BART-component values (recompute is the cleaner quantity), and update
  any identical()-based expectation? Confirm no consumer - notably bartCause -
  reads `$bart_train` directly rather than through extract()/fitted(), since the
  opt-in makes that field absent.

Q4 (caching). Is a memoized posterior-mean vector (n doubles) worth holding for
  repeated fitted() calls, or keep the recompute fully stateless?
