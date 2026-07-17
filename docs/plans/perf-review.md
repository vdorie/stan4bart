# perf-review

scope: RUNTIME performance review (compile time explicitly out of scope, VD
  2026-07-15 TODO). Profile representative fits and rank optimization
  candidates BEYOND the BART core. The dbarts engine (its own repo, frozen
  dbarts.h) is not this package's to optimize; where a bottleneck lands inside
  dbarts it is recorded as a routed finding, not dug into.

budget: MEMO. Measures where per-iteration and post-processing wall-clock go
  across continuous/binary x with/without random effects x n = 1e3..1e5
  single-chain, then prices candidates. No code changed.

## Verdict up front

For the large-n-with-random-effects case (single-chain large-n is the stated
user priority), per-iteration wall splits between two co-owners whose ratio is
set by WALNUTS tuning quality:

  - dbarts BART draw (init.cpp:688): 54-72% of per-iter compute.
  - WALNUTS parametric step (init.cpp:624): 27-45%, of which the hand-written
    log-density+adjoint `ParametricModel::eval` is ~24-40% ALONE.

The single biggest runtime lever is the WALNUTS leapfrog-step count per
transition, which is set by the frozen step-size/mass tuning and is extremely
tuning-sensitive: at n=1e5 the SAME cont/ri model samples at 69 ms/iter after a
20-iter warmup but 28 ms/iter after a >=100-iter warmup (2.2x, measured). Under
adequate warmup BART dominates (~72%); under short warmup the parametric step
does (~65%). The parametric target `eval` is the largest stan4bart-owned kernel
either way. Everything else (R<->C marshalling, per-iteration copies, results
assembly, most post-processing) is second-order at realistic draw counts, with
two cheap exceptions worth taking (fitted()'s apply->rowMeans; the ppd rnorm).

The 54-72% that is BART is a dbarts-side finding, routed, not actionable here.

## Method and measurement quality

Installed stan4bart 0.0.14 / dbarts 1.0.0 (both -g -O2, symbols intact), arm64
macOS, quiet host, single chain cores=1, n.trees=75 unless noted. Friedman
signal reusing the reference-tier simulators. Three instruments:
  - per-iteration SAMPLING wall by differencing two `iter` at matched warmup
    (cancels setup+warmup); 3 reps where noise mattered.
  - native `sample` (1 ms, ~13-22k stacks) on the running fit for the C++
    BART-vs-parametric split (top-of-stack leaf self-time).
  - Rprof (memory profiling on) for the R surface: whole-fit setup/assembly and
    the extract/fitted/predict consumers.

Quality caveats: (1) per-iter differencing at n=1e3 is noise-dominated (first-fit
JIT); trust n>=1e4. (2) `sample` is statistical, leaf attribution +-~2%. (3) the
parametric per-iter is tuning-dependent (see Verdict); two profiler regimes
(warmup 20 and 150) bracket it, production warmup=iter/2 sits near the
well-warmed end. (4) single-host arm64/NEON; ratios move with n.trees (more
trees -> larger BART share), q (RE levels), and family.

## Where wall-clock goes

### Per-iteration SAMPLING wall (ms), by workload x n

    workload      n=1e3   n=1e4   n=1e5
    cont/none     ~0*     1.23    19.2      (BART only, K=1 fixed effect)
    cont/ri       0.22    2.98    28-31**   (+ WALNUTS over q RE levels)
    cont/nc2      0.30    3.02    33.5
    binary/ri     0.31    3.67    52.1***

    *  n=1e3 cont/none is noise (first-fit); ** well-warmed (warmup>=100);
       65 ms at warmup=20. *** binary adds per-iter probit-latent draw +
       setResponse (init.cpp:708-712), measured at warmup=20.

Scales ~linearly in n (BART slightly super-linear, consistent with the engine's
O(n log n) partition/sort). Adding random effects roughly doubles-to-triples
per-iter at large n; without RE the parametric step is negligible (few
dimensions, few leapfrog steps).

### C++ split (native `sample`, n=1e5 cont/ri, q=1000 RE levels, n.trees=75)

Well-warmed (warmup=150), sampling phase, ~12.7k compute stacks:

    component                                     share   owner
    ------------------------------------------    -----   ------
    BART draw  dbarts_sampler_run (init.cpp:688)   72%    dbarts
      misc_computeIndexedSufficientStatistics       21%   dbarts
      misc_setIndexedVectorToConstant              18.5%  dbarts
      misc_partitionIndices/Range_neon (SIMD)       20%   dbarts
      bartcore::Chain::run + refreshSubtree         12%   dbarts
    WALNUTS step  paramSampler->run (init.cpp:624)  27%   stan4bart
      ParametricModel::eval (logp + adjoint)        24%   stan4bart
      Eigen X*beta / X^T g matvec                    2%   stan4bart (K=1: tiny)
      parametricMean + re-emit eval                 ~1%   stan4bart
    R glue / GC (xcopyReal, RunGenCollect)          ~2%   R

Under-warmed (warmup=20) the same fit reads BART 54% / WALNUTS 45%, eval 40%:
short warmup leaves the step size mis-tuned, WALNUTS takes many more leapfrog
steps, and each step is one O(n) `eval`. eval self-time is the top single leaf
in both regimes.

### q- and warmup-sensitivity (n=1e5 cont/ri, sampling per-iter, ms)

    RE levels q:  none(0)  ri(20)  ri(200)  ri(1000)  nc2(2000)
    per-iter:      22.5     55.1    34.4      65.4      33.5      (warmup=20)

    warmup:        20      100      400
    per-iter:      69.5    31.6     28.1                          (ri, q=1000)

Non-monotone in q but monotone-improving in warmup: the driver is leapfrog
COUNT (set by frozen tuning/geometry), not q or n per se. `eval` cost per call
is ~O(n) and roughly fixed across these; the count is what swings 2-3x.

### Setup (marshalling) and results assembly

Rprof of a whole n=1e5 cont/ri fit (iter=140): the `.Call` sampling is 84.6%;
the R-side setup+assembly is ~15% at this short iter, dominated by mkReTrms /
KhatriRao (sparse RE design), factor(), and package_samples' array(sapply(...))
build of the n x draws block (277 MB transient here; stan4bart.R:343-346). BOTH
are FIXED per-fit costs, so their share -> 0 at realistic draw counts (iter in
the thousands). R<->C marshalling and the per-iteration C-side copies
(setOffset/setResponse/getParametricMean memcpys, ~5 x n doubles/iter;
_platform_memmove ~3% in the profile) are not a wall-clock bottleneck.

### R post-processing (n=1e4, 1000 draws, 2 chains; object 109 MB)

    extract(ev)          0.22 s      fitted(ev)          0.18 s
    extract(ppd)         0.40 s      fitted(ppd)         0.53 s
    extract(indiv.bart)  0.002 s     as.matrix           0.001 s
    extract(ranef)       0.001 s     extract(ev, test)   0.025 s

Rprof of a consumer sweep (by self): rnorm 42% (ppd noise draws,
generics.R:549-564), t.default 20% + aperm 6% + array 5% (array permute/alloc),
apply 8% (fitted's average_samples_f apply(x,MARGIN,mean), generics.R:606-616),
%*% 2.4% (X%*%beta in fitted_fixed). ~7 GB transient allocated over the 4x sweep.
Post-proc scales ~linearly in n x draws: extract(ev) at n=1e5/4000 draws ~9 s,
fitted(ppd) ~20 s+. Parametric-only surfaces (fixef/ranef/sigma/Sigma/as.matrix)
are ~instant - they read the small stan block, not the n x draws BART block.

## Ranked opportunities - stan4bart-side

1. WALNUTS leapfrog-count / warmup adequacy (HIGHEST VALUE).
   Payoff: the parametric step is 27-65% of per-iter and is dominated by leapfrog
   count; adequate warmup already halves sampling wall (69->28 ms, 2.2x). The
   lever is making the frozen tuning good with less warmup, or steering large-n-RE
   users to enough warmup. Candidates: a better initial step size / mass seed
   (walnuts_sampler.cpp:277 starts unit mass, step 0.1), a warmup schedule sized
   to the model dimension rather than the outer sweep count, or a documented
   large-n-RE warmup floor.
   Cost: S (defaults/docs) to M (WarmupConfig tuning; partly upstream WALNUTS
   behavior). Risk: MEDIUM - fewer steps can hurt mixing; must gate on the
   distributional-equivalence tiers (walnuts.md) before shipping any tuning change.

2. ParametricModel::eval - scratch reuse + loop fusion.
   Payoff: eval is 24-40% of per-iter compute and runs once per leapfrog step. It
   allocates eta(N), g_eta(N), b(q), g_b(q), g_z_b(q) and per-block temporaries on
   EVERY call (thousands/sec), and makes ~6-8 separate O(N) passes
   (parametric_model.hpp:400-424). Preallocating these as reused members and
   fusing the forward eta+likelihood and the backprop passes could shave ~10-25%
   of eval -> ~3-9% of per-iter at large n.
   Cost: M. Risk: MEDIUM - eval is the FD-gated hot path; reuse buffers WITHOUT
   reordering arithmetic to keep draws bit-identical, and re-run the FD gate
   (test-12) + distributional tiers. Buffer reuse alone is numerically inert.
   LANDED (the buffer-reuse half; measured outcome differed from the estimate):
   naive member-hoisting alone REGRESSED eval ~5% - scratch reached through
   `this` defeats the local-no-alias codegen in the hand-written O(N) loops -
   so the landed form pairs the reuse with __restrict pointers on the three
   hot loops (a true guarantee, the member buffers are disjoint). Net eval win
   ~3-7% (interleaved microbench, arm64 clang), whole-fit per-iter ~3% at
   n=1e4 within noise. Draws bit-identical (identical() TRUE, continuous and
   binary with RE), FD gate 10/10 at the landed baseline. Loop FUSION remains
   unlanded: it reorders arithmetic and would need the distributional tiers -
   the measured ceiling above suggests it is not worth that gate.

3. fitted() apply -> rowMeans (CHEAPEST WIN).
   Payoff: average_samples_f uses apply(x, MARGIN, mean) (generics.R:606-616);
   apply is ~8% of the post-proc sweep and ~2-3x slower than rowMeans/.colMeans
   on the combined 2D matrix. fitted(ev) ~0.18 s -> ~0.07 s at n=1e4, linear to
   larger. Cost: S (a few lines). Risk: VERY LOW.

4. fitted(ppd) - skip the noise for the mean (continuous).
   Payoff: fitted(ppd) draws n x draws x chains normals then averages them to ~0;
   E[ppd]=E[ev], so the mean can short-circuit to the ev path (0.53 s -> ~0.18 s,
   skipping the 42% rnorm). Cost: S. Risk: LOW but SEMANTIC - only the ppd MEAN;
   extract(ppd) intervals still need the draws. Optional / discuss.
   RESOLVED 2026-07-16 (VD): NO-GO on the short-circuit. fitted(ppd) is
   implemented as extract(ppd)-then-average on the same RNG stream, so under a
   fixed seed it reproduces an average of extract(ppd) bitwise - a load-bearing
   reproducibility contract the short-circuit would break irreparably. Landed
   instead: a once-per-session fitted(ppd) message pointing mean-only users at
   type = "ev" (skipped for weighted binomial, where the ppd mean is the
   weights times the ev mean), plus the docs note in man/generics.Rd. bartCause
   was checked and does not call fitted(ppd) programmatically.

5. extract(ev/ppd) - reduce transient arrays.
   Payoff: extract materializes three full n x draws x chains arrays
   (indiv.bart/fixef/ranef) then sums (generics.R:471-539); ~7 GB transient over
   the measured sweep. Accumulate in place / trim aperm+array copies to cut wall
   and peak RSS at large n. Cost: M. Risk: LOW-MEDIUM (array-layout care).

6. package_samples assembly peak (minor).
   array(sapply(...)) double-copies the bart_train block at assembly
   (stan4bart.R:343-346); peaks ~2x the block at n=1e5 x 4000 x 4 chains. Cost:
   S-M. Risk: LOW. Small wall; matters only for peak RSS at extreme scale, and is
   already sidestepped by the existing store="trees" opt-in.

## dbarts-side findings (routed, NOT actionable in this package)

- BART is 54-72% of per-iteration compute at n=1e5; the hot leaves are
  misc_computeIndexedSufficientStatisticsFast (~21%), misc_setIndexedVectorToConstant
  (~18.5%), and the misc_partition{Indices,Range}_neon SIMD (~20%), all inside
  bartcore change-move / refreshSubtree. Per-sweep BART cost scales ~O(n log n)
  in n and linearly in n.trees; any large-n speedup is dominated by these. Route
  to the dbarts engine. (setIndexedVectorToConstant at ~18.5% is a broadcast/fill
  that may be reducible - flag for dbarts.)

## Incidental (correctness, not perf)

- predict(fit, newdata, "ev") on a fit WITH random effects throws
  "non-numeric argument to binary operator" in the indiv.bart+fixef+ranef+offset
  sum (generics.R predict newdata path, offset handling). Reproduced at
  n=1e4/keepTrees=TRUE. Extract on stored test data (sample="test") is fine.
  Robustness bug in the stan4bart R surface, unrelated to this review's scope.

## Recommendation

Single highest-value item: attack the WALNUTS leapfrog-count / warmup-tuning
lever (#1) - it is where the parametric wall actually lives and a 2.2x is on the
table from tuning alone, but it is a sampler-quality change and must clear the
distributional gate. Pair it with the two cheap, low-risk wins that need no
gate: fitted()'s apply->rowMeans (#3) and eval scratch-reuse (#2, buffer reuse
half only, which is numerically inert). Treat the BART 54-72% as a dbarts matter.
Everything in setup/marshalling/assembly is fixed-cost and amortizes; do not
spend there.
