# callback-estimators

scope: DESIGN MEMO, no implementation. Determines whether/how a per-iteration
  callback can let bartCause store draws x estimator arrays instead of the
  n x draws x 2 counterfactual fit surfaces, grounded in what bartCause actually
  computes. Covers the estimand taxonomy (which outputs genuinely need the full
  individual matrix post-fit), the callback contract (what stan4bart must hand a
  callback), the storage mode (callback results as the only individual-level
  record, composing with store="trees"), the bartCause adoption surface, and the
  measured overhead/exactness. The BART engine and dbarts.h are FROZEN; nothing
  here needs a dbarts change. The callback seam and keep_fits=FALSE ALREADY EXIST
  in this worktree (init.cpp / stan4bart.R / generics.R); the work is a contract
  refinement plus a reducing callback written in bartCause.

budget: memo only. If pursued: stan4bart side is S (a documented contract + one
  optional convenience so the callback need not re-derive ev in R); bartCause
  side is M (a reducing-callback constructor + routing the average-estimand
  summary paths to extract(fit,"callback") when present). No C++ sampler change.

## Verdict up front

It works, for the estimands that are the common case, and the mechanism is
mostly already here. Every bartCause AVERAGE-estimand summary (print, summary
target=pate/sate/cate incl. quant/hpd intervals, fitted(pate/sate/cate),
group.effects) reduces to per-draw scalars (draws x 1) or per-draw group scalars
(draws x G): none of them needs the n x draws matrix to survive the fit. The
ONLY user-facing outputs that genuinely need the individual surfaces post-fit are
the ones that RETURN individual quantities - extract/fitted/predict of
icate/ite/mu.*/y.*/p.score/p.weights - plus refit() (which recomputes averages
under a new common-support subset from the stored icate). plot_indiv/plot_support
need only per-individual posterior MEANS (n x 1), a reduction, and the already
stored sd.obs/sd.cf.

The existing callback already sees enough for the Gaussian CATE family: test-11-
callback.R reconstructs the exact per-draw ev for train and test from
(yhat_train, yhat_test, stan_pars) plus the fit data closed over in callbackEnv,
verified bit-equal to extract(fit,"ev"), and runs with keep_fits=FALSE storing
NOTHING but the callback array. Measured: a reducing callback that emits per-draw
ate/att/atc is bitwise-identical to the post-hoc reduction (max abs diff 0) and
adds +1.1% wall time. The gaps are narrow and NAMED below: binary needs the
response-scale (pnorm) mean, not BART-only; SATE/PATE add PPD noise (bitwise
divergence, statistically fine); common-support cutting is a post-hoc
data-dependent subset that a single-pass per-draw average cannot honor; and
TMLE/p.weight cannot move to the response callback at all (they need the
treatment model's per-draw p.score, which lives in a separate, unaligned fit).

Recommendation: land a documented "reducing callback" contract + store="trees"
composition; bartCause routes the default (method.rsp="bart", commonSup="none")
average estimands through it and keeps surfaces only when the user asks for
individual draws or common-support cutting. See Recommendation.

## Goal

Answer, grounded in bartCause's code: (a) for each estimand/summary, does it need
the full ITE matrix, per-draw averages (draws x 1), per-draw group averages
(draws x G), or a running summary - and which USER-FACING outputs require
individual draws to exist post-fit; (b) what per-iteration quantities a callback
needs to compute obs-minus-counterfactual averages, and whether the existing
callback sees them; (c) a storage mode where the callback's per-draw results are
the only individual-level record, composing with store="trees"; (d) what changes
in bartCause and what cannot switch; (e) overhead, warmup/multi-chain behavior,
and the exactness relationship to the C0 harness.

## Context read in code (file:line)

How bartCause drives stan4bart (grep: bartCause calls stan4bart::stan4bart, never
a callback today):
  - Response model: responseFit.R:104,132-136 route to getStan4BartResponseFit
    when `parametric` is supplied; responseFit.R:56 defaults chains=10;
    responseFit.R:49 passes `treatment = <name>` to stan4bart (NOT an explicit
    test set) - stan4bart's first-class `treatment` arg auto-builds the
    counterfactual (flipped-treatment) test design, which is why
    extract(sample="test") returns mu.hat.cf.
  - responseFit.R:65-79: mu.hat.obs = extract(fit, sample="train"), mu.hat.cf =
    extract(fit, sample="test"); both are n x draws x chains ev surfaces. These
    two surfaces are the n x draws x 2 the author cites.
  - responseFit.R:81-84: sd.obs/sd.cf = per-obs posterior SD of each surface (an
    n x 1 REDUCTION), fed to getCommonSupportSubset. Common support is thus
    already reduced to n-vectors at fit time and stored - it does not need the
    matrix post-fit; but computing it DOES consume the matrix during the fit.
  - Treatment model: treatmentFit.R:104 uses stan4bart when `parametric` given;
    treatmentFit.R:181-195 keeps `samples` = n x draws p.score (a THIRD surface),
    reduced to a posterior-mean p.score for covariate use; the full samples are
    retained only for extract(p.score) and TMLE/p.weight.

How the average estimands are computed (the reduction is already a per-draw
average):
  - generics.R:392-409 extract(pate/sate/cate): forms samples.indiv.diff =
    (surface_obs - surface_cf) * trtSign - pate uses PPD draws y.obs.ppd - y.cf,
    sate uses y.obs - y.cf, cate uses mu.hat.obs - mu.hat.cf - then
    getEstimateSamples -> averageDifferences.
  - generics.R:32-69 averageDifferences: subsets to att(trt>0)/atc(trt<=0)/ate AND
    commonSup.sub, then means over the observation dim -> a per-draw (x chain)
    scalar. So extract(fit,"cate") is ALREADY draws x 1; the full matrix is only
    an intermediate.
  - summary.R:228-244: the estimate functions request samples.icate =
    extract("icate","all") (FULL n x draws) OR samples.cate/sate = draws x 1.
  - summary.R:148-165 getCATEEstimate.bart: given samples.icate, immediately
    rowMeans -> draws x 1; est=mean, sd=sd; quant/hpd interval on the draws x 1.
    => reducible to a per-draw CATE scalar.
  - summary.R:64-86 getPATEEstimate.bart.var.exp: needs per-draw MEAN of icate AND
    per-draw VARIANCE of icate across people (apply(.,1,var)); both are per-draw
    scalars. => reducible to draws x 2.
  - summary.R:88-108 getPATEEstimate.bart.ppd: samples.cate=rowMeans (draws x 1),
    sd = sqrt(var(samples.cate) + 2*mean(sigma^2)/n) [gaussian] or a per-draw mean
    of samples.obs*(1-obs)+cf*(1-cf) [binary]. => reducible to draws x 1 + sigma +
    (binary) a per-draw scalar.
  - summary.R:118-146 getSATEEstimate.bart: samples.iscate=(y-mu.cf)*(2trt-1),
    rowMeans -> draws x 1, sd adds a per-draw noise term. => reducible.
  - summary.R:269-314 group.effects: the same reductions per group level ->
    draws x G, then recombined; still no n x draws survivor.
  - summary.R:50-62,110-116,140-146,159-165 the quant/hpd intervals operate on the
    draws x 1 sequence of per-draw averages - NOT the matrix.

User-facing surfaces that DO need individual draws post-fit:
  - generics.R:299-440 extract(icate/ite/mu.obs/mu.cf/mu.0/mu.1/y.cf/y.0/y.1/
    p.score/p.weights): returns n x draws individual arrays.
  - generics.R:250-297 fitted(...): averages the above over draws -> n x 1.
  - generics.R:87-248 predict.bartcFit: individual predictions on newdata;
    already REQUIRES keepTrees (generics.R:105-108,149-152) - i.e. it recomputes
    from trees, it does not read the stored surfaces.
  - generics.R:465-611 refit.bartcFit: recomputes est from extract("icate")
    (generics.R:503) under a new common-support rule - needs the full icate.
  - plot.R:70-83 plot_indiv -> fitted(type) = n x 1; plot.R:146 plot_support ->
    fitted("icate") = n x 1. Both need only the per-individual posterior mean, or
    the matrix long enough to reduce it.
  - responseFit.R:297-462 getPWeightEstimates/getTMLEEstimates: per-obs influence-
    function estimators over the FULL per-draw mu.hat.0/mu.hat.1/p.score vectors
    (+ iterative fluctuation for TMLE), producing a per-draw (est,se). These
    consume the matrix per draw AND the treatment model's per-draw p.score.

The existing callback seam (this worktree):
  - init.cpp:569-591: callbackClosure = Rf_lang4(callback, yhat_train, yhat_test,
    stan_pars); yhat_test only allocated if numTestObservations>0.
  - init.cpp:688-694: yhat_train = the BART-only training fit, offset SUBTRACTED
    (pure tree component); init.cpp:715-717 copies it and yhat_test (also the pure
    BART test component - test-11 verifies the reconstruction below).
  - init.cpp:718 copyOutParameters(REAL(stan_pars), keepFits ? -1 : 0);
    parametric_sampler.hpp:96 copies the full num_pars row, NAMED (init.cpp:581-
    587) so the callback selects ^beta|gamma / b. by name.
  - init.cpp:714-776: the callback fires every iter with NO isWarmup guard - it
    runs during warmup too; results are memcpy'd into a callbackResultLength x
    numIter block, shape/dimnames taken from the first return.
  - init.cpp:596,812-820: numStorageSamples = keepFits ? numIter : 1; with
    keep_fits=FALSE the run returns ONLY the callback list (no bart/stan stores).
  - stan4bart.R:227-229: callbackEnv = list2env(result, parent=baseenv()) - the
    callback's parent.frame() exposes the WHOLE fit: frame (with y and treatment),
    X, reTrms$Zt, test$X, test$reTrms$Zt, X_means, weights, offset.
  - stan4bart.R:418-437: package_samples stacks per-chain callback blocks into
    result$sample$callback (K x draws x chains) and result$warmup$callback.
  - generics.R:263-313 extract(type="callback"): returns object$callback,
    combine_chains-aware; warmup behind include_warmup. So callback results are a
    first-class extractable surface already.
  - generics.R:525-542 extract("ev"): ev = indiv.bart + indiv.fixef + indiv.ranef
    (+ default offset), and for binary result <- pnorm(...). This is exactly what
    a callback must reproduce to match mu.hat.obs/mu.hat.cf.
  - test-11-callback.R:12-73: a callback reconstructs test ev from yhat.test +
    test$X%*%fixef + crossprod(test$reTrms$Zt, ranef), asserted EQUAL to
    extract(fit,"indiv.*"/"ev", sample="test"); lines 76-99: with keep_fits=FALSE
    the same fit has bart_train/bart_test/bart_varcount/stan/warmup all NULL. This
    is a working proof of the proposed storage mode.

store="trees" recompute seam (C3-landed), for the occasional individual-draw ask:
  - generics.R:239-261 bart_needs_recompute / recompute_bart_block: replays kept
    trees through C_stan4bart_predictBART over the train/test design to reproduce
    the n x draws x chains block the fit did not store.
  - bart-train-recompute.md: recompute is ~1e-13 (not bitwise). Its two
    historical blockers are both FIXED and gated (see its Landings block):
    serialization now retains state.bart with a lazy pointer rebuild, and the
    setState-at-scale rejection was the dbarts empty-leaf-veto bug, fixed at
    dbarts 5b98ea5 with the harness proving past-threshold restores. So
    "recompute individual draws on demand" is robust at the scales where the
    surface memory hurts most.

## The estimand taxonomy

Quantity classes: FULL = n x draws matrix must exist post-fit; d1 = per-draw
scalar (draws x chains); dG = per-draw group scalars (draws x G x chains);
run = a running per-observation moment reducible to an n-vector (sd.obs/sd.cf);
n1 = per-individual reduction over draws (posterior mean, n x 1).

| bartCause output (path) | quantity the summary needs | class | needs FULL post-fit? |
|---|---|---|---|
| print.bartcFit (summary.R:433-448 -> fitted cate/pate) | per-draw average | d1 | NO |
| summary target=cate (summary.R:148-165) | per-draw CATE mean | d1 | NO |
| summary target=pate var.exp (summary.R:64-86) | per-draw mean + within-draw var | d1 x2 | NO |
| summary target=pate ppd (summary.R:88-108, default) | per-draw CATE + mean sigma^2 (+ binary per-draw var) | d1 (+sigma) | NO |
| summary target=sate (summary.R:118-146) | per-draw SCATE + noise term | d1 | NO |
| quant/hpd intervals (summary.R:50-62,110-165) | draws-x-1 sequence of averages | d1/dG | NO |
| group.effects (summary.R:269-314) | per-draw GROUP averages | dG | NO |
| fitted/extract(pate/sate/cate) (generics.R:392-409) | per-draw average | d1/dG | NO |
| common support sd.obs/sd.cf (responseFit.R:81-84) | per-obs posterior SD | run/n1 | NO (stored at fit time) |
| plot_indiv, plot_support (plot.R:70-83,146) | per-individual posterior mean | n1 | needs matrix OR n1 |
| extract/fitted(icate/ite/mu.*/y.*) (generics.R:411-440) | individual draws | FULL | YES |
| extract(p.score/p.weights) (generics.R:422-423) | individual p.score draws | FULL | YES (treatment surface) |
| predict.bartcFit (generics.R:87-248) | individual, from newdata | trees | already needs keepTrees |
| refit (new commonSup) (generics.R:503) | full icate re-subset | FULL | YES |
| method.rsp=tmle/p.weight (responseFit.R:297-462) | per-obs IF over per-draw mu0/mu1/p.score | FULL + cross-model | YES (see gap) |

Verdict for (a): the ONLY outputs that genuinely require the n x draws surfaces
to survive the fit are the individual-returning ones (extract/fitted/predict of
icate/ite/mu.*/y.*/p.score/p.weights), refit under a changed common-support rule,
and the tmle/p.weight influence-function estimators. plot_indiv/plot_support want
n1 posterior means; common support is already n-vector-reduced and stored.
Everything print/summary shows is d1/dG - fully callback-reducible.

## The callback contract: what it sees vs needs (gap analysis)

What a reducing callback must compute per draw, for method.rsp="bart":
  mu.obs_i and mu.cf_i (response scale) -> icate_i = (mu.obs_i - mu.cf_i)*(2z_i-1)
  -> emit mean over {ate: all; att: z>0; atc: z<=0} and, for var.exp, the
  within-draw variance; for group.effects, per-group means; for sate, use y_i in
  place of mu.obs_i with a PPD term; sigma is available (aux row / getSigma).

Does the existing callback see enough?
  - yhat_train, yhat_test: the BART components for the observed and counterfactual
    designs. Since bartCause flips only the treatment (a BART covariate) and keeps
    the parametric design identical, indiv.fixef and indiv.ranef CANCEL in the
    Gaussian mu-difference: icate = yhat_train - yhat_test directly. VERIFIED: a
    callback emitting mean(yhat_train - yhat_test)*trtSign matched the post-hoc
    reduction to max abs diff 0.
  - stan_pars (named) + callbackEnv (frame/X/reTrms/test$*/X_means/weights): enough
    to RECONSTRUCT the full ev when the cancellation does not hold. test-11 does
    exactly this and matches extract("ev") bit-for-bit.
  - So for GAUSSIAN CATE/ATE/ATT/ATC the existing seam is already sufficient.

The gaps (what is missing or awkward), each NAMED:
  1. BINARY / nonlinear link. bartCause needs mu on the PROBABILITY scale
     (extract does pnorm, generics.R:541-542; getPATEEstimate.bart.ppd's binary
     branch assumes mu in (0,1)). pnorm(bart+parametric) does NOT cancel the
     parametric offset, so BART-only differencing is wrong; the callback must form
     the FULL latent mean = yhat + fixef + ranef for BOTH designs, pnorm each,
     then difference. It CAN (callbackEnv has X/reTrms), but this is the re-derive-
     ev path, not the cheap cancellation. Contract clarity needed: document that
     yhat_* are pure BART (latent scale for binary), and optionally hand the
     callback the assembled response-scale ev (see Options b).
  2. yhat_train offset-subtracted vs yhat_test not. init.cpp:693-694 subtracts the
     offset from train only; test-11 shows yhat_test is nonetheless the pure BART
     test component (its ev reconstruction matches), so in practice both are BART-
     only and symmetric - but this is undocumented and fragile. Pin it in the
     contract.
  3. SATE/PATE need y and PPD noise. y and weights are in callbackEnv (frame), so
     SATE's y - mu.cf is computable; but PPD noise drawn IN the callback consumes
     RNG mid-sampler -> not bitwise-equal to bartCause's post-hoc sampleFromPPD
     (which replays a stored seed, generics.R:357-376). Statistically equivalent;
     for the AVERAGE the noise is one aggregate N(0, sigma^2*.../n) per draw.
  4. Common-support cutting is fundamentally incompatible with a single-pass per-
     draw average. The subset depends on sd.obs/sd.cf, which are not known until
     all draws are in; a per-draw scalar already averaged over all obs cannot have
     arbitrary observations removed afterward. The callback CAN accumulate running
     per-obs moments (Welford in callbackEnv) to produce sd.obs/sd.cf, but the
     average it emits is fixed to the pre-committed subset (all / att / atc). So
     the callback path serves commonSup.rule="none" (the default) exactly; a
     requested cut needs the surfaces (or accept no re-subsetting).
  5. TMLE/p.weight cannot move to the response callback. They need the treatment
     model's per-draw p.score (responseFit.R:418-423,700-708) aligned to the
     response draw; the two models are independent MCMC runs whose draws are only
     index-matched post hoc, and the callback runs inside the response sampler with
     no access to the treatment chain. Plus the iterative fluctuation is per-draw
     expensive. These stay on stored surfaces.

## The storage mode: callback as the only individual-level record

Design: bartCause constructs a closure that closes over trt/y/weights/group.by/
the design pieces via callbackEnv, and per draw emits a small named vector:
  {ate, att, atc} (or the estimand actually requested), plus optional rows for
  var.exp within-draw variance, per-group averages (dG), and a running-moment
  side-accumulation in callbackEnv for sd.obs/sd.cf. stan4bart stores this as
  result$callback (K x draws x chains, generics.R:427-431); bartCause reads it via
  extract(fit,"callback") and feeds summary/print/fitted.

Compose with store="trees" (not store="fits"): keep the trees C-side so the rare
extract/fitted/predict of individual icate/mu.*/y.* recomputes on demand
(recompute_bart_block, generics.R:253) - the common summaries come from the
callback, the occasional individual ask replays trees. The recompute path's two
historical blockers are both fixed and gated (bart-train-recompute.md,
Landings), so this composition is safe at the scales where the surface
memory actually bites; the residual judgment is latency-only - "keep
surfaces when the user will ask for individuals repeatedly" remains the
right guidance in the Rd.

keep_fits=FALSE + callback (already working, test-11:76-99) is the extreme: no
surfaces, no trees - only the callback array survives. Use when the user wants
only averages and never individuals.

## Memory arithmetic (bartCause-realistic; TWO surfaces)

bartCause's response fit stores mu.hat.obs (train) AND mu.hat.cf (counterfactual
test) - both n x draws x chains doubles - so it pays the sample-storage.md
bart_train cost TWICE, plus the treatment model's n x draws p.score.

Per surface = n x (draws*chains) x 8 B. 4 chains x 1000 sample draws = 4000 draws
(warmup now default-off, sample-storage.md 275853e).

  scale            | one surface | obs+cf surfaces | callback (K=8) | reduction
  -----------------|-------------|-----------------|----------------|----------
  n=10,000         |   320 MB    |     640 MB      |    0.25 MB     |  ~2500x
  n=50,000         |   1.6 GB    |     3.2 GB      |    0.25 MB     |
  n=100,000        |   3.2 GB    |     6.4 GB      |    0.25 MB     |

The callback array is draws*chains*K*8 = 4000*8*8 ~ 0.25 MB, INDEPENDENT of n.
The p.score treatment surface (another 320 MB - 3.2 GB) reduces to its posterior
mean (n x 1) unless extract(p.score)/tmle/p.weight is needed. So for the default
method.rsp="bart", ate/att/atc, commonSup="none" workflow, the callback replaces
~0.6-6.4 GB of surfaces with ~0.25 MB.

## Overhead and correctness (measured)

Measured on the installed stan4bart 0.0.14 (callback + keep_fits + store present),
n=5000, 500 sample draws, 1 chain, random intercept, reducing callback that
reconstructs ranef via Matrix::crossprod and emits ate/att/atc
(/Users/vdorie/.claude/jobs/b073bb28/tmp/overhead2.R):
  - Overhead: +1.1% wall time (median of 5 interleaved reps), ~0.03 ms/draw. The
    callback's R cost is largely offset because keep_fits=FALSE skips the storage
    writes. Even a full ev reconstruction per draw is cheap relative to a WALNUTS
    + BART iteration.
  - Exactness: the per-draw ate from the callback vs the post-hoc reduction of the
    stored ev surfaces = max abs diff 0 (BITWISE) for the mu-level (CATE) path -
    better than store="trees" recompute (~1e-13). SATE/PATE lose bitwise equality
    only through PPD noise ordering (statistically equivalent).
  - Storage: keep_fits=FALSE + callback -> bart_train/bart_test/stan all NULL,
    object 1.4 MB vs 7.9 MB at n=2000 (overhead.R).
  - Warmup: the callback fires during warmup too (init.cpp:714, no guard);
    bartCause should discard result$warmup$callback (or a `callback_warmup=FALSE`
    convenience should gate it). Warmup callback overhead is real but small.
  - Multi-chain: chains run as independent processes; each emits its own
    K x draws block, stacked to K x draws x chains (stan4bart.R:427-431).
    Merging = treat chains as extra draws (combineChains), identical to every
    other per-draw quantity - no special handling.
  - C0 harness relationship: the callback computes the SAME reduction on the SAME
    full-precision per-draw ev the harness compares; for CATE it is exact, so it
    introduces no new tolerance. Only PPD-noise estimands (sate/pate) need a
    statistical, not bitwise, gate - which the harness already applies to draws.

## Options with costs

(a) Do nothing new; document that the existing callback + keep_fits=FALSE is the
    mechanism, and let bartCause write the reducing callback entirely on its side
    (reconstructing ev in R via callbackEnv, as test-11 does). Cost: bartCause M;
    stan4bart 0. Risk: the ev reconstruction and the yhat_train/yhat_test offset
    asymmetry are undocumented contract - fragile across stan4bart versions.

(b) Add a thin stan4bart convenience so the callback does not re-derive ev in R:
    optionally hand the callback the assembled response-scale ev for train and
    test (pnorm applied for binary), i.e. what extract() computes internally, as
    two extra named arguments or a flag. Cost: stan4bart S (reuse the extract ev
    assembly at the per-draw point). Benefit: bartCause's callback becomes a pure
    reducer (difference, subset, average) with no design-matrix algebra, and
    binary is correct by construction. This removes gaps 1 and 2.

(c) Formalize the reducing-callback contract in stan4bart docs (what yhat_* are,
    the scale, the offset convention, warmup firing, callbackEnv contents) and add
    a `callback_warmup=FALSE` gate. Cost: stan4bart S. Pairs with (b).

(d) bartCause side: a constructor `makeAverageCallback(estimand, group.by,
    weights, style)` returning the closure, plus routing getATEEstimates
    (summary.R:167-317) to read extract(fit,"callback") when the fit carried one,
    falling back to the surface path otherwise. Keep the surface path for
    individual-returning outputs, refit, common-support cutting, tmle/p.weight.
    Cost: bartCause M.

## Recommendation

Do (b)+(c)+(d): a small stan4bart contract refinement (hand the callback the
assembled ev + document the seam + gate warmup) and a bartCause reducing-callback
constructor that the DEFAULT average-estimand path uses. Gate bartCause's use on:
method.rsp="bart" AND estimand in {ate,att,atc} AND commonSup.rule="none" AND the
user is not going to request individual draws. Under those conditions store only
the callback array (K x draws x chains, ~0.25 MB) and drop the obs/cf surfaces,
composing with store="trees" so a later extract/fitted/predict of individuals
recomputes from trees (once the recompute blockers are cleared; until then keep
surfaces whenever an individual ask is anticipated). Keep the surface path
verbatim for: any individual-returning extract/fitted/predict, refit under a
changed common-support rule, common-support cutting, and tmle/p.weight.

This captures the author's "callbacks for sample-average estimators" vision
exactly where it is sound (the common default), at ~2500x memory reduction and
~1% overhead, with bitwise-exact averages, while being honest that the individual-
effect use the author calls unavoidable genuinely is - it just moves from an
always-stored surface to an on-demand tree replay.

## Commit partition sketch

C0  Contract doc + test: pin what yhat_train/yhat_test are (pure BART, latent for
    binary, offset convention), stan_pars naming, callbackEnv contents, warmup
    firing; a test asserting a reducing callback's ate == post-hoc reduction
    (bitwise) and keep_fits=FALSE stores only the callback. stan4bart, S.
C1  ev-to-callback convenience (Option b): optionally pass assembled response-scale
    ev(train,test) to the callback; `callback_warmup=FALSE` gate. stan4bart, S.
C2  bartCause makeAverageCallback constructor (ate/att/atc, group dG rows, var.exp
    within-draw variance row, running sd side-accumulation). M.
C3  bartCause routing: getBartResponseFit passes the callback + store="trees"/
    keep_fits under the gate; getATEEstimates/fitted read extract(fit,"callback")
    when present; surface fallback otherwise. M.
C4  Docs + a bartc arg (e.g. store="callback"/"surfaces") + the per-estimand
    capability table for users. S.

## What changes in bartCause / what CANNOT switch (per estimand)

CAN switch to callback (store draws x K, drop surfaces):
  - ate/att/atc point + interval, pate (ppd and var.exp), sate, cate, and their
    quant/hpd intervals - all d1.
  - group.effects averages - dG.
  - common support READOUT sd.obs/sd.cf via running moments (rule stays "none" for
    the average; see below).

CANNOT switch (must keep surfaces or trees):
  - extract/fitted(icate/ite/mu.obs/mu.cf/mu.0/mu.1/y.cf/y.0/y.1) - individual
    draws by definition; served by store="trees" recompute or retained surfaces.
  - predict.bartcFit - already tree-based (keepTrees), unaffected.
  - refit under a NEW common-support rule - needs full icate to re-subset.
  - common-support CUTTING of an average (rule != "none") - post-hoc data-dependent
    subset a single-pass average cannot honor.
  - extract(p.score)/p.weights and method.rsp in {tmle,p.weight} - need the
    treatment model's per-draw p.score (separate, unaligned fit) and per-obs
    influence functions.

## Open questions for the author

Q1. Default vs opt-in in bartCause. Make store="callback" the default when the
    gate conditions hold (method.rsp="bart", ate/att/atc, commonSup="none"), auto-
    falling back to surfaces the moment an individual-returning extract is called
    on a fit that dropped them? That auto-fallback needs store="trees" to be
    robust at scale - which it now IS (both recompute blockers fixed and
    gated; bart-train-recompute.md, Landings). Fork: (i) default callback +
    trees, now unblocked; (ii) default callback + keep surfaces when
    "individuals likely" is signaled by the user, no dependency but coarser;
    (iii) opt-in only. Cost: (i) waits on 2 dbarts-adjacent fixes; (ii)/(iii) ship
    now, less automatic.

Q2. ev convenience (Option b) vs pure-R reconstruction (Option a). Handing the
    callback the assembled ev makes binary correct by construction and removes the
    design-matrix algebra from bartCause's closure, at the cost of a small
    stan4bart-side coupling to the extract assembly. Or keep the callback minimal
    (raw yhat + stan_pars) and make bartCause own the reconstruction (works today,
    per test-11, but re-implements pnorm/offset logic that could drift from
    extract). Which coupling do you prefer to own?

Q3. Common support under callback storage. Options: (i) forbid commonSup.rule !=
    "none" with callback storage (simplest, matches the default); (ii) run a
    cheap pre-pass to fix the subset before the real fit, then the callback
    averages over the frozen subset (two fits, or a warmup-derived subset -
    statistically awkward); (iii) always keep surfaces when a cut is requested.
    (iii) is the honest default; is losing callback storage for common-support
    users acceptable, given they are the memory-hungry large-n case?

Q4. SATE/PATE PPD noise. Draw the aggregate per-draw noise IN the callback (one
    N(0, sigma^2*.../n) variate, cheap, breaks bitwise repro of pate/sate) or keep
    sate/pate on the surface/tree path and let only the mu-level cate/ate family
    use the callback? The former is a wider capture; the latter preserves the
    stored-seed bitwise reproducibility bartCause has today.

Q5. Treatment-model surface. Worth a parallel reducing callback on the propensity
    stan4bart fit to drop its n x draws p.score to a posterior mean when only the
    mean is used (method.rsp="bart", p.scoreAsCovariate)? Or leave it, since it is
    one surface vs the response's two and is needed whole for tmle/p.weight/
    extract(p.score) anyway?
