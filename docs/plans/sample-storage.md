# sample-storage

scope: DESIGN MEMO, no implementation. Measures what stan4bart keeps in the
  returned fit object per posterior draw, verifies the author's recollection
  that full UNTRANSFORMED parameters are stored though consumers only read the
  transformed scale, ranks the storage components at realistic scales, and
  prices the slimming options. The BART engine and dbarts.h are FROZEN (the C
  API is append-only); options that would touch them are flagged dbarts-side.

budget: memo only. If pursued, the recommended parametric-store options (a)+(b)
  are S-M and R+C++ local (walnuts_sampler.cpp / parametric_model.hpp row
  emission + one offset recompute); the real memory lever (c, bart_train) is
  M (warmup-default-off, mostly R) to L (recompute-from-trees / float storage,
  architectural, part dbarts-side).

## Goal

Answer three questions with measured numbers, not recollection:
1. What lands in the fit object per draw, and is the raw/unconstrained block
   actually stored (author's recollection)?
2. At realistic scales, what fraction is the parametric ("stan") store vs the
   BART-side stores? The slimming only pays if its target is a meaningful share.
3. For each slimming option: bytes saved by formula and at scale, consumer cost,
   implementation size.

## Verdict up front

The recollection is CORRECT but the target is small. Every draw stores the raw
unconstrained parameters (z_beta, z_b, z_T) plus the intermediate constrained
ones (rho, zeta, tau, aux_unscaled) at full length, exactly duplicating the
footprint of the transformed block (beta, b, theta_L, aux) that consumers read.
Dropping the raw block halves the parametric store. BUT the parametric store is
only 1-4% of the object; `bart_train` is 96-99%. So option (a) alone saves
~0.5-2% of total memory. The memory lever is `bart_train`, not the raw params.

## Context read in code (file:line)

Per-draw parametric store, WALNUTS row layout (walnuts_sampler.cpp:250-254, 283-293):
  - `sample_names = {"lp__","accept_stat__"}` + `sampler_names = {"stepsize__",
    "treedepth__","n_leapfrog__","divergent__","energy__"}` => 7 diagnostic rows
    at the top of every row; `sample_writer_offset = 7`.
  - `run()` calls `model.eval(position, logp, grad, row + 7)`: `position` is the
    UNCONSTRAINED point, and `eval` writes the CONSTRAINED write_array block into
    `row+7`; `row[0..6]` get the diagnostics (5 of the 7 are constant-zero
    placeholders: accept_stat__=1, treedepth__=0, n_leapfrog__=0, divergent__=0,
    energy__=0; only lp__ and stepsize__ carry a value).
  - constrained block order = `constrainedParamNames()` (parametric_model.hpp:173-191)
    = z_beta(K), z_b(q), z_T(len_z_T), rho(len_rho), zeta(len_concentration),
    tau(t), [aux_unscaled, aux] (if !binary), beta(K), b(q), theta_L(len_theta_L).
    The emission that fills them is eval()'s `if (constrained != nullptr)` block
    (parametric_model.hpp:355-371): it writes z_beta/z_b/z_T straight from the
    unconstrained position, then rho/zeta/tau/au/aux, then beta/b/T.
  - store = `double_writer` (double_writer.hpp:30-50): a flat
    `new double[num_pars * num_samples]`; `createStanResultsExpr`
    (parametric_control.cpp:132-151) memcpys it verbatim into the R `stan` matrix
    with those row names. No filtering: raw rows cross into R.

BART-side per-draw stores (bart_util.hpp:18-58, bart_util.cpp:9-76):
  `IterableBartResults` holds `trainingSamples(numObservations*numSamples)`,
  `testSamples(numTestObs*numSamples)`, `variableCountSamples(numPredictors*
  numSamples)` (uint32 -> R integer, 4 B), `sigmaSamples(numSamples)`,
  `kSamples` (only if k is modeled; default not). `createBartResultsExpr` emits
  the list `sigma/train/test/varcount[/k]`.

Storage sizing + keep switches (init.cpp:592-603, 672-673, 778-779, 785-820):
  `numStorageSamples = keepFits ? numIter : 1`; both the BART results and
  `sample_writer.resize(num_pars, numStorageSamples)` are sized to it; pointers
  advance only `if (keepFits)`. `resultsType` (BOTH/BART/STAN) can suppress one
  side; the R driver always asks "both".

Warmup is stored too (stan4bart_fit.R:51-55; stan4bart.R:325-390, 415-416):
  the worker runs warmup and sample as two `keep_fits` runs and
  `package_samples` builds `result$warmup$*` arrays. With default warmup=iter/2,
  the warmup stores equal the sample stores -> the object is ~2x the sample
  draws. Warmup is read only via `include_warmup=` (opt-in; off by default).

What R consumers actually read (generics.R):
  - `as.array`/`as.matrix`: `grep("^(?:gamma|beta|b|aux)\\.")` (generics.R:32) +
    Sigma derived from `theta_L.` (generics.R:70-99).
  - `extract`: fixef `^beta|gamma` (226,234), ranef `^b\.` (227,242),
    Sigma `^theta_L\.` (228,263), sigma `aux.1` (229,277); `fitted_random`
    reads `^b\.` (generics.R:557). `fitted`/`predict` sit on top of these.
  - NONE of these read z_beta/z_b/z_T/rho/zeta/tau/aux_unscaled. The ONLY surface
    that surfaces the raw rows is `extract(fit, "stan")` (generics.R:283-284),
    which returns `object$stan` whole, and `as.matrix` on the fit does NOT (it
    filters). The raw names appear in R exactly once more, at stan4bart.R:243,
    to compute `upar_names` for `check_sampler_diagnostics` -- which is a
    documented no-op (stan4bart.R:255-269). So the raw block is dead weight to
    every exported computation; only a literal `extract(fit,"stan")` dump sees it.

## Per-draw formula (verified empirically)

Parametric "stan" store, doubles per draw (num_pars):

    num_pars = 7                                   # diagnostics (5 are const-0)
             + K + q + len_z_T + len_rho + len_concentration + t + [!binary]  # RAW/intermediate: z_beta,z_b,z_T,rho,zeta,tau,aux_unscaled  (never read)
             + [!binary] + K + q + len_theta_L     # READ: aux(sigma), beta, b, theta_L(->Sigma)

  where [!binary] = 1 for gaussian, 0 for binary. Store bytes =
  num_pars * n_samples * n_chains * 8, and an equal block again for warmup by
  default. The RAW block and the READ block have the SAME dominant terms (K, q),
  so raw ~= read: dropping raw halves (num_pars - 7).

Empirical check (n=800, K=2, q=40, t=1, random intercept; walnuts build):
  dim(fit$stan)[1] = 95 = 7 diag + 44 raw + 44 read; formula predicts 95.
  Package's own split agrees: par_names$upar = 44, par_names$tpar = 44.
  object.size(fit$stan)/(prod(dim)*8) = 1.047 (dimnames overhead, ~0 at scale).

BART-side per draw: bart_train = n doubles; bart_test = n_test doubles (0 with no
  test data); bart_varcount = n_predictors int (4 B); sigma = 1; k = 1 if modeled.

## Measurements (projected to 4 chains x 1000 SAMPLE draws)

Fits run at reduced iters (storage is exactly linear in draws x chains, verified
at 8.00 B/elt for bart_train, ~8 B/elt for stan) and projected. Scale (a):
n=10,000, 50 groups, (1|g), K=2 -> num_pars=115. Scale (b): n=50,000, two 500-level
factors, (1|g1)+(1|g2), K=2 -> num_pars=2017. No test data; k not modeled.

  component        | scale (a) MB |  %  | scale (b) MB |  %
  -----------------|--------------|-----|--------------|-----
  bart_train       |    305.2     |98.8 |   1525.9     |96.1
  stan (param)     |      3.5     | 1.1 |     61.6     | 3.9
  bart_varcount    |      0.06    | 0.0 |      0.08    | 0.0
  test-fit store   |      0       |  -  |      0       |  -   (none unless test= given)
  -----------------|--------------|-----|--------------|-----
  TOTAL (sample)   |    308.8     | 100 |   1587.5     | 100
  + warmup (deflt) |   +308.8     |     |  +1587.5     |      (warmup stores == sample stores, measured ratio 1.00)
  = draws in object|    617.5     |     |   3175.0     |
  retained data*   |    ~5        |     |    ~30       |      (frame/X/bartData/reTrms; fixed, not draw-scaled)

  *measured at the fit's small draw count (object.size(whole) - draw stores);
   fixed size, so a ~1-2% share once draws reach 1000.

Within the stan store, the RAW block is:
  scale (a): 54 of 115 rows = 47%  -> 1.65 MB sample (3.30 MB with warmup) = 0.53% of total
  scale (b): 1005 of 2017 rows = 50% -> 30.67 MB sample (61.3 MB w/ warmup) = 1.93% of total

## Ranking

At both scales `bart_train` dominates the object (96-99%); the parametric store
is 1-4% and its raw half is 0.5-2% of the total. The author's slimming target is
real (raw params ARE stored, ~half the parametric store) but it is NOT where the
memory is -- `bart_train`, and the default-kept warmup copy of it, are.

## Options (cost priced at scale (b))

(a) Store transformed-only in the stan store. Drop z_beta(K), z_b(q), z_T(len_z_T),
    rho(len_rho), zeta(len_concentration), tau(t), aux_unscaled([!binary]); keep
    aux(sigma), beta, b, theta_L -- which is complete for every consumer (Sigma is
    derived from theta_L; rho/zeta/tau are only its unreconstructed ingredients).
    - Saved: the RAW block = K+q+len_z_T+len_rho+len_concentration+t+[!binary]
      doubles/draw. Scale (b): 1005 * 1000 * 4 * 8 = 30.7 MB sample; x2 warmup =
      61.3 MB (~1.9% of the object). Scale (a): 3.3 MB (~0.5%).
    - Consumer cost: `extract(fit,"stan")` stops returning the z_*/rho/zeta/tau/
      aux_unscaled rows -- the only visible change. No internal consumer, no test,
      and not as.matrix/as.array/fixef/ranef/sigma/Sigma read them. The C++ gradient
      gate (tests/testthat/test-12) builds its own unconstrained vector and is
      unaffected. `par_names$upar` would go empty (feeds only the no-op check).
    - Impl: S-M. In eval()'s constrained emission (parametric_model.hpp:355-371)
      write only {aux,beta,b,theta_L}; trim constrainedParamNames() (180-189); and
      shift betaConstrainedOffset()/auxConstrainedOffset() (160-167) that
      getParametricMean/getSigma index (walnuts_sampler.cpp:301-313). One coherent
      offset change; no BART/ABI impact.

(b) Drop the diagnostics rows. All 7 are unread by any exported surface (5 are
    constant-zero placeholders kept only for Stan-era row names; check_sampler_
    diagnostics is a no-op). Keep lp__ optionally for user MCMC inspection.
    - Saved: 7 doubles/draw = 0.21 MB sample, 0.43 MB with warmup (constant in the
      model size). ~0.01% of the object. Essentially free but essentially nothing.
    - Consumer cost: none (name-filtered everywhere). Impl: S. Low value alone;
      fold into (a) since both edit the same emission + name list.

(c) BART-side -- where the memory is. bart_train is 96-99%.
    - keep_fits=FALSE: EXISTS (init.cpp:595). Sets storage to 1 draw, so
      bart_train/stan are not accumulated -- usable only with a reducing callback;
      not a transparent slimming (fitted/extract then have nothing to read).
    - Thinning: EXISTS via `skip`/n.thin -- reduces n_samples linearly.
    - Warmup-default-off: NEW, the biggest transparent win. Warmup fits are stored
      whenever keep_fits && warmup>0 and read only through opt-in include_warmup.
      Not storing them by default halves the object: scale (b) saves ~1587 MB (50%).
      Impl: M, mostly R (have the worker store only the last warmup draw, or gate
      warmup accumulation on a keep_warmup flag); changes include_warmup default
      availability. No dbarts change.
    - Recompute bart_train from kept trees: NEW, deepest win. With keepTrees the
      trees live C-side and `predict` already recomputes fits (used for test data,
      generics.R:669); calling it on the TRAINING X would let bart_train be lazy
      instead of stored, removing the 96-99% component. Cost: requires keepTrees=TRUE
      (default FALSE), stored-tree memory, and recompute time on access; changes the
      fit-object contract. Impl: L; uses existing dbarts.h predict, so NO ABI change,
      but architectural on the R side.
    - Single-precision bart_train: halves the dominant component but R has no float
      vector and dbarts_results.train is double* (ABI, append-only). dbarts-SIDE,
      L; not recommended.

(d) Never-read blocks to delete outright: same as (a)+(b). `par_names$upar` and the
    no-op check_sampler_diagnostics are dead (names only, negligible bytes) and can
    go with (a). No other stored block is fully unread (warmup and bart_test are
    opt-in, not dead).

## Recommendation

Do (a)+(b) as an always-on, zero-consumer-cost correction: the raw + diagnostic
rows are genuinely never read by any computed surface, and storing beta/b/aux/
theta_L is complete. Total savings at scale (b): (1005+7) doubles/draw x 1000 x 4
x 8 x 2 (warmup) = ~61.8 MB, i.e. ~1.9% of the ~3.2 GB object. Be clear-eyed that
this is small.

If runtime memory is the actual objective, (a)+(b) are not enough -- pursue (c):
warmup-default-off is a clean ~50% cut (scale (b): ~1.6 GB) and is R-local;
recompute-bart_train-from-kept-trees is the only route that attacks the dominant
component itself. Neither requires a dbarts.h change.

Sequence: land (a)+(b) (cheap, correct, removes the author's specific concern),
then decide (c)-warmup as the real lever.

## Open questions for the author

1. extract(fit,"stan") is the one surface that would change under (a): it would no
   longer return z_beta/z_b/z_T/rho/zeta/tau/aux_unscaled. Is that dump a supported
   contract anyone relies on (raw MCMC diagnostics, external tooling), or an
   internal escape hatch we can narrow to the transformed rows? If it must keep raw
   rows, (a) is off the table and the parametric store cannot be slimmed at all.

2. Given raw-param slimming buys only ~2% while default-kept warmup is ~50%: is the
   real ask "stop storing untransformed params" (do (a), small) or "make the fit
   object smaller" (do (c)-warmup, large but changes include_warmup's default
   availability and any workflow that inspects warmup draws)? These are different
   projects; (a) does not move the needle on total footprint.

3. Would you accept keepTrees becoming the default so bart_train can be recomputed
   on demand rather than stored? That is the only option that touches the 96-99%
   component without a precision change, but it trades ~all of bart_train's memory
   for stored-tree memory + recompute latency on every fitted/extract call, and
   changes when predict is available. Worth it, or is stored bart_train the
   intended fast-access design?
