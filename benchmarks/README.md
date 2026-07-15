# C0 baselines (docs/plans/walnuts.md)

Nothing here ships with the package (add to .Rbuildignore before C5 if not
already covered). Everything runs against the *installed* stan4bart; install
first with `R CMD INSTALL .`.

## PREREQUISITE: the bart_util.hpp structSize fix

At C0 recording time the tree as of 967a1c6 was found BROKEN against
dbarts 1.0.0: `IterableBartResults::current` (a `dbarts_results`) never set
`structSize`, so `dbarts_sampler_run` skipped every output field
(`DBARTS_RESULTS_HAS` false for all of them), the train-fit buffers stayed
zero, and init.cpp fed Stan `stanOffset = 0 - parametricMean`. The Gibbs
coupling then decoheres: beta random-walks hundreds of prior sds out
(observed z_beta ~ 200 on a y = 2x + noise fixture), the stored "BART" fits
are exactly `-X*beta`, `predict`-vs-`extract` tests fail (test-02, test-06),
and chains progressively lock into permanently-divergent zero-acceptance
states. Zero-initializing alone makes the skip deterministic but does NOT
restore output - the constructor must also set
`current.structSize = sizeof(dbarts_results);`. Both land in the commit
these baselines are named after; on the fixed build the full testthat suite
is green (including the NOT_CRAN accuracy tests) and every chain in every
tier mixes cleanly at default `stan_args`.

## R/record-posterior-baselines.R - the C0(b) distributional oracle

Records fixef/ranef/sigma/Sigma posterior means, sds, and MCSEs from the
Stan-era sampler on five reference datasets spanning the WALNUTS gradient's
tiers: continuous nc=1, nc=2, nc>=3 (the decov "onion"), weighted
continuous, and binary two-block multilevel. Data regenerate
deterministically from fixed seeds (see `TIERS`); no external data files.
MCSEs are hand-rolled (no coda dependency): the max of pooled batch means
(b = n^(2/3)) and the between-chain estimator - the latter is essential
here because the slow BART-vs-parametric tradeoff carries autocorrelation
far past any within-chain batch (observed 3-8x the sqrt(n)-batch estimate
on ranef summaries).

    Rscript benchmarks/R/record-posterior-baselines.R record \
      benchmarks/baselines/posterior-baselines-75d7970.rds

Recording saves after every tier and merges into an existing outfile, so a
subset can be re-recorded. `rerun-check` is the C0 stability gate: refits
one tier's SAME data at a different MCMC seed and reports every gated
summary's mean shift against k=2 times the combined MCSE. With ~100
summaries a handful of ratios in (2, 3) is the expected tail of calibrated
MCSEs, not instability; the verdict line prints the expected exceedance
count.

    Rscript benchmarks/R/record-posterior-baselines.R rerun-check \
      continuous_nc3 4242 benchmarks/baselines/posterior-baselines-75d7970.rds

Gate results at recording (2026-07-15): continuous_nc3 at seed 4242 -
94/98 within 2x (expected ~4.5 exceedances, max ratio 2.85);
binary_multilevel at seed 5151 - 53/55 (expected ~2.5, max 2.33). PASS.

## R/bench-perf.R - the C0(a) perf baselines (written, not yet run)

Measures wall-clock `R CMD INSTALL --preclean` time, peak sampling RSS, and
per-iteration sampling wall-time (differencing two `iter` values to remove
fixed setup cost) on three of the same reference fits (continuous_nc2,
binary_multilevel, weighted_continuous). NOT YET RUN - needs a quiet
machine, a later grant. Re-run identically at C4 and diff by hand.

    Rscript benchmarks/R/bench-perf.R record benchmarks/baselines/bench-perf-75d7970.csv

## baselines/ - recorded output

`posterior-baselines-<commit>.rds`: `$manifest` (commit, package_state -
NOTE: includes the structSize fix, date, package/R versions, MCMC control,
chain-health thresholds) and `$tiers` (per-tier simulated data, formula
shape, healthy-chain mask, summary data frame with mean/sd/mcse/mcse_bm/
mcse_bc per gated quantity). See baselines/MANIFEST for provenance.
