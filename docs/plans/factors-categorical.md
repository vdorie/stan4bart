# factors-categorical

scope: adopt dbarts's categorical splits for factors in the BART
  component, replacing the port-era factors = "indicators" pin.
  DECISION RECORD (VD 2026-07-16): adopt outright - no user argument,
  no indicator fallback. Categorical handling applies ONLY to the BART
  component. Factors in the fixed-effects part keep the standard
  model-matrix contrasts (the parametric design is not dbarts's
  business). Factors in the random-effects structure keep the lme4
  machinery untouched, including its new-level semantics: unknown
  grouping levels at predict remain possible and get prior draws or
  zero-ev random effects per sample_new_levels. Unknown BART-factor
  levels, by contrast, are an error.

budget: S-M. One commit. Draw-moving for factor-bearing fits only
  (a different tree prior over factor structure - one categorical
  column split by level subsets instead of one indicator column per
  level); bitwise-neutral for factor-free fits.

## Context (file:line ground truth)

- The pin: lme4_functions.R:178 `dbarts::dbartsData(bartform, bartfr,
  factors = "indicators")`, with the comment marking it as port
  scaffolding. dbarts >= 1.0 defaults factors = "categorical" in every
  entry point (dbarts R/data.R:307).
- Level bookkeeping: lme4_functions.R:106-109 stores per-factor
  training levels as attr(terms, "levels.bart"); test_data.R:44-53
  re-levels newdata through model.frame(xlev = orig.bart.levs) - a NEW
  level in newdata errors there today (R's standard "factor has new
  levels" error), in both modes. This is the desired policy; keep it
  and test it.
- Test matrix: test_data.R:135 `dbarts::makeTestModelMatrix(
  object$bartData, mf.bart)` - dbarts's own builder over the stored
  bartData; under categorical it codes test columns over the training
  factor.levels attribute. Fit-time test data flows through the same
  frames (stan4bart.R:119).
- The C boundary is NOT affected: src/init.cpp:195/374 pass the
  dbartsData S4 object itself to dbarts_sampler_create; categorical
  metadata rides as attributes on bartData@x exactly as dbarts's own
  front end builds them. No dbarts.h change.
- Varcount surface: package_samples names varcount rows by
  colnames(bartData@x) (stan4bart.R:251); under categorical a factor
  is ONE column, so extract(type = "varcount") reports per FACTOR, not
  per level. n_bart_vars consumers: generics.R:403, 893.
- Name-collision guard: stan4bart.R:51 checks `"b" %in%
  colnames(bartData@x)`; under categorical, colnames are bare variable
  names (no .level suffixes), so a factor literally named "b" now
  trips the same guard a numeric "b" always did. No change needed;
  noted for awareness.

## Design

Flip lme4_functions.R:178 to factors = "categorical" and rewrite its
comment to state the modeling intent (single categorical column split
by level subsets; selection prior sees one variable per factor;
varcounts per factor). Everything else is verification, tests, and
docs:

- Factor-free fits must be BITWISE identical before/after the flip
  (all-numeric frames produce the same matrix under either builder).
  This is the commit's hard gate.
- Factor-bearing fits: draws move by design. Verify structure instead
  of values: varcount has one row per bart variable with the factor's
  bare name; extract/fitted on fit-time test data agree with
  predict-on-the-same-rows (encoding consistency through
  makeTestModelMatrix and the store = "trees" recompute seam, which
  reads bartData@x/@x.test - generics.R:217-218, 246, 258, 723);
  new BART-factor level in newdata errors informatively; new grouping
  level in the random effects still follows sample_new_levels (prior
  draws / zero-ev), unchanged.
- The distributional-equivalence tiers are all-continuous and are
  unaffected; they are not a gate for this commit (the sampler path is
  untouched - this is data encoding).
- Docs: man/stan4bart.Rd notes factor handling in the bart() part
  (categorical splits, per-factor varcount, new-level error). NEWS.md
  entry: behavior change for factor-bearing fits, per-factor
  varcounts, and that fixed/random-effect factor handling is
  unchanged.
- test-04-factor_levels.R expectations reworked for categorical
  (bare-name columns, per-factor varcount, new-level error), plus the
  factor-free bitwise case wherever it fits best.

## Gates

Full testthat suite green (NOT_CRAN=true, max_fails Inf; 338+ grows
with new cases); factor-free seeded fit bitwise identical() to the
pre-flip build (draws, varcount, adaptation$step_size); R CMD check
--no-tests --no-manual --no-vignettes on the built tarball Status OK.

## Landing

Changed: R/lme4_functions.R:178 (the flip, comment rewritten to state
modeling intent); R/generics.R's recompute_bart_block (a fix, see
below); tests/testthat/test-04-factor_levels.R (reworked + new cases);
tests/testthat/test-08-glFormula.R (one hardcoded per-level varcount
name updated to categorical's bare-name convention); tests/testthat/
test-15-store-trees.R (new factor-bearing store = "trees" case);
man/stan4bart.Rd and NEWS.md (factor-handling docs).

Surprise: the flip exposed a real bug, not anticipated by this plan's
Context section. Under `factors = "categorical"`, `bartData@x` is a
`dbartsMixedMatrix` (dbarts's sparse/dense training representation)
INSTEAD OF a plain matrix - not only when a factor is present, but for
every fit, including factor-free ones (`factors = "indicators"` always
produced a plain matrix; confirmed directly against dbarts::dbartsData
with and without a factor column). recompute_bart_block (the store =
"trees" recompute seam) passed `bartData@x` straight into predictBART's
C boundary, which requires a real matrix and errored ("x.test must be
of type real") for EVERY store = "trees" fit's train-sample recompute
under categorical splits - this broke the pre-existing (factor-free)
test-15 suite outright, not just new factor-bearing cases. bartData@x.test
was unaffected (already densified by makeTestModelMatrix), so test-sample
recompute worked throughout. Fixed with a one-line R-side conversion
(`if (inherits(X, "dbartsMixedMatrix")) X <- as.matrix(X)`, dbarts's own
dense-conversion method, encoding factor columns the same 0-indexed-double
way x.test already does) - no src/ or dbarts.h change. Also manually
verified (not a permanent test) that saveRDS/readRDS round-trips a
factor-bearing store = "trees" fit correctly, since data.bart.light@x's
zero-row placeholder (stan4bart_fit.R:677) is also now a MixedMatrix.

Gate results:

- Bitwise factor-free gate: seeded continuous + random-effects fit
  (n=200, chains=2, warmup=5, iter=15, n.trees=10, seed=99), reverted
  the flip, reinstalled, ran and saved as the true "before"; reapplied
  the flip, reinstalled, ran again. `identical()` TRUE component by
  component: $stan, $bart_train, $bart_varcount, adaptation$step_size.
- test-04 (factor-bearing, X1 with empty levels, keepTrees): varcount
  named per factor ("X1", not "X1.a".."X1.e" - confirmed against the
  pre-flip build, which produced "X1.a".."X1.e" for the 5 observed
  levels only); predict(df.test) on fit-time test rows identical() to
  extract(sample = "test"); new BART factor level at predict errors
  with "factor X1 has new levels z"; new random-effects grouping level
  (g.1, a random-slope block) predicts finite values under both
  sample_new_levels settings, with the same pre-existing "longer object
  length is not a multiple of shorter object length" warning test-01
  documents (levelfun, unrelated to this change, left alone).
- test-15 factor-bearing addition: store = "trees" vs store = "fits"
  on a factor-bearing bart() component - $stan and $bart_varcount
  identical(); extract(indiv.bart/ev, train/test) and fitted(indiv.bart/
  ev) equal to the file's existing TOL = 1e-10.
- Full suite: 65 test_that blocks, 358 passed expectations, 0 failures,
  3 warnings (all expected: test-01's pre-existing one plus 2 from the
  new test-04 sample_new_levels loop).
- R CMD check --no-tests --no-manual --no-vignettes on the built
  0.0-14 tarball: Status OK, no NOTEs.

Off-spec noticed, not touched: none found beyond the recompute_bart_block
fix above, which is squarely in-scope (R-only, on the seam the plan's
Design section names).
