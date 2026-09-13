# re-scale-mixing

The random-effect standard deviation is the one parametric scalar the WALNUTS
block does not mix. This file records why, what was measured, what landed, and
what is still open.

## The bar

dbarts is retiring `rbart_vi` in favour of stan4bart, and set a release bar for
the hand-off, stated on one design: Gaussian random intercept, n = 2000, K = 20
equal groups, Friedman f over 10 predictors, tau = 1, sigma = 1, one chain,
1000 warmup and 1000 kept draws, 200 trees. On that design the random-effect
standard deviation must reach lag-1 autocorrelation below 0.8 and ESS at least
100 per 1000 kept draws.

Reproduced at 979a91a against dbarts 8931adac: lag-1 autocorrelation 0.946 to
0.974 and ESS 2.1 to 28.3 over five seeds. The forest and the group contrasts
mix acceptably over the same draws (contrast ESS 106-148). dbarts's conjugate
Gibbs draw of the same quantity on the same data reaches lag-1 autocorrelation
0.14-0.47.

Every figure in this file was re-measured against dbarts 8931adac, which is
past the commit that retired `rbart_vi`. The draws are bit-identical to the
recording made against c585ba2c, so that retirement is draw-neutral here.

## Diagnosis

The scale is sampled on the log scale already (`ParametricModel::eval`,
`tau = exp(tau_free)`) and the block is already non-centered
(`b_level = T_i * z_b_level`), so neither spelling is the fault. Four candidates
were measured.

**The parameterization (the fault).** With 100 observations per group the effects
are pinned by the likelihood, so the non-centered product `b = tau * scale *
dispersion * z_b` puts the log scale and the standardized effects on a ridge:
over the kept draws, `cor(log tau, log sum z_b^2) = -0.978` (seed 7, the figures
below are that chain). The ridge is exact by construction rather than incidental.
Every block's Cholesky factor is homogeneous of degree one in its own scale -
`bl.s = tau[i] * re_scale[i] * dispersion` in `ParametricModel::eval`, and every
entry of `T` is a multiple of `bl.s` - so scaling `tau_i` by r and dividing that
block's `z_b` by r leaves `b`, and therefore the linear predictor and the
likelihood, untouched. Only the priors and the Jacobian resist the move.

In the coordinates the sampler actually integrates - the unconstrained position
divided by the square root of its adapted diagonal inverse mass - the block's
covariance has two long directions and twenty short ones:

| direction | length | acf1 | ESS | loading |
| --- | --- | --- | --- | --- |
| 1 | 6.48 | 0.972 | 6 | 0.70 z_b, 0.30 tau |
| 2 | 4.86 | 0.968 | 13 | 0.96 z_b, near-uniform |
| 3-22 | 0.42 to 1.07 | median 0.107 | 8 to 875 | contrasts |

Direction 1 is the ridge; it carries 90 percent of the scale coordinate's
variance and spreads the rest over all twenty standardized effects, so it is a
rotated direction - a diagonal metric can equalize the coordinates and still
leave it. Direction 2 is the shared level of the standardized effects
(correlation -0.93 with `mean(z_b)`), which is a different defect and is treated
below. The sampler mixes every direction of the parametric block except those
two, and it stops short of them because the NUTS U-turn test is dominated by the
twenty short ones: mean 10.5 leapfrog steps at step size 0.364 is an integration
time of 3.8, against a length scale of 6.5.

**Adaptation (cannot help, and is not separately at fault).** The metric is
measurably wrong on the scale coordinate - whitened by the frozen inverse mass,
the tau coordinate has marginal sd 3.76 where a calibrated metric would give 1 -
but that is downstream of the failure rather than upstream of it: the marginal
variance of a direction the chain never crosses is exactly what warmup cannot
estimate. Lengthening warmup confirms it. Warmup 4000 leaves lag-1
autocorrelation at 0.933 on seed 7 and 0.978 on seed 42, with whitened tau sd
2.23 and 4.18 - no convergence toward 1 and no consistent gain. The NUTS doubling
cap is not binding either - the vendored default is 5, and raising it to 6, 7, 8,
10, or 12 reproduces the same draws bit for bit, so trajectories end at a U-turn
well below the cap. `adapt_delta` moves it a little and not enough: 0.99 buys
step size 0.108 and 31 leapfrog steps for lag-1 autocorrelation 0.896.

**The two-block alternation (not the scale's fault - but it is the level's).**
Substituting a fixed forest for the BART block (`offset_type = "bart"` with a
frozen offset) leaves the parametric target static from sweep to sweep. It
separates the two failures cleanly, on one seed and one design:

| quantity | live fit | frozen forest |
| --- | --- | --- |
| scale, acf1 / ESS | 0.963 / 13 | 0.954 / 22 |
| shared level, acf1 / ESS | 0.980 / 4 | 0.093 / 727 |
| sigma, ESS | 13 | 633 |
| group contrasts, ESS | 149 | 815 |

The scale barely moves, so the alternation is not what it is paying. The shared
level and sigma go from unusable to clean, so the alternation is the whole of
what they are paying.

**The prior (not the fault).** Widening `decov()`'s scale from 1 to 10 moves
lag-1 autocorrelation from 0.963 to 0.943; `shape = 0.5, scale = 5` gives 0.961;
concentrating it (`shape = 2, scale = 0.5`) gives 0.953. The prior modulates the
ridge's length and does not remove it.

## What landed

`skip`'s `"stan"` element - documented in man/stan4bart.Rd and reaching
`StanControl::skip` since the port, but consumed by nothing - now takes that many
parametric transitions per BART sweep, keeping the last (`run` in src/init.cpp).
The scale's ridge is the only direction that needs them, and they cost one
parametric transition each and no BART work. `mean_leapfrog` divides by the
transition count rather than the sweep count so it stays per-transition, which is
what its warning threshold means.

The extra transitions apply to the sampling phase only. Under a uniform loop the
parametric block outruns the still-growing forest during early warmup: with
`skip = 8` the shared level of the group intercepts starts warmup at 8.2 rather
than 2.0, the scale tracks it (Sigma is the spread of the effects about zero, not
about their own mean), and the pair needs some 4000 sweeps to relax. On the bar's
design one of two seeds was still at tau_mean 5.3 after 1000 warmup and 1000 kept
draws. Warmup 4000 clears it, and so does confining the extra transitions to
sampling, which is also draw-neutral for warmup at every skip.

At the default `skip = 1` the change is a no-op by construction, and the
posterior baselines and the exactness gate confirm it.

The scale's ridge move landed after it, and is what the default now relies on;
the "Open" section below records where that leaves this file.

## Result against the bar

Bar design, 1000 warmup and 1000 kept, ten seeds:

| skip | acf1 median | acf1 worst | ESS median | ESS worst | tau mean | seconds |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.956 | 0.974 | 13.7 | 2.1 | 0.990 | 2.30 |
| 4 | 0.883 | 0.902 | 54.8 | 9.2 | 1.006 | 2.94 |
| 8 | 0.770 | 0.821 | 114.7 | 29.8 | 1.004 | 3.61 |
| 16 | 0.625 | 0.703 | 222.3 | 133.0 | 1.010 | 5.19 |

`skip = c(bart = 1, stan = 16)` clears the bar on all ten seeds at 2.3x wall
time. `skip = 8` clears it on the median and on eight of ten seeds; its two
misses are one chain at lag-1 autocorrelation 0.821 and one whose lag-1 reads
0.785 but whose ESS estimate reads 29.8. Read the autocorrelation, not the ESS,
when they disagree: at 1000 draws `ess_basic` is the noisier of the two on a
chain with structure past lag 1 (at `skip = 24` the median autocorrelation falls
to 0.543 and one seed's ESS still reads 17).

The scale's posterior mean is right at every skip, so nothing here is buying
mixing with bias. The other two designs move the same way and how far the skip
has to go tracks how strongly the groups are identified, which is what the ridge
picture predicts. From benchmarks/baselines/re-scale-BASELINE.csv, worst of five
seeds:

| design | skip 1 | skip 8 | skip 16 |
| --- | --- | --- | --- |
| gaussian_k20 (100 per group) | 0.974 | 0.781 | 0.671 |
| gaussian_k5 (400 per group) | 0.994 | 0.964 | 0.930 |
| probit_k20 (100 per group, latent) | 0.938 | 0.697 | 0.482 |

The 5-group Gaussian fit - four times the observations per group, so the longest
ridge of the three - is not reached by any skip worth paying for. Skip buys the
bar's design and no more; the ridge move below is what closes the rest. On the
ESS half of the bar probit_k20 needs 16 as well: its worst seed reads 42.6 at
skip 8 against a worst acf1 of 0.697.

## The shared level does not mix, at any skip

The scale is the quantity the bar names, and it is not the only one failing. The
shared level of the group intercepts - `mean(b)` over the groups - has lag-1
autocorrelation 0.97 to 0.99 and ESS 1.5 to 19 per 1000 kept draws on every
design and at every skip in the baseline. The extra parametric transitions do
not touch it, and neither does the `"bart"` element:

| seed | skip | scale acf1 / ESS | level acf1 / ESS | seconds |
| --- | --- | --- | --- | --- |
| 7 | bart 1, stan 1 | 0.963 / 13 | 0.980 / 3.9 | 2.5 |
| 7 | bart 1, stan 8 | 0.742 / 104 | 0.987 / 2.3 | 3.2 |
| 7 | bart 1, stan 16 | 0.601 / 272 | 0.977 / 20.2 | 4.4 |
| 7 | bart 8, stan 1 | 0.947 / 17 | 0.984 / 2.5 | 15.5 |
| 7 | bart 8, stan 16 | 0.633 / 183 | 0.982 / 5.5 | 18.0 |

That is the expected shape once the frozen-forest table above is read: the level
is confounded with the forest's overall level, so it is the alternation that has
to cross it and no amount of work inside either block will. Eight BART sweeps per
stored draw costs six times the wall clock and buys nothing, which is the
signature of a between-block ridge rather than a within-block one.

The consequence for a user is bounded but real. The group contrasts mix (ESS 76
to 250 across the baseline) and the scale's posterior mean is right at every
skip, so `Sigma` and the relative ordering of the groups are trustworthy. What is
not trustworthy at 1000 draws is the common component of `ranef`: its Monte Carlo
error is roughly `tau / sqrt(K * ESS)`, which on the bar's design is the same
order as the level itself.

The fix is a joint move that shifts a scalar between the forest and the
intercepts, leaving the linear predictor fixed - the location analogue of the
scale's ridge move. It needs a way to add a constant to every leaf of the forest,
which the dbarts flat C API does not expose. Not scheduled; recorded here so the
level columns in the baseline are not read as noise.

## Open

**The scale is drawn along the ridge, and the default is the move, not `skip`.**
The closed-form slice move this file argued for landed: once per sweep, per
random-effect block, during sampling, rebuilding the frozen sampler at the moved
position rather than editing the vendored headers. The ridge is exactly
traversable, so the conditional along it carries no likelihood term and a slice
sampler draws the scale nearly independently at O(q) with no gradient - which is
the direction WALNUTS cannot travel and nothing else needs.

On the bar's design at the DEFAULT `skip`, worst of the harness's five seeds,
lag-1 autocorrelation goes from 0.963 to 0.090 and ESS from 7.0 to 526 per 1000
kept draws, at a wall-time ratio of 1.01 on the bartCause configuration. The
derivation as implemented, the corrections it makes to the design's formula, and
the full bar are in docs/design/re-scale-and-grouped-cost.md, section 6. It is
on by default; `stan_args = list(ridge_move = FALSE)` turns it off.

**`skip` stays the escape hatch, and its default stays 1.** Nothing above about
its cost has changed - a WALNUTS transition is 0.079 ms against a 0.19 ms BART
sweep, so `stan = 16` still roughly doubles the fit - but the reason to reach for
it has narrowed to what the move does not reach. The measured cost of raising it
is still:

| design | skip 1 | skip 16 | ratio | what it buys |
| --- | --- | --- | --- | --- |
| gaussian_k20 | 2.30 s | 4.63 s | 2.0x | scale ESS 14 -> 250 |
| probit_k20 | 2.51 s | 4.80 s | 1.9x | scale ESS 35 -> 227 |
| no random effects (`t == 0`) | 2.32 s | 2.85 s | 1.2x | sigma ESS 65 -> 220 |

The `t == 0` row is the one the move does not touch at all: there is no ridge
there, and the extra transitions still triple sigma's ESS, because sigma pays
the two-block alternation whether or not there is a random-effect block. That
is the surviving case for the escape hatch.

**What the move does not reach.** Two things, both of them the two-block
alternation rather than the parametric block's geometry.

The shared level of the group intercepts is the first, and it has its own
section above. The move leaves it exactly where it was - lag-1 0.97 to 0.99 on
every design - which is the confirmation the ridge picture predicted: the level
is confounded with the forest's overall level, and a move that holds the linear
predictor fixed by construction cannot cross it. Its fix is still a joint move
shifting a scalar between the forest and the intercepts, which still needs a way
to add a constant to every leaf that the dbarts flat C API does not expose. Not
scheduled.

Five-group designs are the second. On `gaussian_k5` the move takes worst-seed
lag-1 from 0.986 to 0.150 - so the ridge is crossed there as exactly as
anywhere - but one of the five seeds reads an ESS of 10.5 against that lag-1,
and on the group-sd harness's `few_large` case (five groups of four hundred,
four chains) the group sd clears the autocorrelation half of the bar on all
three seeds and the effective-sample-size half on none, one of them at an R-hat
of 1.82. With five groups the scale is weakly identified and the four chains
disagree about it; what remains is between-chain spread, not within-chain
autocorrelation, and the ESS estimator charges for the first. That is a
different defect from the one this file diagnoses and it is not the ridge.

**Warmup is unchanged.** The move applies after the freeze. The second local
patch this file held in reserve - a position setter on `AdaptiveWalnuts`, so a
warmup sampler could be moved without resetting Adam and the mass estimator -
was not needed: the sampling-only move clears the bar on every seed, so warmup
adaptation is not what was holding the scale back. The vendored tree still
carries exactly one local patch.

**The centered alternative** - sampling the effects directly and the scale from
its own conditional - is the textbook answer for well-identified groups and is a
rewrite of `make_theta_L`, `make_b`, the hand adjoint, the raw block layout, the
gradient-gate fixtures, and every baseline. It also reintroduces the funnel for
weakly identified designs, so it wants a switch rather than a swap. The move
makes it unnecessary for the quantity the bar names.

## The WALNUTS refresh, as a control

The vendored sampler was refreshed to upstream head at 2fa75b7, four of whose
changes move numerics: an aliased Welford update in the mass estimator, a
doubling count that was double-reported on capped or aborted warmup draws, a
non-finite acceptance statistic that is now a rejection rather than an input to
Adam, and a reused gradient at the selected position. Every draw moved. The
question this section answers is whether the scale's mixing moved with them, so
that later work is measured against the refreshed build and not against a
picture the refresh had already changed.

It did not. The bar's design, one chain, 1000 warmup and 1000 kept, 200 trees,
the same ten seeds on both builds - 20260907, 7, 42, 1234, 99, 2, 13, 314, 2718,
4242, the harness's own five followed by five more. Every column but the last is
a median or a worst over the ten; the seconds column is their mean:

| skip | build | acf1 median | acf1 worst | ESS median | ESS worst | tau mean | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | before | 0.959 | 0.977 | 13.0 | 3.0 | 1.019 | 2.04 |
| 1 | after | 0.960 | 0.993 | 13.3 | 1.3 | 1.040 | 2.03 |
| 16 | before | 0.651 | 0.759 | 176.0 | 10.6 | 1.001 | 4.15 |
| 16 | after | 0.639 | 0.797 | 184.2 | 44.8 | 1.009 | 3.95 |

The medians agree to about a hundredth in autocorrelation and a few percent in
effective sample size; the worst-seed columns move more, in both directions,
which is what a worst-of-ten statistic does when every draw has changed. Mean
leapfrog steps per transition are unchanged at the median, 10.8 before and 10.6
after. At skip 16 they fall from 10.8 to 9.8, which is the largest move in the
table and still not a cost. The ten seeds used here are not the ten the table
above under "Result against the bar" was recorded on, which were not written
down; on this set skip 16 does not clear the worst-seed ESS floor on either
build, so read the two builds against each other rather than against that
table's verdict.

The four-design group-sd measurement behind `docs/mixing-group-sd.md` moves the
same way. Group sd lag-one before and after, three seeds per design:

| design | acf1 before | acf1 after | ESS/1000 before | ESS/1000 after |
| --- | --- | --- | --- | --- |
| bar_reference | 0.93 - 0.96 | 0.95 - 0.96 | 12 - 19 | 10 - 19 |
| many_small | 0.88 - 0.91 | 0.87 - 0.91 | 24 - 39 | 19 - 45 |
| few_large | 0.96 | 0.96 - 0.97 | 1 - 12 | 2 - 21 |
| weak_signal | 0.62 - 0.73 | 0.59 - 0.73 | 18 - 164 | 24 - 92 |

The bar is met in no design on either build. The one fit that cleared it before
- weak_signal at the first seed, on an effective sample size of 164 - reads 24
after, and another seed of the same design goes 18 to 72. Nothing in the
autocorrelation column moves by more than 0.06. That contrast is the useful
finding here: on this parameter the effective-sample-size estimate is not stable
enough at 1000 draws to carry a verdict on its own, and the autocorrelation is.

## Two sampler settings that had never been measured

`max_hamiltonian_error` and `min_micro_steps` reach the sampler through this
package's `SamplingConfig` and had never been varied. Neither is reachable from
R: each probe row is a separate build with the `SamplingConfigBuilder` call in
src/walnuts_sampler.cpp given the setting, so reproducing a row means rebuilding.
Both were probed on the bar's design at the default skip and the harness's own
five seeds - 20260907, 7, 42, 1234, 99 - against the refreshed build, read the
same way as the table above. The expectation was that neither would matter. One
of them does.

| setting | acf1 median | acf1 worst | ESS median | ESS worst | leapfrog | seconds |
| --- | --- | --- | --- | --- | --- | --- |
| default (0.5, 1) | 0.957 | 0.963 | 13.5 | 7.0 | 10.7 | 2.05 |
| max_hamiltonian_error 0.1 | 0.964 | 0.986 | 14.7 | 5.3 | 26.5 | 2.42 |
| max_hamiltonian_error 2.0 | 0.938 | 0.960 | 25.0 | 14.7 | 9.1 | 2.02 |
| min_micro_steps 2 | 0.943 | 0.953 | 25.6 | 20.1 | 14.8 | 2.14 |
| min_micro_steps 4 | 0.821 | 0.893 | 64.4 | 54.5 | 23.9 | 2.37 |
| min_micro_steps 8 | 0.887 | 0.901 | 49.7 | 31.1 | 13.6 | 2.20 |

`max_hamiltonian_error` behaves as expected: tightening it to 0.1 buys nothing
and costs two and a half times the gradient evaluations, and loosening it to 2.0
is a small gain inside seed noise. `min_micro_steps` is not what was expected.
Four micro steps per macro step takes the median autocorrelation from 0.957 to
0.821 and the median effective sample size from 13.5 to 64.4, worst seed 7.0 to
54.5, for 1.16x the wall time - against the 2x that `skip = 16` costs for a
comparable gain. It still does not clear the bar, and it is not a monotone
lever: eight micro steps is worse than four on every column, so this is a single
five-seed probe of a setting with an interior optimum, not a curve with a known
shape.

That is enough to say the setting is worth a proper look and not enough to move
a default on. What it changes for the ridge move is the comparison: the move now
has to beat a cheaper alternative than `skip` alone.

## Harness

benchmarks/R/bench-re-scale.R records the table above over three grouped designs
(20-group and 5-group Gaussian random intercepts, 20-group probit), a skip grid,
and five seeds, and prints the verdict against the bar. It also records the
shared level's own autocorrelation and ESS, so the second defect is visible in
the baseline rather than only here. It replaces the `rbart_vi` comparison the bar
came from, which cannot survive that function's removal from dbarts; the same
removal took the third comparator out of inst/tinytest/test-02-binary.R, whose
surviving glmer and plain-forest arms carry the intent.

Its `ridge` mode is the scale move's own acceptance bar: the default skip, the
move off and on, every design, the five seeds, so the two columns are read
against each other.

The two sections above drive the same harness from a caller that names the
design, the skip and the seeds, rather than through its `record` mode:

    source("benchmarks/R/bench-re-scale.R")
    record_re_scale(outfile, "gaussian_k20", c(1L, 16L), SEEDS_10)

with SEEDS_10 as listed there, and `1L` alone for a probe row. Each row is
deterministic in its seed, so a re-run of any of them on the same build
reproduces it exactly.
