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

**The default is still 1**, so the bar is met by an argument and not by the
package as shipped. A user who does not set `skip` still gets a scale that does
not mix. The measured cost of raising it to 16 is smaller than it looks, and it
is not confined to models with random effects:

| design | skip 1 | skip 16 | ratio | what it buys |
| --- | --- | --- | --- | --- |
| gaussian_k20 | 2.30 s | 4.63 s | 2.0x | scale ESS 14 -> 250, clears the bar |
| probit_k20 | 2.51 s | 4.80 s | 1.9x | scale ESS 35 -> 227, clears the bar |
| no random effects (`t == 0`) | 2.32 s | 2.85 s | 1.2x | sigma ESS 65 -> 220 |

The `t == 0` row is worth noting for a different reason than it first suggests:
there is no ridge there, so the extra transitions should be pure waste, and they
are not - they cost 1.2x and triple sigma's ESS, because sigma pays the two-block
alternation whether or not there is a random-effect block.

Recommendation nonetheless: **leave the default at 1.** The wall-time
decomposition in docs/design/re-scale-and-grouped-cost.md is what decides it. A
WALNUTS transition costs 0.079 ms against a 0.19 ms BART sweep, so `stan = 16`
adds 1.2 ms per sweep and roughly doubles the fit - on the axis that is already
the live complaint about this package against the dbarts route it replaced.
Buying the scale's mixing at 2x wall time is the wrong trade when the same
mixing is available for free along the ridge itself (that design's section 4).
`skip` stays the escape hatch a user can reach for today; the default moves, if
at all, only if the ridge move is abandoned.

**The fix that does not cost 8x**, worked out with its acceptance bar in
docs/design/re-scale-and-grouped-cost.md. The ridge is exactly traversable in
closed form. Every block's Cholesky factor is homogeneous of degree one in its own scale
(`bl.s = tau[i] * re_scale[i] * dispersion`, and every entry of `T` is a multiple
of `bl.s`), so scaling `tau_i` by r and dividing that block's `z_b` by r leaves
`b`, and therefore the linear predictor and the likelihood, untouched. Along that
curve the target is a one-dimensional density in closed form - the standardized
effects' normal prior, the scale's gamma prior, and the change-of-variables
Jacobian, no likelihood term - so a slice sampler draws the scale nearly
independently at O(q) per sweep with no gradient evaluation. It is the direction
WALNUTS cannot travel and nothing else needs.

The cost is that WALNUTS holds its position privately: neither `AdaptiveWalnuts`
nor `WalnutsSampler` exposes a way to write `theta_`. Reading it is already
solved - `LatestDraw::on_sample` is handed `theta_` after every transition, so
the handler holds the current position. The ask on the vendored headers is
therefore ONE accessor, a position setter, not two; adding it ends the "vendored
verbatim at commit 5854be8" claim in LICENSE.note.

The alternative, rebuilding a `WalnutsSampler` after each draw, is nearly but not
quite available without touching the vendored code. Of the five sampling tuning
values, three (`max_trajectory_doublings`, `max_step_halvings`,
`max_hamiltonian_error`) come from this package's own `SamplingConfig` and two
have getters (`inverse_mass_matrix_diagonal`, `macro_time`). The sixth,
`min_micro_steps`, comes from the adapter's own estimator and `WalnutsSampler`
exposes no getter for it, so a rebuild silently loses it. That path also reaches
only the sampling phase.

On top of either, the density, the slice sampler, and their gate. Call it 150-250
lines and a design note, and it moves every draw.

**The centered alternative** - sampling the effects directly and the scale from
its own conditional - is the textbook answer for well-identified groups and is a
rewrite of `make_theta_L`, `make_b`, the hand adjoint, the raw block layout, the
gradient-gate fixtures, and every baseline. It also reintroduces the funnel for
weakly identified designs, so it wants a switch rather than a swap.

**The shared-level defect** has its own section above; what belongs here is the
warmup half of it. `y ~ bart(...) + (1 | g)` has no fixed-effect column at all -
the design matrix is empty, K = 0 - so nothing but the forest and the random
intercepts can carry the response mean, and early in warmup the forest has not
grown. The intercepts take it and give it back slowly. It is visible at
`skip = 1` too (the shared level is 1.96 at the first stored warmup draw and 0.2
by the end); `skip = 8` in warmup makes it eight times worse, which is why the
extra transitions are confined to sampling rather than fixed at the root. Fixing
it at the root - seeding the forest with the response mean, or carrying an
explicit intercept - would also make a uniform skip loop safe.

## Harness

benchmarks/R/bench-re-scale.R records the table above over three grouped designs
(20-group and 5-group Gaussian random intercepts, 20-group probit), a skip grid,
and five seeds, and prints the verdict against the bar. It also records the
shared level's own autocorrelation and ESS, so the second defect is visible in
the baseline rather than only here. It replaces the `rbart_vi` comparison the bar
came from, which cannot survive that function's removal from dbarts; the same
removal took the third comparator out of inst/tinytest/test-02-binary.R, whose
surviving glmer and plain-forest arms carry the intent.
