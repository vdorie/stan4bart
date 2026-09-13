# The random-effect scale: mixing and cost

Two complaints arrived together about grouped fits, and the first job of this
design is to say that they are two problems and not one.

The first is mixing. dbarts is retiring `rbart_vi` in favour of this package and
set a release bar for the hand-off: on a Gaussian random intercept, n = 2000,
K = 20 equal groups, Friedman f over 10 predictors, tau = 1, sigma = 1, one
chain, 1000 warmup and 1000 kept draws, 200 trees, the random-effect standard
deviation must reach lag-1 autocorrelation below 0.8 and ESS at least 100 per
1000 kept draws. It reaches 0.95-0.97 and ESS 2-28.

The second is wall time. bartCause's `group.by` route now runs through this
package, and its side-by-side against the old dbarts route is several times
slower.

The measurements below say the two have almost nothing to do with each other.
The mixing failure is a geometry problem in the parametric block that costs
almost no time; the wall-time gap is four unrelated costs, of which the
parametric block is the smallest. A remedy aimed at one can make the other
worse, and the remedy that landed first - more parametric transitions per sweep,
`skip`'s `"stan"` element - does exactly that. This design recommends the one
move that fixes the mixing without paying wall time for it, and lists the
wall-time items separately.

Measurement conditions: arm64 macOS, R 4.6.1, dbarts 8931adac, stan4bart
979a91a plus the working tree. Wall times were NOT taken on a quiet machine;
read the ratios, not the seconds. The mixing figures are exact and reproduce
bit-for-bit.

## 1. The mixing root cause

Fully measured in docs/plans/re-scale-mixing.md; the result in one paragraph.

The random-effect block is non-centered - `b_level = T_i * z_b_level` in
`ParametricModel::eval` - and the scale is already sampled on the log scale
(`tau[i] = exp(tau_f_p[i])`). Neither spelling is the fault. The fault is that
every block's Cholesky factor is homogeneous of degree one in its own scale
(`bl.s = tau[i] * re_scale[i] * dispersion`, and every entry of `T` is a
multiple of `bl.s`), so scaling `tau_i` by r and dividing that block's `z_b` by
r leaves `b`, the linear predictor, and the likelihood exactly untouched. With
100 observations per group the effects are pinned by the likelihood, so the
posterior concentrates on that curve: `cor(log tau, log sum z_b^2) = -0.978`.

In the coordinates the sampler integrates - the unconstrained position divided
by the square root of its adapted diagonal inverse mass - the block has one
direction of length 6.5 carrying 90 percent of the scale's variance, one of 4.9
that is the shared level of the standardized effects, and twenty of length 0.4
to 1.1. The long directions are rotated, spread over the scale and all twenty
standardized effects, so a diagonal metric cannot rescale them. Their lag-1
autocorrelations are 0.972 and 0.968; the median over the other twenty is 0.107.

The sampler stops short of the ridge because the NUTS U-turn test is dominated
by the twenty short directions: 10.5 leapfrog steps at step size 0.364 is an
integration time of 3.8, against a length scale of 6.5. **The trajectories are
too short, not too long.** That is the fact that separates this from the wall
time, and section 3 turns it into a measurement.

Ruled out with measurements, each in the plan file: adaptation (warmup 4000
leaves lag-1 at 0.933 and 0.978 on two seeds, and the metric cannot converge on
a direction the chain never crosses), the NUTS doubling cap (raising it from 5
to 12 reproduces the draws bit for bit), `adapt_delta` (0.99 buys 0.896), the
prior (widening `decov`'s scale to 10 buys 0.943), and the two-block alternation
(freezing the forest leaves the scale at 0.954 while taking the shared level
from 0.980 to 0.093 and sigma's ESS from 13 to 633).

## 2. The wall-time root causes

The bartCause configuration: n = 1000, K = 20, four chains, 500 kept draws,
otherwise bartCause defaults.

| what ran | seconds | ratio |
| --- | --- | --- |
| dbarts route (`use.ranef = FALSE`), 500/500/4 | 1.62 | 1.0 |
| stan4bart route, as the user writes it | 28.94 | 17.9 |
| stan4bart route, configuration matched | 5.96 | 3.7 |
| stan4bart route, matched + `cores = 4` | 3.98 | 2.5 |

Four separate costs, largest first.

**(a) bartCause drops the sample and chain arguments on this route, 5x.** The
call it builds is
`stan4bart(..., verbose = -1L, chains = 10L)`: `n.samples`, `n.burn` and
`n.chains` never reach the fit. It runs 10 chains at this package's own
`iter = 2000, warmup = 1000` - 20000 sweeps where the user asked for 4000. The
dbarts route honours all three. This is a bartCause defect, not one of ours, and
it is the single largest factor whenever it bites; it is also why the observed
gap is bigger than the engine gap. Reported upward, not fixed here.

**(b) Chain parallelism, up to 2x and negative below the crossover.** dbarts
threads chains inside its engine with pthreads at no startup cost
(`n.threads = guessNumCores()` by default, 0.84 s serial to 0.51 s threaded
here). This package builds a PSOCK or FORK cluster in
`stan4bart_fit.R`, and `cores` defaults to `getOption("mc.cores", 1L)`, so four
chains run serially unless the caller asks otherwise. Measured cluster startup
is 1.25 s flat. At the bartCause fit size `cores = 4` is a net **loss** on the
bare fit (1.52 s against 1.25 s at iter 1000); the crossover is well above it
(4.38 s against 9.12 s at iter 8000). bartCause never sets `cores`.

**(c) The counterfactual test surface, 3x on the fit.** bartCause hands the
flipped treatment frame in as `test`, so the fit carries 2000 rows rather than
1000 and evaluates the forest on the test half every stored draw: the bare fit
in bartCause's shape goes from 1.21 s to 3.78 s. Both routes pay a counterfactual
surface, so only the difference is ours, but it is where most of the matched
5.96 s sits.

**(d) The parametric block itself, about 1.4x per sweep - the smallest.** By
differencing the skip grid on one chain: a WALNUTS transition costs 0.079 ms and
a BART sweep 0.19 ms, so the parametric block is 29 percent of per-sweep cost.
dbarts draws the same scale by conjugate Gibbs at O(K). Deleting our parametric
block outright would buy 1.4x and no more.

The item to take from this table is that **`skip` is the wrong instrument.** At
0.079 ms per transition, `skip = c(bart = 1, stan = 16)` adds 1.2 ms per sweep
against a 0.19 ms sweep - it multiplies the fit's cost by about 2, on the axis
that is already the complaint, to fix a problem that costs nothing to have.

## 3. They are not the same phenomenon

The question worth pre-empting is whether the ridge is also what drives the
gradient-evaluation count, in which case one fix would serve both. It is not.

| model | unconstrained dim | leapfrog per transition |
| --- | --- | --- |
| `y ~ bart(...) + x1` (no random effects) | 2 | 8.6 |
| `y ~ bart(...) + (1 \| g)` | 22 | 10.4 |

Eleven times the dimension, and the trajectory grows by a fifth. The
random-effect block is not inflating the step count, and it cannot be: the
U-turn fires on the twenty short directions, which is precisely the mechanism
that leaves the ridge uncrossed. If the ridge were inflating trajectories the
sign would be the other way and the scale would mix.

The same reading holds across designs. In the baseline, `gaussian_k5` - 400
observations per group, so the longest ridge of the three and the worst-mixing
scale (lag-1 0.99) - has the highest leapfrog count at 13.8-15.8 against 10.5 for
`gaussian_k20`. The ridge does lengthen trajectories a little, but a factor of
1.4 on a term that is 29 percent of the sweep, while the autocorrelation it
tracks moves from 0.97 to 0.99. One problem is a factor of a few in time; the
other is two orders of magnitude in ESS.

Fixing the scale at its posterior mean was the proposed direct test. It needs a
worktree hack - there is no interface for it - and the two measurements above
already answer the question in the direction the hack would, so it was not
built. Recorded as untested.

## 4. The remedy

This section is the proposal as it was written. It has landed; section 6 is the
landing note, and it corrects two things stated here.

The ridge is exactly traversable in closed form, and the conditional along it is
log-concave. That is the whole of the recommendation.

Move along the homogeneity curve: `tau_i -> tau_i * e^u` with that block's `z_b`
divided by `e^u`. `b`, the linear predictor and the likelihood are invariant, so
the target for `u` carries no likelihood term at all - only the standardized
effects' normal prior, the scale's own prior, and the change-of-variables
Jacobian. Writing `t` for the new scale, `q` for the block's coordinate count
and `B = ||b_block||^2 / (re_scale * dispersion)^2` (invariant along the curve),
and reading the prior off `ParametricModel::eval` (`tau ~ Gamma(shape, 1)` with
its log Jacobian):

```
log pi(t) = (shape - q) log t - t - B / (2 t^2),    t > 0
```

Each term is concave in `t` for `shape <= q`, which holds for every reachable
model (`shape` defaults to 1, `q` is the block size). So the conditional is
log-concave and one-dimensional: a stepping-out slice sampler draws it exactly,
with no tuning and no gradient, at O(q) to form `B` plus a handful of scalar
evaluations. Against a 0.19 ms BART sweep this is free. It is also, in spirit,
what `rbart_vi` does - which is why its lag-1 autocorrelation on this design is
0.06.

### The alternatives, and why not

**A dense or per-coordinate mass matrix.** A dense metric would rotate the ridge
into a coordinate and fix it. It costs an O(d^2) estimator and an O(d^2)
solve per leapfrog step in a block whose dimension grows with the number of
group levels (q = 20 here, thousands in a large fit), and the warmup cannot
estimate the very direction at issue - the metric is already wrong on the scale
coordinate by a factor of 3.8 in sd and lengthening warmup does not converge it.
Treating scale coordinates specially with a hand-set mass is the same guess in
cheaper clothing: it fixes a length, not a rotation, and the ridge is rotated.

**A different parameterization of the scale coordinate.** Log scale is what is
already sampled; sampling `tau^2`, or `log tau` with a different Jacobian, moves
the ridge and does not remove it, because the ridge is a property of the
`tau * z_b` product and not of the chart the scale is read in. The prior sweep
is the empirical version of this argument: four priors, lag-1 0.94 to 0.98.

**Centering the random-effect block** - sampling the effects directly and the
scale from its own conditional - is the textbook answer for well-identified
groups. It reintroduces the funnel for weakly identified ones, so it wants a
switch rather than a swap, and it is a rewrite of `make_theta_L`, `make_b`, the
hand adjoint, the raw block layout, the gradient-gate fixtures and every
baseline.

**`skip`'s `"stan"` element**, which has landed, reaches the bar at 8 to 16
transitions per sweep. Section 2 is the argument against making it the default:
it buys the mixing on the wall-time axis, at about 2x. It stays as the
user-facing escape hatch and as the thing that makes the bar reachable today.

**Chain-level parallelism** belongs on the wall-time list, not this one. It
does not touch mixing.

### What it costs

The density and the slice sampler are perhaps 60 lines. WALNUTS holds its
position privately - `WalnutsSampler::theta_` has no setter - but **the move
needs nothing added to the vendored headers.** Rebuild the sampler instead.
`WalnutsSampler` is copyable and its constructor takes the position and every
tuning value: the base generator, the handler, the log density, `theta`,
`inv_mass`, `macro_time`, `max_nuts_depth`, `max_step_halvings`,
`min_micro_steps` and `max_error`. Reading the position back is already solved,
since `LatestDraw::on_sample` is handed `theta_` after every transition, so
constructing a fresh sampler at the slice-updated position is the setter.

Every value the rebuild needs is reachable. `inv_mass`, `macro_time` and
`max_error` have getters on the sampler; `max_nuts_depth` and
`max_step_halvings` come from this package's own `SamplingConfig`; and
`min_micro_steps`, which used to be the gap in this route, is a public getter on
`AdaptiveWalnuts`. The base generator is held by reference throughout, so a
rebuilt sampler draws from where the old one left off rather than restarting it.
It is not bit-identical to a persistent sampler: the vendored `Random` wrapper
holds its own `std::normal_distribution`, which caches a spare variate, and a
rebuild constructs a fresh one. The stream stays valid; it stops being the same
sequence, so the move has to be gated on distributions rather than on draws.

Three things the route costs, none of them a blocker:

- `WalnutsSampler::freeze` currently drops the adapter, so `min_micro_steps`
  has to be read off it and kept before the adapter goes. That is a
  package-side change of one line and one member.
- The rebuilt constructor evaluates the log density and gradient at the handed
  position. That evaluation is already paid once a sweep, because the vendored
  sampler caches those values across transitions and this package has to
  refresh them after swapping the target; a sweep that rebuilds does not need
  to refresh as well. That holds in the sampling phase only - warmup still
  needs the refresh, for the reason in the next item.
- The route reaches the sampling phase only. Warmup runs `AdaptiveWalnuts`,
  which has no constructor from a running state, so rebuilding it would reset
  Adam and the mass estimator. The move applies after the freeze, which is what
  the bar is stated on.

With gates and baselines, call it 150-250 lines and a full re-record: it moves
every draw of every model with a random-effect block. That is not small and not
local, so it is a scheduled item and not this pass's edit.

### The acceptance bar, pre-registered

On the bar's design (Gaussian, n = 2000, K = 20, tau = 1, one chain, 1000 warmup
and 1000 kept, 200 trees), at the DEFAULT `skip`, worst of five seeds:

- scale lag-1 autocorrelation below 0.8
- scale ESS at least 100 per 1000 kept draws
- scale posterior mean within Monte Carlo error of the value the current
  sampler converges to, so that mixing is not bought with bias

On the bartCause configuration (n = 1000, K = 20, four chains, 500 kept, the
counterfactual test surface, `cores = 1`):

- total wall time within **1.1x** of the same fit with the move disabled

The wall-time bar is stated against this package's own fit and not against the
dbarts route's 1.62 s, because the four costs in section 2 are independent of
this remedy and three of them are not ours to move: a remedy that is free should
be gated on being free, not on carrying the other three. Closing the gap to the
dbarts route is the separate list below. If the slice move cannot hold 1.1x it
is not the cheap fix this design claims and the recommendation fails.

Secondary, not gated: `gaussian_k5` (400 observations per group, the longest
ridge) should also clear 0.8, since the remedy is exact along the curve and does
not degrade as the ridge lengthens. If it does not, the diagnosis is incomplete.

## 5. The wall-time list, separately

Not part of the remedy above; recorded so the two are not conflated again.

1. bartCause forwards neither `n.samples` nor `n.burn` nor `n.chains` on the
   stan4bart route, and never sets `cores`. Largest factor, and it is theirs.
2. Chain parallelism costs 1.25 s of cluster startup and is a net loss below
   roughly 2000 iterations x 4 chains. In-engine threading, the way dbarts does
   it, would remove both the startup and the serial default; it is a real piece
   of work because the parametric block reaches R.
   Short of that, `cores` could pick itself from the problem size rather than
   defaulting to 1, which is small and would take the matched bartCause case
   from 5.96 s to 3.98 s.
3. The test surface triples the fit. Worth confirming that the test rows are not
   being evaluated more often than the stored draws need.
4. The parametric block at 29 percent of per-sweep cost is the floor, and it is
   the price of the model this package fits.

## 6. What landed

The slice move, on by default, drawn once per sweep for every random-effect
block during sampling. Warmup is unchanged.

### The derivation, as implemented

The move itself is section 4's: `tau_i -> tau_i * e^u` with that block's `z_b`
divided by `e^u`. That leaves every entry of the block's Cholesky factor scaled
by `e^u` and every standardized effect divided by it, so `b`, the linear
predictor and the likelihood are exactly invariant - for a scalar block, for a
correlated block under `decov`, for several random-effect terms at once (each
has its own `tau` and its own `z_b` segment, so the moves are separate
conditionals), and for the binary family, whose latents enter only as the
response of that same invariant likelihood.

The conditional was re-derived from `ParametricModel::eval` rather than taken
from section 4, and it differs from what section 4 states. Write `x = log tau_i`
- which IS the unconstrained coordinate the sampler carries - and `v` for the
block's standardized effects in `R^{q_i}`, `q_i = p_i * l_i` counting every
coordinate of every level of the block. Change variables `(x, v) -> (x, w)` with
`w = e^x v`, the quantity the curve holds fixed. The Jacobian of `v -> w` at
fixed `x` is `e^{q_i x}` times the identity, so the density in the new chart
carries a factor `e^{-q_i x}`. The only terms of the log density that involve
`x` or `v` are the standardized effects' `N(0, 1)` prior, `tau`'s
`Gamma(shape_i, 1)` prior, and the log Jacobian of `tau = e^x`. Collecting them,

```
log p(x | w, rest) = (shape_i - q_i) x - e^x - (A_i / 2) e^{-2x},
A_i = ||w||^2 = tau_i^2 ||z_b(block i)||^2
```

with `(shape_i - 1) x - e^x` the Gamma prior, `+x` its Jacobian, `-q_i x` the
chart's, and the last term the normal prior at `v = e^{-x} w`. The second
derivative is `-e^x - 2 A_i e^{-2x} < 0`, so the conditional is log-concave on
all of `R` unconditionally - there is no `shape <= q` side condition, which
section 4 needed only because it wrote the density in `t`. `re_scale` and
`dispersion` cancel: they multiply `tau_i` and `z_b` in one product.

Two corrections to section 4:

- Its `log pi(t) = (shape - q) log t - t - B / (2 t^2)` is a density in `t`
  carrying the exponent of the density in `log t`. In `t` the exponent is
  `shape - q - 1`; the Jacobian `dt/dx` is missing. What is implemented is the
  `x` form above, where `shape_i - q_i` is right.
- Its invariant `B = ||b_block||^2 / (re_scale * dispersion)^2` equals `A_i`
  only for a scalar block. For a correlated block `T` is not a multiple of an
  orthogonal matrix, so `||b||^2` is not `s^2 ||z_b||^2`. `A_i` is also the
  cheaper of the two: it needs no `b`.

A stepping-out slice sampler draws it (Neal 2003, interval stepped out with a
split budget and then shrunk), at `O(q_i)` to form `A_i` plus a handful of
scalar evaluations and no gradient. The interval's width is a function of the
block's geometry only - a width read off the current position would cost the
draw its reversibility.

The invariance the derivation rests on is checked in `ridgeMove` itself, under
`NDEBUG`, by comparing the linear predictor across the move: a shipped build
pays nothing for it, and a build configured with `-UNDEBUG` runs it on every
move of every fit. Over five hundred moves of a two-block Gaussian fit and a
two-block binary one the two linear predictors agree to 8.9e-15 absolute, and
the check's band catches a relative error in the scaling above about 1e-9.

### Where it sits in the sweep

At the top of the parametric block: after the BART draw has set the offset and,
for a binary response, the latent response, and before the WALNUTS transitions.
The conditional carries no likelihood term, so the move is stationary for the
parametric conditional wherever it is applied; that position is chosen because
it is where the rebuild's own log density and gradient evaluation is the one the
stale-target cache already owed. A sweep that moves therefore pays no
evaluation it was not already paying, and the mean leapfrog count per transition
is unmoved (10.65 to 10.84 on the bar's design).

Section 4's route is what was built: nothing was added to the vendored headers.
The sampler is rebuilt at the moved position from the tuning captured at the
freeze - `min_micro_steps` read off `AdaptiveWalnuts` before the adapter is
dropped, `inv_mass`, `macro_time` and `max_error` off the frozen sampler, the
two caps off this package's own `SamplingConfig`. The base generator is held by
reference and continues; the rebuilt sampler's own `normal_distribution` is
fresh and caches a spare variate, so the move is gated on distributions and not
on draws. `LICENSE.note`'s vendored-verbatim claim is untouched, and the
vendored tree still carries exactly one local patch.

### The bar

Bar design (Gaussian, n = 2000, K = 20, tau = 1, one chain, 1000 warmup and 1000
kept, 200 trees), at the DEFAULT `skip`, the five harness seeds, the move off
and on:

| gate | move off | move on | bar |
| --- | --- | --- | --- |
| scale lag-1, worst seed | 0.963 | 0.090 | below 0.8 |
| scale lag-1, median | 0.957 | 0.059 | - |
| scale ESS per 1000, worst seed | 7.0 | 526 | at least 100 |
| scale ESS per 1000, median | 13.5 | 739 | - |
| scale posterior mean, worst seed shift | - | 0.78 combined MCSE | within MCSE |
| seconds | 2.15 | 2.08 | - |

The seconds row was not taken on a quiet machine and does not reproduce to the
figure: read the two columns against each other, not the absolute times. Every
other row is the harness's `ridge` mode and reproduces exactly. PASS on all
three clauses. Per seed the two arms' posterior means for the scale
differ by 0.01 to 0.78 of two combined Monte Carlo errors, and the move cuts
that Monte Carlo error by four to eight times, seed for seed (0.029-0.050 down
to 0.005-0.007).

Wall time on the bartCause configuration (n = 1000, K = 20, four chains, 500
kept, the counterfactual test surface, `cores = 1`), eight paired repetitions:
median 6.36 s with the move off against 6.43 s with it on, a ratio of **1.01**;
means 6.51 s and 6.54 s, 1.004. The bar is 1.1x. PASS. The machine was not
quiet, which is why the comparison is paired and read on the median.

Secondary, not gated: `gaussian_k5` (400 observations per group, the longest
ridge) goes from worst-seed lag-1 0.986 to 0.150, so it clears 0.8 on every
seed, and the diagnosis is not incomplete. Its ESS half does not follow on one
of the five seeds, which reads 10.5 against a lag-1 of 0.150 - the
effective-sample-size estimator at 1000 draws is the noisier of the two, and
that seed's chain carries structure past lag one that the scale's own
conditional does not produce. The third design, `probit_k20`, goes from
worst-seed 0.947 / 5.9 to 0.248 / 232 and clears both halves.

Bias is the group-sd harness's `bias` mode - the four grouped designs at three
seeds each, the same data and the same MCMC seed in both arms, at the package
defaults. All 36 posterior-mean differences (group sd, residual sd, fixed
effect) fall within two combined Monte Carlo errors, worst ratio 0.94, and 16 of
the 36 differences are positive, so there is no systematic direction. A second
check on a design this harness does not cover - a binary response with a random
intercept block and a correlated intercept-and-slope block, three seeds - puts
all 15 of its comparisons inside the same band, worst 0.79.

The full tinytest suite passes (551 expectations, up from 546 by the five the
move's own test adds), the posterior baseline gate passes on all five tiers, and
the tree-replay exactness gate passes on all three. No test pins draws of a
random-effect model, so nothing needed the off switch or a regenerated value.

### The current state

The maintainer's ruling of 2026-09-13 is what was built: the move rebuilds the
frozen `WalnutsSampler` each sweep rather than editing the vendored headers. The
second local patch the plan held in reserve - a position setter on
`AdaptiveWalnuts`, so the move could run in warmup too - was **not needed**: the
sampling-only move clears the bar at the default `skip` on every seed, so warmup
adaptation is not what holds the scale back and the vendored tree keeps its one
patch.

The move is **on by default**. `stan_args = list(ridge_move = FALSE)` turns it
off, which is how the two columns above were measured and is the only reason to
reach for it. `skip`'s `"stan"` element stays where it was, default 1, as the
escape hatch for anything the move does not reach.
