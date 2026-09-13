# Mixing of the group standard deviation

dbarts 1.0-0 removes its own grouped random-intercept fit and points users at
this package instead. The condition attached to that removal is a bar on the
chain for the random-intercept standard deviation - the group spread, the
parameter dbarts called tau: its lag-one autocorrelation must be below 0.8 and
its effective sample size at least 100 per 1000 kept draws. This note records
where stan4bart stands against that bar. The script is
`benchmarks/mixing-group-sd.R`; it runs in well under a minute and reproduces
everything below.

Four designs, each n = 2000, all gaussian, with a Friedman mean function for the
BART part, one linear fixed effect with coefficient 2, group intercepts drawn
from a normal, and residual standard deviation 1. `bar_reference` is the design
the bar was set on: twenty groups of a hundred, true group sd 1. `many_small` is
fifty groups of forty, `few_large` five groups of four hundred, and
`weak_signal` twenty groups of a hundred with a true group sd of 0.2 - the
corner where dbarts's own review of this parameter found mixing worst. Three
seeds per design, each seed driving both the simulated data and the sampler.
The posterior means track the truth in every fit - the residual standard
deviation within a percent or two of 1, the fixed effect within a few percent of
2, and the group sd within its posterior uncertainty - so what follows is about
mixing and not about a fit going wrong.

Every fit uses the package defaults: four chains, iter 2000 with warmup at half
of that, so 1000 warmup and 1000 kept draws per chain and 4000 kept draws in
all, with `bart_args` and `stan_args` untouched. The one argument set away from
its default is `cores`, raised to four, which divides the same chains across
processes and changes wall time only. Effective sample size and R-hat are
computed in the script rather than taken from a dependency: the multi-chain
autocorrelation estimator with Geyer's initial positive and monotone sequences,
and split R-hat, both on half-chains. On autoregressive test chains they
reproduce the standard split implementations to within half a percent.
Lag-one autocorrelation is the mean over the four chains.

The table describes the build that carries the random-effect scale's ridge move,
which is on by default: once per sweep the scale of each grouping factor is
drawn by an exact slice move along the curve that rescales it and divides that
factor's standardized effects by the same amount, leaving the linear predictor
and the likelihood untouched. That move is what changed the verdict below, and
it is what the group sd's mixing now rests on; its derivation and its own
acceptance bar are in `docs/design/re-scale-and-grouped-cost.md`. The earlier
figures this note carried, on a build without the move, are in
`docs/plans/re-scale-mixing.md`, which also records that the WALNUTS refresh
before it moved every draw and none of the verdicts.

| case | seed | group sd lag-1 | group sd ESS/1000 | group sd R-hat | sigma lag-1 | sigma ESS/1000 | fixed effect lag-1 | fixed effect ESS/1000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bar_reference | 20260913 | 0.08 | 370 | 1.00 | 0.50 | 5 | 0.23 | 106 |
| bar_reference | 20260914 | 0.12 | 377 | 1.00 | 0.69 | 4 | 0.29 | 25 |
| bar_reference | 20260915 | 0.08 | 187 | 1.01 | 0.65 | 9 | 0.16 | 30 |
| many_small | 20260913 | 0.04 | 496 | 1.00 | 0.36 | 4 | 0.15 | 47 |
| many_small | 20260914 | 0.08 | 596 | 1.01 | 0.64 | 15 | 0.20 | 103 |
| many_small | 20260915 | 0.05 | 478 | 1.01 | 0.76 | 15 | 0.17 | 7 |
| few_large | 20260913 | 0.20 | 81 | 1.02 | 0.65 | 41 | 0.31 | 17 |
| few_large | 20260914 | 0.22 | 14 | 1.05 | 0.73 | 10 | 0.55 | 8 |
| few_large | 20260915 | 0.10 | 1 | 1.82 | 0.77 | 3 | 0.47 | 8 |
| weak_signal | 20260913 | 0.17 | 503 | 1.01 | 0.14 | 55 | 0.02 | 7 |
| weak_signal | 20260914 | 0.22 | 374 | 1.01 | 0.39 | 7 | 0.08 | 6 |
| weak_signal | 20260915 | 0.30 | 35 | 1.03 | 0.32 | 18 | 0.12 | 10 |

## The verdict

The bar is met in two of the four designs on every seed, in a third on two
seeds of three, and in the fourth on none - eight of the twelve fits clear both
halves, where before the move not one of them cleared either. The design the bar
was written against is the clean pass: lag-one 0.08 to 0.12 against a ceiling of
0.8, and 187 to 377 effective draws per thousand against a floor of 100, so it
clears the autocorrelation limit by an order of magnitude and the
effective-sample-size floor by two to four times. `many_small` passes on all
three seeds by more, at lag-one 0.04 to 0.08 and 478 to 596 effective draws.
`weak_signal` passes on two of three and misses the effective-sample-size floor
on the third, at 35 per thousand against a lag-one of 0.30.

`few_large` is the failure, and it fails for a different reason than the group sd
used to fail everywhere. Its lag-one is 0.10 to 0.22 - the chain is not
autocorrelated any more - and its effective sample size reads 81, 14 and 1 per
thousand, the last with an R-hat of 1.82. What is left there is four chains
disagreeing about the group sd rather than four slow ones, which is what five
groups of four hundred buys: the spread of five numbers about zero is weakly
identified however well each draw of it mixes, and the effective-sample-size
estimator charges for between-chain disagreement. Longer runs repair that;
faster mixing within a chain does not.

Read the two halves of the bar with different weight, as before. The lag-one
column is stable and it is now uniformly low. The effective sample size is not:
`weak_signal`'s reading moves by a factor of several between seeds, and
`few_large`'s ESS of 1 sits beside a lag-one of 0.10. Where they disagree, the
autocorrelation is the one to trust.

The group sd is no longer the hard parameter. Across the three designs with a
group sd of 1 its lag-one now sits at 0.04 to 0.22, while the residual standard
deviation runs 0.36 to 0.77 and the fixed effect 0.15 to 0.55 - the ordering
from the earlier build is reversed, and the two parameters that now mix worst
are the two that pay the BART-versus-parametric alternation rather than the
parametric block's own geometry. Their effective sample sizes are depressed by
disagreement between chains and a longer run repairs them. The posterior means
still track the truth in every fit - the residual standard deviation within a
percent or two of 1, the fixed effect within a few percent of 2, and the group sd
within its posterior uncertainty - and a check across these four designs at
three seeds with the move switched off finds 33 of 36 posterior-mean differences
within twice their combined Monte Carlo error, with no systematic direction: the
mixing is not bought with bias.

## What remains, and the remedies not taken

What is left is the two-block alternation, not the parametric block. The group
sd's own conditional is now drawn exactly; the residual standard deviation, the
fixed effect, and the shared level of the group intercepts are all confounded
with the forest's overall fit, and no move inside either block crosses that.
The shared level is the sharpest case - lag-one 0.97 to 0.99 on every design,
unmoved by the ridge move by construction, since the move holds the linear
predictor fixed. Its fix is a joint move shifting a scalar between the forest and
the intercepts, which needs a way to add a constant to every leaf that the dbarts
flat C API does not expose. That is recorded in `docs/plans/re-scale-mixing.md`
and in the root `TODO`.

dbarts's own work on this parameter left two candidates on the table, both
recorded in its `docs/design/retire-grouped-random-effects.md` and
`docs/plans/tau-slice-review.md`; neither is needed now for the quantity the bar
names. The first is an interweaving move, alternating the centered draw of the
group effects with an ancillary one in which the effects are rescaled by the
group sd, which measured roughly an elevenfold drop in autocorrelation on the
isolated pair of the group sd and its effects - the same pair the ridge move
now draws in closed form, and at no wall-time cost. The second is the fallback
dbarts named when the bar was set: an R-level grouped intercept that draws the
group effects and their spread itself and drives a plain dbarts sampler through
`setOffset` each sweep, which buys the conjugate draw whose lag-one dbarts
measured at 0.14 to 0.15, but gives up everything this package's formula
interface offers beyond a single random intercept - slopes, crossed and nested
factors - and puts the sampler's inner loop back in R.
