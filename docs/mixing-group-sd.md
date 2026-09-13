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

| case | seed | group sd lag-1 | group sd ESS/1000 | group sd R-hat | sigma lag-1 | sigma ESS/1000 | fixed effect lag-1 | fixed effect ESS/1000 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bar_reference | 20260913 | 0.95 | 18 | 1.08 | 0.43 | 10 | 0.24 | 10 |
| bar_reference | 20260914 | 0.96 | 19 | 1.02 | 0.66 | 6 | 0.33 | 9 |
| bar_reference | 20260915 | 0.93 | 12 | 1.07 | 0.54 | 9 | 0.17 | 29 |
| many_small | 20260913 | 0.88 | 39 | 1.01 | 0.40 | 21 | 0.17 | 107 |
| many_small | 20260914 | 0.91 | 26 | 1.04 | 0.60 | 15 | 0.23 | 12 |
| many_small | 20260915 | 0.90 | 24 | 1.04 | 0.76 | 6 | 0.11 | 10 |
| few_large | 20260913 | 0.96 | 9 | 1.11 | 0.64 | 3 | 0.31 | 11 |
| few_large | 20260914 | 0.96 | 12 | 1.05 | 0.67 | 32 | 0.54 | 12 |
| few_large | 20260915 | 0.96 | 1 | 3.43 | 0.77 | 3 | 0.39 | 4 |
| weak_signal | 20260913 | 0.62 | 164 | 1.01 | 0.07 | 10 | -0.03 | 13 |
| weak_signal | 20260914 | 0.73 | 70 | 1.03 | 0.39 | 8 | 0.16 | 33 |
| weak_signal | 20260915 | 0.65 | 18 | 1.05 | 0.38 | 4 | 0.15 | 18 |

## The verdict

The bar is met in none of the four cases. Exactly one of the twelve fits clears
it - the weak-signal design at the first seed, at lag-one 0.62 and 164 effective
draws per thousand - and every other fit misses at least one half of the bar,
most of them both halves. The reference design the bar was written against is
the plainest failure: lag-one 0.93 to 0.96 against a ceiling of 0.8, and 12 to
19 effective draws per thousand against a floor of 100, so it misses the
autocorrelation limit by about 0.15 and the effective-sample-size floor by a
factor of five to eight. `few_large` is worse - lag-one 0.96 on all three seeds
and effective draws of 9, 12 and 1 per thousand, the last with an R-hat of 3.4,
which is four chains that have not found the same posterior rather than four
slow ones. `many_small` is the least bad of the three designs with a group sd of
1 and still fails everywhere, at lag-one 0.88 to 0.91 and 24 to 39 effective
draws, missing the effective-sample-size floor by between two and a half and
four times. `weak_signal` is the only design that comes near: its lag-one
clears on all three seeds, but its effective sample size swings over an order
of magnitude across them, 164 then 70 then 18, so it passes once and misses by
factors of 1.4 and 5.6 on the other two. Group count and group size both matter
and point the same way - the fewer groups there are, the worse the group sd
mixes - but no design with a group sd large enough to be worth estimating
reaches the bar on any seed.

The group sd is the hard parameter, and it is hard in a way the other two are
not. Across the three designs with a group sd of 1 its lag-one autocorrelation
sits at 0.88 to 0.96, while the residual standard deviation runs 0.40 to 0.77
and the fixed effect 0.11 to 0.54. The effective sample sizes for those two are
also low, sometimes lower than the group sd's, but for a different reason: they
are depressed by disagreement between chains, and a longer run repairs it.
Quadrupling warmup and kept draws on the reference design lifts the residual
standard deviation from 10 effective draws per thousand to 15 and drops its
R-hat from 1.07 to 1.02, while the group sd's lag-one stays at 0.96 and its
effective sample size per thousand does not improve. That is the distinction
that matters for the bar: the other parameters are under-warmed at the defaults,
and the group sd is autocorrelated from one draw to the next, which no amount of
extra sampling fixes.

## Candidate remedies, not implemented

dbarts's own work on this parameter leaves two candidates on the table, both
recorded in its `docs/design/retire-grouped-random-effects.md` and
`docs/plans/tau-slice-review.md`. The first is an interweaving move, alternating
the centered draw of the group effects with an ancillary one in which the
effects are rescaled by the group sd, which measured roughly an elevenfold drop
in autocorrelation on the isolated pair of the group sd and its effects, but
costs a new sampling block, changes every draw this package makes, and addresses
only part of the problem, since dbarts attributed most of its own weak-signal
difficulty to the forest and the group effects competing to explain the same
group-level structure rather than to that pair. The second is the fallback
dbarts named when the bar was set: an R-level grouped intercept that draws the
group effects and their spread itself and drives a plain dbarts sampler through
`setOffset` each sweep, which buys the conjugate draw whose lag-one dbarts
measured at 0.14 to 0.15, but gives up everything this package's formula
interface offers beyond a single random intercept - slopes, crossed and nested
factors - and puts the sampler's inner loop back in R.
