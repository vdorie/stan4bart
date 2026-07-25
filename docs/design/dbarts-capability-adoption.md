# Adopting the bartcore capability set: design

Status: SURVEY, 2026-07-25. Not an arc; a ranked map of what the dbarts bartcore
program made reachable from this package, what each would cost, and which items
need dbarts to grow first. Individual items graduate to their own
`docs/plans/<item>.md` when picked up.

## 1. The seam

Every single-forest feature bartcore added is configured through attributes
parked on the control and model SEXPs plus a family string:
`bartcore.n.categories` (ordinal), `bartcore.dispersion` (nbinom),
`bartcore.survival` (aft), `bartcore.variance` (heteroscedastic),
`bartcore.hazard.periods`, `monotone`, `resid.df`, `interaction.max.order`,
`interaction.forbidden`, `block.of.column`, `block.tree.counts`. All three
SEXPs cross `dbarts_sampler_create(control, model, data, family)` untouched.

So most of the bartcore program is reachable from here with R-side work only.
Two things are not:

- MULTI-FOREST models (BCF, multinomial) are built by their own creation path
  inside dbarts, not exposed through `dbarts.h` at all.
- Any PER-OBSERVATION quantity coming back OUT of the sampler.
  `dbarts_results` carries a scalar `sigma` per draw and nothing per
  observation but the fits.

dbarts 1.0-0 exports `dbartsSpec()` (docs/design/consumer-spec-surface.md in
that repo), which resolves the control/model/data triple and family token the
way `dbarts()` does. Adopting it is the prerequisite for most of what follows:
today `R/stan4bart_fit.R` hand-builds a `dbartsModel` positionally and calls
`dbarts:::parsePriors` through a `quoteInNamespace` shim, which reaches the
priors and nothing else.

## 2. Ranked

1. ADOPT dbartsSpec. Replaces the `:::` shim, and makes every control/model
   attribute reachable in one change instead of one per feature. Also widens
   `bart_args`, which today forwards only names matching `dbartsControl`
   formals plus hand-wired k/power/base/split.probs.
2. ORDINAL responses. Structurally the existing binary path: `getLatents` ->
   the parametric block's response, with dbarts drawing the cutpoints
   internally. Multilevel ordinal regression is bread-and-butter for this
   package's audience and this is the cheapest genuinely new family.
3. `log_lik` / `loo` / `waic`. `dbarts_results.logLikelihood` is
   per-observation for gaussian, binary, and aft, and the BART sampler's
   offset IS the parametric mean, so what it reports is the joint conditional
   likelihood - already computed, currently discarded. This package has no loo
   support at all and its users arrive from rstanarm expecting it. Storage is
   the size of the BART train block, so it belongs behind the same opt-in
   discipline as `store=`.
4. PASS-THROUGHS to wire, document, and test: `storage = "single"`,
   `missing = "incorporate"` (MIA lets BART-only covariates carry NAs; the
   parametric side still needs complete data), monotone, sparse columns,
   warm starts, and `interactions()` / `blocks()`. Note the last one:
   constraining how group-level covariates enter the forest is a cheap R-side
   mitigation of the forest-ranef confounding dbarts investigated and declined
   to fix in the engine (its forest-ranef-interweaving.md, NO-GO).
5. PER-OBSERVATION VARIANCE CHANNEL - one dbarts change unlocking two families.
   Robust-t residuals and a heteroscedastic variance forest both need
   per-observation precision to flow BART -> parametric, reversing today's
   scalar `setSigma` direction (`init.cpp`, the Gibbs loop). dbarts would
   append a working-weights/variance field to `dbarts_results` (additive,
   MINOR bump, explicitly supported by the ABI); this side would set the
   parametric weights per iteration - the Eigen model already carries a
   `has_weights` vector, fixed at construction today. Do these as one arc.
6. AFT WITH FRAILTY. Cheap IF the imputed censored response is reachable
   through `getLatents` - worth a short check, because multilevel survival is
   an unoccupied niche. `dbartsSpec(family = "aft", survival = )` already
   builds the specification.
7. MULTILEVEL BCF. The largest research payoff and the one bartCause would
   consume, but it needs a genuinely new `dbarts.h` entry point for the
   two-forest sampler. Decision-gated on that.

Count families (nbinom) sit below all of these: the log link makes the
parametric block a weighted gaussian on a Polya-Gamma working response, which
is real WALNUTS-side work rather than wiring.

## 3. What not to do

- Do not chase within-chain threading or GPU offload for the BART block.
  dbarts measured both dead on the real engine (~1.1x, capped by a ~47% serial
  fraction); multi-chain remains the way to use cores.
- Do not make `storage = "single"` a default. It is opt-in, machine-dependent,
  and worth ~10-17% only at large n.
- Do not absorb the random-intercept-only case. dbarts now has grouped random
  effects in-engine (`rbart_vi`), so this package's documented niche is
  anything past a random intercept: random slopes, multiple grouping factors,
  covariance priors, fixed-effect priors. Worth stating in the README so users
  pick correctly, rather than competing for it.
