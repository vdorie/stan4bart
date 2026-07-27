# dbarts-spec-adoption

scope: build the BART component's control/model/data triple with
  dbarts 1.0-0's exported `dbartsSpec()` instead of hand-assembling a
  `dbartsModel` around a `dbarts:::parsePriors` call, and widen
  `bart_args` from "dbartsControl formals plus hand-wired
  k/power/base/split.probs" to the whole model-level half of a dbarts
  specification. Prerequisite for ordinal-outcomes, log-lik-loo,
  dbarts-passthroughs, and aft-frailty; the `:::` shim already broke
  once at the 1.0 port. Survey: docs/design/dbarts-capability-adoption.md
  section 2 item 1.

budget: S. One commit. Draw-neutral - bitwise, verified below.

## What was there

`stan4bart_fit.R` built the control by intersecting `bart_args` against
`dbartsControl`'s formals, then called `dbarts:::parsePriors` through
`quoteInNamespace`, splicing `normal(k = )` and `cgm(power=, base=,
split.probs=)` into the call by position (`which(as.character(prior_call)
== "normal")`), and finally constructed `new("dbartsModel", ...)`
positionally with a hand-chosen `node.scale` of 3.0 or 0.5. Nothing past
the priors was reachable, and `control@binary` was poked directly.

## What landed

`dbartsSpec(data, control, resid.prior = fixed(1), family = <token>,
parentEnv = evalEnv)`, with every `bart_args` name matching one of its
formals forwarded into the call. It resolves the family, the response
coding, `control@binary`, `node.scale`, the priors, and the monotone /
interaction / block attributes the same way `dbarts()` does, so a family
can never resolve two ways between the two entry points.

Four decisions worth recording:

- RESERVED NAMES. `sigma`, `resid.prior`, `resid.dist`, and `variance`
  are refused with an error rather than forwarded. The parametric block
  draws the residual sd and pushes it into the forest with `setSigma`
  every sweep; dbarts therefore holds sigma fixed at 1 and has no
  residual model of its own to configure. A variance forest or a
  student residual distribution would be silently wrong here rather
  than merely untested - both are the per-observation-variance arc,
  which is dbarts-gated on a working-weights field in `dbarts_results`.
  `family` is likewise refused unless it names the token the response
  already implies; the parametric block reads probit latents, so
  "logistic" and the rest are not wired into the Gibbs loop.
- PRIOR VOCABULARY. dbarts's prior constructors are not exported - they
  are a list `parsePriors` binds into its own evaluation environment -
  so a bare `cgm`/`normal`/`chisq` in `bart_args` would not resolve at
  the point `stan4bart` evaluates it. `stan4bart.R` already bound `chi`
  to a function returning its own `match.call()` for exactly this
  reason; that is now a loop over `names(dbarts::dbartsPriors)`, so any
  prior expression survives to dbartsSpec as a call and is resolved
  there. Consequence, unchanged from `chi`: a prior must be spelled
  INLINE in the `stan4bart` call, since `bart_args` is evaluated from
  the matched call.
- SHORTHANDS KEPT. `k`, `power`, `base`, `split.probs` still write into
  `node.prior` / `tree.prior`. Supplying both a shorthand and the prior
  it writes into is now an error instead of a silent win for the
  shorthand.
- FAMILY TOKEN TO C. `dbartsSpec` parks the resolved family on
  `model@family`, and that is what `dbarts_sampler_create`'s fourth
  argument wants; `init.cpp` read it off the model SEXP in place of the
  `""` shape-dispatch request, at both creation sites (the Gibbs
  sampler and the stored-tree sampler `predict` restores). No behavior
  change today - `""` is correct for gaussian and probit, and wrong
  only for aft and logistic - but the header names that trap and the
  token is now in hand. `"auto"` and an absent slot both fall back to
  `""`, so a model restored from a pre-1.0 saved fit still creates.

## Not in scope, and why

`missing = "incorporate"` is listed with the other knobs in the survey
but is not one: `dbartsSpec` does not take it, `dbartsData` does, and
the lme4 ingestion drops NA rows through `model.frame`'s na.action
before `dbartsData` is called and again through the union
`na.action.all`. Forwarding the name would do nothing. Left with
dbarts-passthroughs, with the reason recorded in the TODO.

## Gates

- Object equivalence, ahead of the edit: old triple vs `dbartsSpec`'s
  across 9 configurations (continuous/binary-numeric/binary-factor x
  defaults, n.trees, k, `chi` k, cgm parameters, weights, n.cuts), slot
  by slot plus attribute sets. Every slot identical except
  `model@family`, which resolves from the class default "auto" to
  "gaussian"/"probit".
- Bitwise draws: seeded fits (n=120, 2 chains, warmup 15, iter 40) for
  continuous with k/power/base, binary, and weighted continuous, built
  before and after the change and compared with `identical()` on
  `$bart_train`, `$stan`, and `$sigma` - TRUE for all three.
- `compare-posterior.R compare-all posterior-baselines-75d7970.rds`:
  ALL 5 TIERS PASS (continuous nc1/nc2/nc3, weighted, binary), 8/8
  healthy chains each.
- `recompute-exactness.R`: 3/3 tiers PASS, max |stored - recomputed|
  6.2e-14 - the store="trees" path exercises the stored-sampler
  creation site that now passes the family token.
- Suite: 530 expectations, 0 failures (517 before, 13 added to
  test-09-bartArgs.R: the k/node.prior and power/tree.prior
  equivalences, the two shorthand-conflict errors, monotone reaching
  dbarts and its proposal.probs and unknown-variable errors, the
  interactions() error, the four reserved names, and the family
  mismatch).
