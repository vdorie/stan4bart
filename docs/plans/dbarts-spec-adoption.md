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

  SHADOWING, and how it differs from dbarts. dbarts binds the
  vocabulary to evaluate ONE PRIOR ARGUMENT at a time -
  `parsePriors` evaluates only `matchedCall$tree.prior`,
  `$node.prior`, `$resid.prior`, `$resid.dist` in it (and `xbart`
  only `$node.prior`, with four names bound) - so the shadow never
  reaches an argument whose vocabulary is anything else. This package
  has no such seam: `bart_args` is one list expression carrying both
  priors and plain values, evaluated whole. Binding all eight names
  unconditionally would therefore shadow a caller's own `fixed` or
  `normal` used as a VALUE (`bart_args = list(monotone = fixed)`),
  which the previous single `chi` binding was too obscure a name to
  hit. Resolved by binding only names appearing in CALL position in
  the expression (`called_names`, misc.R): a prior name called as a
  function has no competing reading, and a bare symbol keeps the
  caller's binding. This also narrows the pre-existing `chi` shadow.
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

## The prior evaluation environment, and one asymmetry left alone

Worth recording, since adopting dbartsSpec is what first made any of
this reachable from here. `parsePriors` binds `control`, `data`,
`num.vars`, and `numvars` into its evaluation environment alongside the
prior constructors, so a prior expression can reference them.

Two of the four are load-bearing in dbarts itself, not just available:
`bart2()` and `rbart_vi()` carry `split.probs = 1 / num.vars` as a
DEFAULT ARGUMENT, and `bart()` carries `splitprobs = 1 / numvars`;
those defaults are lazily evaluated in exactly this environment, so
removing either binding breaks the default equiprobable split prior.
Both are documented too - `?dbartsPriors` Details for `num.vars`, and
`?bart`/`?bart2` under `splitprobs` ("'numvars' and 'num.vars' symbols
will be rebound before execution to the number of columns in the model
matrix"). `control` and `data` are the only two that are neither
documented nor referenced by any prior default, which is the inverse of
the collision risk - the names that carry real collision risk are the
ones nothing appears to use.

It is live through this package as of this arc, via the prior arguments
but NOT via the shorthands. Verified:

    bart_args = list(tree.prior = cgm(split.probs = rep(1/num.vars, num.vars)))   # fits
    bart_args = list(split.probs = rep(1/num.vars, num.vars))   # object 'num.vars' not found

The reason is structural: the `cgm` stub returns its own `match.call()`
without forcing arguments, so the inner expression survives unevaluated
to `parsePriors`, while `k`/`power`/`base`/`split.probs` are values
`bart_args` evaluated in the caller's frame before any prior call is
built. NOT levelled, deliberately: the only way to make the shorthand
see `num.vars` is to bind it in `defn_env`, which re-widens exactly the
shadow the call-position fix above narrowed, in service of a dbarts
feature whose future is undecided (VD 2026-07-27: the bindings were
intentional, current usage unknown). Documented instead - `?stan4bart`
now says the prior arguments are expressions and the shorthands are
values.

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
