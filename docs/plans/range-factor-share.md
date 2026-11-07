# range-factor-share

Status: built, in review. Branch range-factor-share from main 27ef163. Rulings: dbarts dec-B616 (the share), dec-A265
call 41 (the guard); context dec-B441 to dec-B447 and dbarts docs/design/stan4bart-response-range.md.

## What changes

1. `range_fixed_block` (R/initial_fit.R) records, for each BART factor it enters as indicators, which kept model
   columns are constant within the factor's categories on the rows in, and makes a share for that factor only where
   the rank condition below holds. `range_factor_share` then replaces those columns' coefficients, after the fit and
   before the model's own part is formed, by the slopes of the least-squares fit, one point per category, of the
   fit's category totals on the columns' category values. The fit itself, its predicted random effects, its residual
   sd and every other coefficient are unchanged. Both estimated routes ("mixed" and "lm") and a supplied pair's
   start go through it; a binary start has no forest columns and is untouched.
2. `glFormula` (R/lme4_functions.R) checks for a BART design with no term before calling `dbartsData`, so
   `stan4bart(y ~ x + (1 | g))` stops with stan4bart's "no bart component detected" again, not with dbarts's
   "predictors must be specified for regression tree analysis". The check after the call is gone: dbartsData now
   refuses every design it would have caught.
3. Help (`bart_range`): the share, its limits, and the cost of a many-category factor. NEWS: one sentence in the
   0.0-14 `bart_range` entry, the feature being new in 0.0-14.

## The share's argument

Notation. The initial fit reads the n rows in (kept, positive weight), with weights w_i and target z_i (the response
less any offset). The fixed block is X = [1, the model's own columns, the forest's columns], a BART factor F entering
as indicators I_2, ..., I_K of its K categories present among the rows in, the first left out. The fit is weighted
least squares ("lm", or "mixed" without random terms) or lme4's REML mixed model with fixed design X. A pivoted QR at
tolerance 1e-7 keeps the first column of any dependent set, the model's columns before the forest's.

What goes wrong. Let c be a model column constant within F's categories, c_i = v_k for row i in category k. Then
c = v_1 + sum_k (v_k - v_1) I_k lies in the span of the intercept and F's indicators. The QR keeps c (it comes first)
and drops one more of F's indicators, the one whose position makes the set dependent; that position is set by the
order of the category names. With the indicator of category m dropped, the fit's part constant within the categories
is t_k = a + b v_k + g_k with g_1 = g_m = 0, so b = (t_m - t_1) / (v_m - v_1): the contrast between the two
categories left without an indicator. The own part, b (c - mean c), and the range read from it follow the names.

What does not depend on the order. Let V be the K-dimensional space of functions of F's category (constant within
each category); it contains 1, c and every I_k, so V lies in the span of X however the names fall. Relabeling the
categories replaces the set {I_2, ..., I_K} by another basis of the same V less its constants, so the span of the
whole block is the same. Weighted least squares depends on X only through that span, and so does REML: replacing X
by X A for an invertible A shifts the restricted criterion by the constant log det(A' A) and leaves its optimum, the
predicted random effects and the fitted values unchanged. So the fit's values f = X beta, its random effects and its
residual sd are the same for every order; only how f is written in the kept columns changes.

The rank condition. Let S be the block's columns constant within F's categories (the intercept, c and any other such
model column, F's indicators, a numeric forest column taking one value per category, the indicators of a forest factor
whose categories each lie within one of F's), and N the rest. Condition (C): no nonzero combination of the columns of
N is constant within F's categories, that is span(N) and V meet only at 0. The code tests it as rank(block) =
K + rank(N), both ranks by the QR at the fit's tolerance under the fit's weights; the inequality <= always holds and
equality is (C). Under (C):

- span(block) = V + span(N), a direct sum, and neither summand depends on the order of the names (V is F's space,
  N's columns are named and placed independently of F's names; a second factor's indicators in N are another basis of
  its own space under its names, as above).
- The kept columns of S span V and the kept columns of N span span(N): a column dropped as dependent on earlier kept
  columns has, by the direct sum, its V part in the kept S columns and its N part in the kept N columns, and the
  counts add up to the rank.
- So f splits in exactly one way as f = t + r, t in V and r in span(N). t, the category totals, is f's V part, the
  same for every order; it is computed as the kept S columns' coefficients times their values at one row of each
  category.

The share. Given the totals t_k and the category values v_k of the shared columns (one column: a line; several
constant within the same F: a plane, one least-squares fit), the share is the slope of the unweighted least-squares
fit of t_k on [1, v_k] over the K categories. It is a function of the K pairs (t_k, v_k), which a relabeling only
permutes, and least squares does not depend on the order of its points: the share is order-free. The intercept of the
line and the deviations d_k = t_k - (a + s v_k) are the forest's part, as the indicators' were; only s (c - mean c)
enters the model's own part. With one categorical column and nothing else in the model, t_k is the category's
weighted mean of z, so the share is the slope of the line through the category means, each category counted once.
Because the totals depend only on f, nothing about the fit changes; because each category is one point, the counting
is that of the ruling: categories, not rows or weight.

Weights. They enter the fit (and so the totals, which are the weighted fit's) and the QR's rank decisions, which with
positive weights are those of the unweighted block. They do not enter the line: a category's point has weight one
whatever its rows' weights. Rows of weight zero are not rows in, so a category held only by such rows is no category.

Random effects. Under "mixed" the totals are those of the fixed part of the mixed fit, net of the predicted random
effects, which REML's invariance above makes order-free. A factor whose levels each lie within one level of a
grouping factor is left out of the fit (an older rule), so it has no indicators and no share; a model column constant
within its categories then keeps the fit's own coefficient, which is identified.

Empty and single-row categories. A level with no row in has no indicator and is no point of the line. A category with
one row in has an indicator that fits that row's fixed part exactly; its total is identified and is one point, the
same weight as any other category. It has as much leverage on the line as a large category, which is what counting
each category once means.

Limits, left as the QR leaves them and named in the help: a column constant only within two factors' joint
categories (constant within neither, so there is no candidate; such a column is in the span of the two factors'
indicators only when it is additive in them, which is when the QR drops an indicator); a column constant within a
factor's categories only in combination with other columns; and any factor failing (C), such as two forest factors
whose categories both lie within those of a coarser grouping that the model does not name (the second factor's
indicators then carry that grouping, which is constant within the first factor's categories). A column constant
within two factors, one nested in the other, takes the finer factor's share: the coarser one fails (C), since the
finer one's indicators are in its N and span its V. Of the factors a column is constant within, only those with
the finest categories can pass (C): if c is constant within F and G and G passes, F's indicators cannot be in G's N
(c less a constant is a combination of them and is constant within G's categories), so F is constant within G's
categories; if F's categories are coarser, G's indicators are in F's N and carry c, and F fails. Two factors with the
same categories under different names both pass, and give the same totals and the same share.

What the check costs. Nothing where no kept model column is constant within an entered factor's categories (one pass
over the model's kept columns per factor). Where one is: one pass over the block and a QR of N. At 20,000 rows and a
500-category factor, "lm", 8.3 s against main's 8.0 s on this machine.

## Tests (inst/tinytest/test-26-factor-share.R, 205 results; test-08-glFormula.R, 4 more)

- The region example: 4 regions of unequal size, 400 rows, income per region. Over all 24 orders of the names the
  share is the same to 1e-10, as is the width; it is within 0.05 of the line through the raw region means (-0.086)
  where the QR's own slope runs from -7.6 to 9.3 by order (a precondition) and the width from 24.2 to 62.6 against
  the response's 26.2 (so on main the cap took over in some orders and not others).
- The definition, for all 24 orders under weights: the coefficient equals the line through the category totals of
  an independent weighted `lm(y ~ x + r)`, the range and the own part equal those formed from it.
- Counted once: the share is not the row-weighted or weight-weighted line (preconditions that those differ by more
  than 0.1); weights raised twenty-fold in one region move the totals and the share follows the new totals.
- The model's other columns untouched; a numeric forest column constant within the categories; the default route.
- Two columns of one factor: one plane, against `lm` on the totals.
- Random effects: "mixed" under weights for 5 orders, against `lmer(y ~ x + r + (1 | g))`'s totals to 1e-4, ranges
  and own parts the same across orders to 1e-6. A factor nested in a grouping factor is left out and takes no share.
- An unused level, a region of weight zero (three categories), a one-row region, and that region's response moved.
- Nested forest factors (states within regions, income per region): the line through the state totals, each state
  once, the same over 3 region orders by 2 state orders; the share is the state factor's.
- Limits: a column additive in two factors (constant within their joint categories) takes no share and its slope
  moves with the second factor's names; two factors crossed under a common coarser grouping take no share, with the
  column constant within both and with it constant within one, the slope then moving with the other's names.
- No such column: the coefficients come back identical, for four formulas under weights and one with no forest
  factor.
- The guard: `y ~ X1 + (1 | g.1)` and `... + bart(1)` stop with "no bart component detected", from glFormula and
  from stan4bart.

Bit for bit, outside the suite: main's `R/initial_fit.R` sourced beside the branch gave `identical()` route results
in 24 of 24 fits with no such column (6 formulas, "mixed" and "lm", with and without weights and offset, test-25's
data) and differed in the 4 of 4 that have one.

Mutants (each made, installed, run, reverted and touched): see Gates.

## Gates

At e36d44b (the code's last commit; this record follows it), dbarts main f79874ea8 (shipped mode), lme4 2.1.0, R 4.6.2,
arm64 macOS:

- Full suite, NOT_CRAN=true, at_home = TRUE, one process: 1919 results, 0 failures, of which this branch adds 209
  (test-26's 205, test-08's 4); test-25, 506 of 506.
- R CMD check --as-cran on the built tarball: 1 NOTE, CRAN incoming feasibility (the Description's arXiv form, the
  Date field, days since the last update), main's own. A first run found a second NOTE, `lm.fit` not imported,
  repaired in e36d44b.
- Lint and format: stan4bart's CI runs none (its workflows are R CMD check, the exactness and posterior gates, rhub
  and the sanitizers). The exactness gate, which runs on every push, passes 4 of 4 tiers within 1e-10.
- Mutants, each installed, run against the test file and reverted: the line weighted by each category's rows
  (`lm.wfit` with the counts) fails 98 of 205 in test-26; the share skipped (the QR's coefficient kept) fails 116 of
  205; the guard removed fails 4 of 40 in test-08.

## Calls made

1. The share is computed after the fit from its coefficients, not by reparameterizing the block (F's indicators
   replaced by a basis orthogonal, over categories, to [1, v]). The two give the same number; the first leaves the
   fit and its QR exactly as on main. Alternative: the reparameterized block.
2. Several model columns constant within one factor's categories take one least-squares fit (a plane), not one line
   each. With one column it is the ruling's line; separate lines would count a correlated pair's slope twice.
   Alternative: a line per column.
3. The share is made only under the rank condition (C), which also decides between nested factors (the finer). The
   first build took the factor with the most categories where the others were nested in it and skipped crossed ones;
   that missed a factor crossed with another under a common coarser grouping, where the share would have been applied
   and still followed the other factor's names. Alternative: that rule, or the share applied wherever a column is
   constant, order-free or not.
4. "Constant within" is to a tolerance of sqrt(machine epsilon) times the column's largest absolute value, well inside
   the QR's 1e-7. Alternative: exact equality, which a column computed per row (log of a merged value) still meets
   but a rounding difference would not.
5. The check is vectorized over the block from category sums. One column at a time added 73 percent to the fit's time
   at 20,000 rows and 500 categories; vectorized, 4 percent.
6. range.info gains no field naming the shared columns. Its `dropped` still names the indicator the QR dropped, which
   follows the names; it describes the QR, not the range. Alternative: a `shared` field.
7. The guard tests the BART terms' labels, the condition dbartsData refuses on; the check after the call is removed
   as unreachable. Alternative: keep both.
8. NEWS: one sentence in the 0.0-14 `bart_range` entry, not a bug-fix line, since the range is new in 0.0-14. The
   guard has no NEWS line: the message is the one 0.0-13 gave.

## Deviations and open points

- The ruling's help line says giving `bart_range` avoids the quick fit's cost. A pair does not: it still makes the
  linear fit, with the indicator block, for the first sweep's start. The help says `"response"` avoids it and a pair
  still makes the linear fit. Recommendation: keep the help as written; if a pair should avoid the cost too, its start
  would have to come from the model's own columns alone, a change for its own ruling.
- Code is about 45 lines against the ruling's 30, the rank condition and the vectorized check being the difference.
- Found, not in scope: `stan4bart(y ~ x + bart(log(v) + f), d)` stops with "object 'v' not found", where
  `bart(u + log(u + 1))` works; a transformed BART term whose variable is not also a plain term seems to fail in the
  frame handed to dbartsData. Seen on the branch; the path is main's, untouched here.

## Landing note

Not landed.
