# the initial fit's share for a model column constant within the categories of
# a BART factor: the slope of the least-squares line through the fit's
# category totals, each category counted once, whatever the order of the
# category names

has_lme4 <- requireNamespace("lme4", quietly = TRUE)

route_of <- function(...) stan4bart:::bart_range_route(...)
share_of <- function(...) stan4bart:::range_factor_share(...)
block_of <- function(object, w = NULL) {
  rows_in <- if (is.null(w)) seq_along(object$y) else which(w > 0)
  stan4bart:::range_fixed_block(object, rows_in,
                                if (is.null(w)) rep_len(1, length(rows_in)) else w[rows_in])
}
# the data under the category names 'names', the categories in the order the
# codes give them
relabel <- function(d, codes, names, column = "r") {
  d[[column]] <- factor(names[codes], levels = sort(names))
  d
}
all_orders <- function(names) {
  if (length(names) <= 1L) return(list(names))
  do.call(c, lapply(seq_along(names), function(i)
    lapply(all_orders(names[-i]), function(rest) c(names[i], rest))))
}
# the slope of the unweighted least-squares line through 'totals' on 'value',
# one point per category
line_slope <- function(totals, value) unname(coef(lm(totals ~ value))[-1L])
# the fit's category totals, computed apart from the package: a weighted
# linear fit with the factor as indicators of every category less the first
# and no column of the model's; the totals are its intercept plus each
# category's coefficient
lm_totals <- function(d, formula, w = NULL, factor = "r") {
  # the weights by value: a name would be looked up in the data first
  d$.w <- if (is.null(w)) 1 else w
  fit <- lm(formula, d, weights = .w)
  cf <- coef(fit)
  levels_in <- levels(droplevels(d[[factor]][if (is.null(w)) TRUE else w > 0]))
  cf[[1L]] + c(0, cf[paste0(factor, levels_in[-1L])])
}

set.seed(616L)
n <- 400L
# 4 regions of unequal size; income is a region's, and the region's own
# deviation from the line through income is large
codes <- sample(4L, n, TRUE, prob = c(0.1, 0.2, 0.3, 0.4))
income_r <- c(1, 2, 4, 7)
d0 <- data.frame(income = income_r[codes], x = runif(n), g = factor(sample(12L, n, TRUE)))
d0$y <- 0.4 * d0$income + c(0, 9, -7, 2)[codes] + 3 * sin(pi * d0$x) +
  rnorm(12L, sd = 2)[d0$g] + rnorm(n)
d0$w <- runif(n, 0.5, 2)
orders <- all_orders(c("A", "B", "C", "D"))

# the region example ----------------------------------------------------------

# the share is the slope of the line through the region means, near it here
# where x is independent of region; the QR's own slope, set by the two regions
# left without an indicator, is far from it and differs by the names' order
local({
  means <- tapply(d0$y, codes, mean)
  target <- line_slope(means, income_r)
  qr_slopes <- shares <- widths <- numeric(length(orders))
  for (i in seq_along(orders)) {
    d <- relabel(d0, codes, orders[[i]])
    object <- stan4bart(y ~ income + bart(r + x), d, iter = 0L)
    columns <- block_of(object)
    fit <- stan4bart:::range_fit(object$y, rep_len(1, n), columns$X)
    qr_slopes[i] <- fit$beta[columns$own]
    route <- route_of(object, "lm")
    shares[i] <- route$coefficients[["income"]]
    widths[i] <- diff(route$range)
  }
  # PRECONDITION: the QR's slope depends on the order and is not the share
  expect_true(diff(range(qr_slopes)) > 1)
  expect_true(min(abs(qr_slopes - target)) > 0.5)
  expect_true(abs(shares[1L] - target) < 0.05)
  expect_true(diff(range(shares)) < 1e-10)
  expect_true(diff(range(widths)) < 1e-10 * widths[1L])
})

# the share's definition and its order-free range --------------------------

# for every order of the names, the coefficient is the line through the
# weighted linear fit's category totals, and the range is that of the
# response less the column times that slope
local({
  reference <- NULL
  for (names in orders) {
    d <- relabel(d0, codes, names)
    object <- stan4bart(y ~ income + bart(r + x), d, weights = w, iter = 0L)
    route <- route_of(object, "lm")
    totals <- lm_totals(d, y ~ x + r, d$w)
    slope <- line_slope(totals, income_r[match(sort(names), names)])
    expect_equal(route$coefficients[["income"]], slope, tolerance = 1e-8)
    own <- (d$income - object$X_means[["income"]]) * slope
    names(own) <- NULL
    expect_equal(route$range, range(d$y - own), tolerance = 1e-8)
    expect_equal(unname(route$own), own, tolerance = 1e-8)
    if (is.null(reference)) reference <- route
    expect_equal(route$range, reference$range, tolerance = 1e-10)
    expect_equal(route$sigma, reference$sigma, tolerance = 1e-10)
  }
})

# counted once: the line is not that of the rows or of the weights
local({
  d <- relabel(d0, codes, orders[[1L]])
  object <- stan4bart(y ~ income + bart(r + x), d, weights = w, iter = 0L)
  totals <- lm_totals(d, y ~ x + r, d$w)
  by_rows <- unname(coef(lm(totals ~ income_r, weights = tabulate(codes)))[2L])
  by_weights <- unname(coef(lm(totals ~ income_r, weights = tapply(d$w, codes, sum)))[2L])
  share <- route_of(object, "lm")$coefficients[["income"]]
  expect_equal(share, line_slope(totals, income_r), tolerance = 1e-8)
  # PRECONDITION: the counts and the weights move the line
  expect_true(abs(by_rows - share) > 0.1 && abs(by_weights - share) > 0.1)
  # weights scaled within one region move the totals and not the counting:
  # the share is the line through the new totals
  w2 <- d$w * ifelse(codes == 4L, 20, 1)
  d$w2 <- w2
  object2 <- stan4bart(y ~ income + bart(r + x), d, weights = w2, iter = 0L)
  expect_equal(route_of(object2, "lm")$coefficients[["income"]],
               line_slope(lm_totals(d, y ~ x + r, w2), income_r), tolerance = 1e-8)
})

# the model's other columns, a numeric forest column constant within the
# categories, and the default route
local({
  d <- relabel(d0, codes, orders[[1L]])
  d$z <- rnorm(n)
  d$pop <- c(5, 1, 3, 2)[codes]
  reference <- NULL
  for (names in orders[c(1L, 8L, 17L, 24L)]) {
    d <- relabel(d, codes, names)
    object <- stan4bart(y ~ income + z + bart(pop + r + x), d, iter = 0L)
    route <- route_of(object, "mixed")
    totals <- lm_totals(d, y ~ z + x + r)
    expect_equal(route$coefficients[["income"]],
                 line_slope(totals, income_r[match(sort(names), names)]), tolerance = 1e-8)
    # the other model column is the fit's, untouched by the share
    expect_equal(route$coefficients[["z"]], unname(coef(lm(y ~ z + x + r, d))[["z"]]),
                 tolerance = 1e-8)
    if (is.null(reference)) reference <- route
    expect_equal(route$range, reference$range, tolerance = 1e-10)
  }
})

# two model columns constant within one factor share one fit, the plane
# through the totals
local({
  d <- d0
  d$density <- c(3, 1, 2, 2.5)[codes]
  for (names in orders[c(1L, 10L, 24L)]) {
    d <- relabel(d, codes, names)
    object <- stan4bart(y ~ income + density + bart(r + x), d, iter = 0L)
    route <- route_of(object, "lm")
    permuted <- match(sort(names), names)
    totals <- lm_totals(d, y ~ x + r)
    plane <- coef(lm(totals ~ income_r[permuted] + c(3, 1, 2, 2.5)[permuted]))
    expect_equal(unname(route$coefficients[c("income", "density")]), unname(plane[-1L]),
                 tolerance = 1e-8)
  }
})

# random effects ------------------------------------------------------------

# under "mixed" the totals are those of the mixed model, net of the predicted
# random effects, and the range is the same for every order
if (has_lme4) local({
  control <- lme4::lmerControl(check.conv.grad = "ignore", check.conv.singular = "ignore",
                               check.conv.hess = "ignore")
  reference <- NULL
  for (names in orders[c(1L, 6L, 11L, 19L, 24L)]) {
    d <- relabel(d0, codes, names)
    object <- stan4bart(y ~ income + (1 | g) + bart(r + x), d, weights = w, iter = 0L)
    route <- route_of(object, "mixed")
    expect_identical(route$info$fit, "lme4")
    oracle <- lme4::lmer(y ~ x + r + (1 | g), d, weights = w, control = control)
    cf <- lme4::fixef(oracle)
    totals <- cf[[1L]] + c(0, cf[paste0("r", sort(names)[-1L])])
    expect_equal(route$coefficients[["income"]],
                 line_slope(totals, income_r[match(sort(names), names)]), tolerance = 1e-4)
    if (is.null(reference)) reference <- route
    expect_equal(route$range, reference$range, tolerance = 1e-6)
    expect_equal(route$own, reference$own, tolerance = 1e-6)
  }
})

# a factor nested in a grouping factor is left out, so nothing is dropped and
# the coefficient is the fit's own
if (has_lme4) local({
  d <- relabel(d0, codes, orders[[1L]])
  d$big <- factor(c(1, 1, 2, 2)[codes])
  object <- stan4bart(y ~ income + (1 | big) + bart(r + x), d, iter = 0L)
  columns <- block_of(object)
  expect_identical(columns$left.out, "r")
  beta <- seq_len(ncol(columns$X)) + 0.5
  expect_identical(share_of(columns, beta), beta)
})

# empty and single-row categories -------------------------------------------

# a level with no row, or none of positive weight, is no category; a category
# of one row counts once, as any other
local({
  d <- relabel(d0, codes, orders[[1L]])
  d$r <- factor(as.character(d$r), levels = c("A", "B", "C", "D", "E"))
  object <- stan4bart(y ~ income + bart(r + x), d, iter = 0L)
  totals <- lm_totals(d, y ~ x + r)
  expect_equal(route_of(object, "lm")$coefficients[["income"]], line_slope(totals, income_r),
               tolerance = 1e-8)

  # region D's rows at weight zero: three categories in
  d$w <- ifelse(codes == 4L, 0, 1)
  object <- stan4bart(y ~ income + bart(r + x), d, weights = w, iter = 0L)
  totals <- lm_totals(d[d$w > 0, ], y ~ x + r)
  expect_equal(route_of(object, "lm")$coefficients[["income"]],
               line_slope(totals, income_r[1:3]), tolerance = 1e-8)

  # a fifth region of one row
  d <- relabel(d0, codes, orders[[1L]])
  d <- rbind(d, d[1L, ])
  d$r <- factor(c(as.character(d$r[-nrow(d)]), "E"))
  d$income[nrow(d)] <- 10
  d$y[nrow(d)] <- 30
  object <- stan4bart(y ~ income + bart(r + x), d, iter = 0L)
  totals <- lm_totals(d, y ~ x + r)
  expect_equal(route_of(object, "lm")$coefficients[["income"]],
               line_slope(totals, c(income_r, 10)), tolerance = 1e-8)
  # the one row's total is its residual after x: with its response moved,
  # its point on the line moves by as much
  d$y[nrow(d)] <- 40
  object <- stan4bart(y ~ income + bart(r + x), d, iter = 0L)
  expect_equal(route_of(object, "lm")$coefficients[["income"]],
               line_slope(lm_totals(d, y ~ x + r), c(income_r, 10)), tolerance = 1e-8)
})

# nested BART factors -------------------------------------------------------

# a column constant within regions, with states nested in regions also in
# the forest, takes the line through the state totals, each state once
local({
  state <- (codes - 1L) * 3L + sample(3L, n, TRUE)
  d <- d0
  d$s <- factor(sprintf("s%02d", state))
  reference <- NULL
  for (names in orders[c(1L, 13L, 24L)]) {
    d <- relabel(d, codes, names)
    for (state_names in list(sprintf("s%02d", 1:12), sprintf("s%02d", 12:1))) {
      d$s <- factor(state_names[state], levels = sort(state_names))
      object <- stan4bart(y ~ income + bart(r + s + x), d, iter = 0L)
      expect_identical(names(block_of(object)$shares), "s")
      route <- route_of(object, "lm")
      totals <- lm_totals(d, y ~ x + s, factor = "s")
      income_s <- tapply(d$income, d$s, unique)
      expect_equal(route$coefficients[["income"]], line_slope(totals, income_s),
                   tolerance = 1e-8)
      if (is.null(reference)) reference <- route
      expect_equal(route$range, reference$range, tolerance = 1e-10)
    }
  }
})

# category values that are not exact in binary --------------------------

# the categories' means of such a column round, and it is still constant
# within them
local({
  for (values in list(c(31.7, 45.3, 52.9, 70.1), c(0.1, 0.3, 0.7, 1.1))) {
    d <- d0
    d$income <- values[codes]
    d$y <- d$y + 0.3 * d$income
    shares <- numeric()
    for (names in orders[c(1L, 8L, 17L, 24L)]) {
      d <- relabel(d, codes, names)
      object <- stan4bart(y ~ income + bart(r + x), d, weights = w, iter = 0L)
      expect_identical(names(block_of(object, d$w)$shares), "r")
      share <- route_of(object, "lm")$coefficients[["income"]]
      expect_equal(share, line_slope(lm_totals(d, y ~ x + r, d$w),
                                     values[match(sort(names), names)]), tolerance = 1e-8)
      shares <- c(shares, share)
    }
    expect_true(diff(range(shares)) < 1e-10 * max(1, abs(shares)))
  }
})

# other columns that combine to a constant ---------------------------------

# where the other columns meet the factor's category functions in a constant
# alone, the share is made; each case is checked against a fit in which the
# factor's indicators are replaced by a basis of the category functions
# orthogonal, one point per category, to the constant and the shared columns,
# so that the shared columns' coefficients are the slopes of the line
local({
  orthogonal_basis <- function(codes, values) {
    K <- max(codes)
    A <- cbind(1, values)
    Q <- qr.Q(qr(cbind(A, diag(K))))[, -seq_len(ncol(A)), drop = FALSE]
    Q[codes, seq_len(K - ncol(A)), drop = FALSE]
  }
  set.seed(3L)
  m <- 600L
  cf <- sample(6L, m, TRUE)
  cg <- sample(5L, m, TRUE)
  vf1 <- c(1, 4, 2, 8, 3, 5.5)
  vf2 <- c(0, 1, 1, 0, 2, 3.1)
  vg <- c(10, 12, 9, 15, 11.3)
  d <- data.frame(a1 = vf1[cf], a2 = vf2[cf], b = vg[cg], x = rnorm(m),
                  pg = c(3, 1, 4, 1, 5)[cg], sex = factor(sample(c("F", "M"), m, TRUE)))
  d$y <- d$a1 + 2 * d$a2 - 0.5 * d$b + rnorm(6L, sd = 2)[cf] + rnorm(5L, sd = 2)[cg] + d$x +
    rnorm(m)
  cases <- list(
    # a second model column constant within a crossed factor in the forest
    list(formula = y ~ a1 + a2 + b + bart(f + gg + x), shares = c("f", "gg"),
         oracle = function(d) lm(y ~ a1 + a2 + b + x + orthogonal_basis(cf, cbind(vf1, vf2)) +
                                   orthogonal_basis(cg, vg), d),
         shared = c("a1", "a2", "b")),
    # a numeric forest column constant within a crossed forest factor
    list(formula = y ~ a1 + bart(f + gg + pg + x), shares = "f",
         oracle = function(d) lm(y ~ a1 + x + pg + orthogonal_basis(cf, vf1) + gg, d),
         shared = "a1"),
    # every level of a model factor, which sum to the intercept
    list(formula = y ~ 0 + sex + a1 + bart(f + x), shares = "f",
         oracle = function(d) lm(y ~ sex + a1 + x + orthogonal_basis(cf, vf1), d),
         shared = "a1"))
  name_sets <- list(list(f = letters[1:6], g = LETTERS[1:5]), list(f = letters[6:1], g = LETTERS[5:1]),
                    list(f = c("c", "a", "f", "b", "e", "d"), g = c("D", "B", "E", "A", "C")))
  for (case in cases) {
    reference <- NULL
    for (names in name_sets) {
      d$f <- factor(names$f[cf], levels = sort(names$f))
      d$gg <- factor(names$g[cg], levels = sort(names$g))
      object <- stan4bart(case$formula, d, iter = 0L)
      expect_identical(names(block_of(object)$shares), case$shares)
      route <- route_of(object, "lm")
      expect_equal(unname(route$coefficients[case$shared]),
                   unname(coef(case$oracle(d))[case$shared]), tolerance = 1e-8)
      if (is.null(reference)) reference <- route
      expect_equal(route$range, reference$range, tolerance = 1e-10)
    }
  }
})

# the limits ----------------------------------------------------------------

# a column constant within two factors' joint categories, as the sum of a
# part of each, is in their indicators' span without being constant within
# either: no share is made, and the coefficient and the range follow the
# order of the names
local({
  other <- sample(3L, n, TRUE)
  d <- d0
  d$income <- income_r[codes] + c(0, 3, 5)[other]
  d$y <- d$y + c(1, -3, 2)[other]
  slopes <- numeric()
  for (names in all_orders(c("t", "u", "v"))) {
    d <- relabel(relabel(d, codes, orders[[1L]]), other, names, "q")
    object <- stan4bart(y ~ income + bart(r + q + x), d, iter = 0L)
    columns <- block_of(object)
    beta <- seq_len(ncol(columns$X)) + 0.5
    expect_identical(share_of(columns, beta), beta)
    slopes <- c(slopes, route_of(object, "lm")$coefficients[["income"]])
  }
  expect_true(diff(range(slopes)) > 0.1)
})

# a column constant within each of two forest factors neither nested in the
# other is left as the QR leaves it
local({
  pair <- sample(2L, n, TRUE)
  first <- ifelse(pair == 1L, sample(1:3, n, TRUE), sample(4:6, n, TRUE))
  second <- ifelse(pair == 1L, sample(1:2, n, TRUE), sample(3:4, n, TRUE))
  d <- data.frame(y = c(1, 4)[pair] + first / 3 + rnorm(n), x = runif(n),
                  income = c(2, 5)[pair], fa = factor(first), fb = factor(second))
  object <- stan4bart(y ~ income + bart(fa + fb + x), d, iter = 0L)
  columns <- block_of(object)
  expect_identical(length(columns$shares), 0L)
  beta <- seq_len(ncol(columns$X)) + 0.5
  expect_identical(share_of(columns, beta), beta)

  # a column constant within the first factor's categories alone, where both
  # factors' categories lie within those of the pair: the second factor's
  # indicators carry the pair, which is constant within the first's
  # categories, so the first takes no share, and the coefficient follows the
  # order of the second's names
  d$income <- c(1, 2, 4, 3, 6, 8)[first]
  slopes <- numeric()
  for (names in list(c("p", "q", "r", "s"), c("s", "r", "q", "p"), c("q", "s", "p", "r"))) {
    d$fb <- factor(names[second], levels = sort(names))
    object <- stan4bart(y ~ income + bart(fa + fb + x), d, iter = 0L)
    expect_identical(length(block_of(object)$shares), 0L)
    slopes <- c(slopes, route_of(object, "lm")$coefficients[["income"]])
  }
  expect_true(diff(range(slopes)) > 0.01)
})

# no such column ------------------------------------------------------------

# where no model column is constant within a factor's categories the share
# leaves the coefficients as they are, bit for bit
local({
  d <- relabel(d0, codes, orders[[1L]])
  d$income <- d$income + rnorm(n, sd = 0.1)
  for (formula in list(y ~ income + bart(r + x), y ~ income + (1 | g) + bart(r + x),
                       y ~ x + bart(r + income), y ~ (1 | g) + bart(r + x))) {
    object <- stan4bart(formula, d, weights = w, iter = 0L)
    columns <- block_of(object, d$w)
    beta <- rnorm(ncol(columns$X))
    expect_identical(share_of(columns, beta), beta)
  }
  # nor where the forest has no factor
  object <- stan4bart(y ~ income + bart(x), d, iter = 0L)
  columns <- block_of(object)
  beta <- rnorm(ncol(columns$X))
  expect_identical(share_of(columns, beta), beta)
})
