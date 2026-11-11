# extract for trees
source(system.file("common", "quiet-mixing.R", package = "stan4bart"), local = TRUE)

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

glCall <- quote(
  stan4bart:::glFormula(
    formula = formula, data = df,
     control = stan4bart:::make_glmerControl(),
     na.action = getOption("na.action", "na.omit")
   )
)

# glFormula correctly identifies component variables
#
# each section below is wrapped in local() because it mutates the shared
# glCall/df scaffolding (glCall$offset, glCall$formula, df$<col>); under
# testthat those mutations were confined to the test_that()'s own
# environment, so local() reproduces that isolation here rather than
# letting one section's edits leak into the next.
local({
  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X3 - X4 - g.1 - g.2 + X3:X4 + (1 + X4 | g.1) + (1 | g.2))

  expect_warning(glResult <- eval(glCall))


  expect_true(setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2")))
  expect_true(setequal(colnames(glResult$bartData@x), c("X1", "X2")))

  expect_true(setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 1:10), "z")))
  expect_true(setequal(colnames(glResult$X), c("X1", "X2", "X3:X4", paste0("X", 5:10), "z")))

  expect_true(setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4")))
  expect_true(setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1")))
})

# glFromula removes offset
local({
  offsetLocalVariable <- rep(c(1:10), length.out = nrow(df))

  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X1 - X2 - g.1 - g.2 + (1 + X4 | g.1) + (1 | g.2))
  glCall$offset <- quote(offsetLocalVariable)

  glResult <- eval(glCall)

  expect_true(setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2")))
  expect_true(setequal(colnames(glResult$bartData@x), c("X1", "X2")))

  expect_true(setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 3:10), "z")))
  expect_true(setequal(colnames(glResult$X), c(paste0("X", 3:10), "z")))

  expect_true(setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4")))
  expect_true(setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1")))

  df$offsetInDataFrame <- offsetLocalVariable
  rm(offsetLocalVariable)
  glCall$offset <- quote(offsetInDataFrame)
  glCall$formula <- quote(y ~ bart(X1 + X2) + . - X1 - X2 - g.1 - g.2 - offsetInDataFrame + (1 + X4 | g.1) + (1 | g.2))
  glResult <- eval(glCall)

  expect_true(setequal(attr(glResult$terms, "varnames.bart"), c("X1", "X2")))
  expect_true(setequal(colnames(glResult$bartData@x), c("X1", "X2")))

  expect_true(setequal(attr(glResult$terms, "varnames.fixed"), c(paste0("X", 3:10), "z")))
  expect_true(setequal(colnames(glResult$X), c(paste0("X", 3:10), "z")))

  expect_true(setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2", "X4")))
  expect_true(setequal(names(glResult$reTrms$Ztlist), c("1 | g.2", "1 + X4 | g.1")))
})

# works when parametric component contains only nested factors
local({
  glCall$formula <- quote(y ~ bart(. - g.1 - g.2 - X4 - z) + (1 | g.1:g.2))
  glResult <- eval(glCall)

  expect_true(setequal(attr(glResult$terms, "varnames.bart"), paste0("X", c(1:3, 5:10))))
  expect_true(setequal(colnames(glResult$bartData@x), paste0("X", c(1:3, 5:10))))

  expect_equal(attr(glResult$terms, "varnames.fixed"), character())
  expect_equal(dim(glResult$X), c(nrow(df), 0L))

  expect_true(setequal(attr(glResult$terms, "varnames.random"), c("g.1", "g.2")))
  expect_equal(names(glResult$reTrms$Ztlist), "1 | g.1:g.2")
})

# works with a factor variable
local({
  df$X10 <- as.factor(LETTERS[1 + round(4 * df$X10, 0)])
  glCall$formula <- quote(y ~ bart(. - g.1 - g.2 - X4 - z) + X10 + (1 | g.1))
  expect_warning(glResult <- eval(glCall))

  expect_true(setequal(attr(glResult$terms, "varnames.bart"), paste0("X", c(1:3, 5:10))))
  # categorical splits: X10 is one bart column, not one indicator per level
  expect_true(setequal(colnames(glResult$bartData@x), paste0("X", c(1:3, 5:10))))

  expect_equal(attr(glResult$terms, "varnames.fixed"), "X10")
  expect_true(setequal(colnames(glResult$X), paste0("X10", levels(df$X10)[-1L])))

  expect_true(setequal(attr(glResult$terms, "varnames.random"), "g.1"))
  expect_equal(names(glResult$reTrms$Ztlist), "1 | g.1")
})

# lme4 is a suggested package: every formula helper is this package's own
# copy, so a build made beside lme4 or reformulas holds none of their
# functions and runs where they are not installed
local({
  ns <- asNamespace("stan4bart")
  home <- vapply(ls(ns, all.names = TRUE), function(name) {
    f <- get(name, envir = ns)
    if (is.function(f) && !is.primitive(f)) environmentName(environment(f)) else ""
  }, "")
  expect_identical(unname(home[c("glFormula", "mkReTrms", "findbars", "subbars")]),
                   rep("stan4bart", 4L))
  expect_identical(names(home)[home %in% c("lme4", "reformulas")], character())
})

# lme4 does the mixed-model fits when it is installed and not turned off
local({
  old <- options(stan4bart.lme4 = NULL)
  on.exit(options(old))
  expect_identical(stan4bart:::lme4_is_used(), requireNamespace("lme4", quietly = TRUE))
  options(stan4bart.lme4 = FALSE)
  expect_false(stan4bart:::lme4_is_used())
})

# a variable that enters bart() only through a transform fits: each fit on
# transformed terms equals the fit on the same columns made by hand, in its
# draws and in predict() on new rows
local({
  set.seed(2)
  n <- 60L
  d <- data.frame(x = rnorm(n), v = rexp(n) + 0.1, f = rnorm(n), g = factor(sample(4L, n, TRUE)))
  d$y <- exp(0.3 * (d$x + log(d$v) + d$f) + rnorm(n, 0, 0.3))
  d$lv <- log(d$v)
  d$v2 <- d$v^2
  d$ly <- log(d$y)
  newd <- d[1:7, ]
  fit_pair <- function(formula, hand) {
    fit <- function(form)
      stan4bart(form, d, chains = 1L, warmup = 5L, iter = 10L, verbose = -1L, seed = 3L,
                bart_args = list(n.trees = 5L, keepTrees = TRUE))
    a <- fit(formula)
    b <- fit(hand)
    expect_equal(unname(extract(a)), unname(extract(b)), info = deparse(formula))
    expect_equal(unname(predict(a, newd)), unname(predict(b, newd)), info = deparse(formula))
  }
  fit_pair(y ~ x + bart(log(v) + f), y ~ x + bart(lv + f))
  fit_pair(y ~ x + bart(log(v)), y ~ x + bart(lv))
  fit_pair(y ~ x + bart(I(v^2) + f), y ~ x + bart(v2 + f))
  fit_pair(y ~ x + v + bart(log(v) + f), y ~ x + v + bart(lv + f))
  fit_pair(y ~ x + (1 | g) + bart(log(v) + f), y ~ x + (1 | g) + bart(lv + f))
  fit_pair(log(y) ~ x + bart(v + f), ly ~ x + bart(v + f))
})

# a formula with no BART predictor is refused in stan4bart's words, before
# dbarts's refusal of an empty design is reached
local({
  for (formula in list(y ~ X1 + (1 | g.1), y ~ X1 + (1 | g.1) + bart(1))) {
    glCall$formula <- formula
    expect_error(eval(glCall), "no bart component detected")
    expect_error(stan4bart(formula, df, iter = 0L), "no bart component detected")
  }
})

# a bart() part whose every column holds one value on the rows kept is refused
# in the same words, with a clause naming the columns; one column that varies
# keeps the part
local({
  df$cst <- 1
  df$cst2 <- 2
  df$f1 <- factor("a")
  for (formula in list(
    y ~ X1 + (1 | g.1) + bart(cst),
    y ~ X1 + (1 | g.1) + bart(cst + cst2),
    y ~ X1 + (1 | g.1) + bart(f1)
  )) {
    glCall$formula <- formula
    expect_error(eval(glCall), "no bart component detected in formula: every bart\\(\\) predictor holds one value")
  }
  glCall$formula <- y ~ X1 + (1 | g.1) + bart(cst + cst2)
  expect_error(eval(glCall), "'cst', 'cst2'")
  glCall$formula <- y ~ X1 + (1 | g.1) + bart(cst + X2)
  expect_silent(eval(glCall))
  # a column constant only on the rows kept
  df$cst <- NULL
  df$X3 <- c(rep(1, 50), rnorm(nrow(df) - 50L))
  glCall$formula <- y ~ X1 + (1 | g.1) + bart(X3)
  expect_silent(eval(glCall))
  df$X1[51:100] <- NA
  expect_error(eval(glCall), "holds one value")
})
