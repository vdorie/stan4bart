# stan4bart na.action

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

df.train <- df[seq_len(floor(0.8 * nrow(df))),]
df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]

df.train$X1[1L:3L] <- NA
df.train$X4[4L] <- NA
df.train$g.1[5L:6L] <- NA

df.test$X1[1L:3L] <- NA
df.test$X4[4L] <- NA
df.test$g.1[5L:6L] <- NA


# na.action = na.omit works correctly
fit <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) +
                       X4 + z + (1 + X4 | g.1) + (1 | g.2),
                 df.train, test = df.test,
                 cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                 bart_args = list(n.trees = 11),
                 na.action = na.omit)

types <- c("ev", "indiv.fixef", "indiv.ranef", "indiv.bart")
for (type in types) {
  value <- fitted(fit, type = type)

  expect_true(length(value) == nrow(df.train) - 6L)
  expect_true(!anyNA(value))

  value <- fitted(fit, type = type, sample = "test")

  expect_true(length(value) == nrow(df.test) - 6L)
  expect_true(!anyNA(value))
}

# na.action = na.exclude works correctly
fit <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) +
                       X4 + z + (1 + X4 | g.1) + (1 | g.2),
                 df.train, test = df.test,
                 cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                 bart_args = list(n.trees = 11),
                 na.action = na.exclude)

value <- fitted(fit, type = "indiv.ranef")
expect_true(length(value) == nrow(df.train))
# 3 NA: two missing group and one missing group predictor
expect_true(sum(is.na(value)) == 3L)

# This should complain as we had to fit with data minus 6 rows and
# adding the three ones for which predictions can be made to the
# test data is too much work.
expect_warning(value <- fitted(fit, type = "indiv.bart"))
expect_true(length(value) == nrow(df.train))
expect_true(sum(is.na(value)) == 6L)

value <- fitted(fit, type = "indiv.fixef")
expect_true(length(value) == nrow(df.train))
expect_true(sum(is.na(value)) == 1L)

value <- fitted(fit)
expect_true(length(value) == nrow(df.train))
expect_true(sum(is.na(value)) == 6L)


value <- fitted(fit, type = "indiv.ranef", sample = "test")
expect_true(length(value) == nrow(df.test))
expect_true(sum(is.na(value)) == 3L)

value <- fitted(fit, type = "indiv.bart", sample = "test")
expect_true(length(value) == nrow(df.test))
# With the test data, we can give it a full data frame containing
# all possible rows, whether or not we would have to drop them
# elsewhere.
expect_true(sum(is.na(value)) == 3L)

value <- fitted(fit, type = "indiv.fixef", sample = "test")
expect_true(length(value) == nrow(df.test))
expect_true(sum(is.na(value)) == 1L)

value <- fitted(fit, sample = "test")
expect_true(length(value) == nrow(df.test))
expect_true(sum(is.na(value)) == 6L)

# a missing value in one component alone -----------------------------------

# The response, every design, the weights and the offset are cut to the rows
# no component dropped, whichever component held the missing value: the fit
# is the fit to the complete rows.
local({
  set.seed(31L)
  n <- 80L
  d <- data.frame(g = factor(sample(8L, n, TRUE)), x = rnorm(n), s = runif(n), u = runif(n),
                  v = runif(n), w = runif(n, 0.5, 2), o = rnorm(n))
  d$y <- 2 * d$x + 3 * d$u + (rnorm(8L) + d$s * rnorm(8L))[d$g] + d$o + rnorm(n)
  dropped <- c(7L, 40L)
  kept <- setdiff(seq_len(n), dropped)

  fit_rows <- function(data, with_wo, na.action = na.omit, dot = FALSE, ...) {
    formula <- if (dot) y ~ x + (1 + s | g) + bart(. - g - x - s - w - o - y - extra)
               else y ~ x + (1 + s | g) + bart(u + v)
    if (with_wo)
      stan4bart(formula, data, weights = w, offset = o, na.action = na.action, ...)
    else
      stan4bart(formula, data, na.action = na.action, ...)
  }
  short <- function(...)
    fit_rows(..., cores = 1L, verbose = -1L, chains = 1L, warmup = 3L, iter = 7L, seed = 9L,
             bart_args = list(n.trees = 5L, keepTrees = TRUE))

  # the forest's column, the random slope's, a fixed column, and a column
  # the forest's formula names only to leave out
  cases <- c(forest = "u", slope = "s", fixed = "x", dot = "extra")
  for (case in names(cases)) for (with_wo in c(FALSE, TRUE)) {
    info <- paste0(case, if (with_wo) ", weights and offset")
    holed <- d
    if (case == "dot") holed$extra <- holed$v
    holed[[cases[[case]]]][dropped] <- NA
    complete <- holed[kept, ]
    dot <- case == "dot"

    aligned <- TRUE
    for (na.action in list(na.omit, na.exclude)) {
      fit0 <- fit_rows(holed, with_wo, na.action, dot, iter = 0L)
      aligned <- aligned && nrow(fit0$X) == length(fit0$y) && ncol(fit0$reTrms$Zt) == length(fit0$y)
      expect_identical(as.vector(attr(fit0$frame, "na.action.all")), dropped, info = info)
      expect_identical(c(length(fit0$y), nrow(fit0$X), ncol(fit0$reTrms$Zt),
                         length(fit0$reTrms$flist[[1L]]), nrow(fit0$bartData@x),
                         length(fit0$bartData@y)),
                       rep(length(kept), 6L), info = info)
      expect_identical(c(length(fit0$weights), length(fit0$offset)),
                       if (with_wo) rep(length(kept), 2L) else c(0L, 0L), info = info)
      expect_identical(fit0$y, d$y[kept], info = info)
      expect_identical(unname(fit0$X[, "x"]), d$x[kept], info = info)
      expect_identical(as.character(fit0$reTrms$flist[[1L]]), as.character(d$g[kept]), info = info)
      expect_identical(unname(fit0$bartData@x[, "u"]), d$u[kept], info = info)
      if (with_wo) {
        expect_identical(fit0$weights, d$w[kept], info = info)
        expect_identical(fit0$offset, d$o[kept], info = info)
      }
      reference0 <- fit_rows(complete, with_wo, na.action, dot, iter = 0L)
      expect_identical(unname(fit0$X), unname(reference0$X), info = info)
      expect_identical(as.matrix(fit0$reTrms$Zt), as.matrix(reference0$reTrms$Zt), info = info)
      expect_identical(unname(fit0$bartData@x), unname(reference0$bartData@x), info = info)
    }

    # pieces on different rows are not fitted: the expectations above say so
    if (!aligned) next

    # the draws are those of the complete rows
    reference <- short(complete, with_wo, dot = dot)
    fit <- short(holed, with_wo, dot = dot)
    expect_identical(unname(fit$stan), unname(reference$stan), info = info)
    expect_identical(unname(fit$bart_train), unname(reference$bart_train), info = info)
    expect_identical(fit$range.bart, reference$range.bart, info = info)
    expect_identical(extract(fit, "ev"), extract(reference, "ev"), info = info)
    expect_identical(residuals(fit), residuals(reference), info = info)
    expect_identical(predict(fit), extract(reference, "ev"), info = info)

    # under na.exclude too, reported on the frame's rows: the expected value,
    # the residuals and the predictive draws are missing on the dropped rows,
    # and a component is missing only where its own value is
    excluded <- short(holed, with_wo, na.exclude, dot)
    expect_identical(unname(excluded$stan), unname(reference$stan), info = info)
    ev <- tryCatch(suppressWarnings(extract(excluded, "ev")), error = function(e) NULL)
    expect_identical(dim(ev), c(n, ncol(extract(reference, "ev"))), info = info)
    if (!identical(nrow(ev), n)) next
    expect_true(all(is.na(ev[dropped, ])), info = info)
    expect_equal(ev[kept, ], extract(reference, "ev"), check.attributes = FALSE, info = info)
    expect_identical(predict(excluded), ev, info = info)
    expect_equal(fitted(excluded), replace(rep_len(NA_real_, n), kept, fitted(reference)),
                 info = info)
    expect_equal(residuals(excluded), replace(rep_len(NA_real_, n), kept, residuals(reference)),
                 info = info)
    ppd <- suppressWarnings(extract(excluded, "ppd"))
    expect_identical(which(is.na(ppd[, 1L])), dropped, info = info)
    parts <- list(indiv.fixef = "fixed", indiv.ranef = "slope", indiv.bart = c("forest", "dot"))
    for (type in names(parts)) {
      part <- suppressWarnings(extract(excluded, type))
      expect_identical(nrow(part), n, info = paste(info, type))
      if (!identical(nrow(part), n)) next
      expect_identical(which(is.na(part[, 1L])),
                       if (case %in% parts[[type]] || type == "indiv.bart") dropped else integer(),
                       info = paste(info, type))
      expect_equal(part[kept, ], extract(reference, type), check.attributes = FALSE,
                   info = paste(info, type))
    }

    # na.fail stops
    expect_error(fit_rows(holed, with_wo, na.fail, dot, iter = 0L), "missing values in object",
                 info = info)
  }

  # pieces on different rows never reach the sampler
  fit0 <- fit_rows(d, TRUE, iter = 0L)
  for (piece in c("X", "weights", "offset")) {
    torn <- fit0
    torn[[piece]] <- if (piece == "X") fit0$X[-1L, , drop = FALSE] else fit0[[piece]][-1L]
    expect_error(stan4bart:::stan4bart_fit(torn, torn$family, numeric(n), 1, range(d$y), -1L, 7L,
                                           3L, 1L, 1L, 1L, 2L, 9L, TRUE, FALSE, NULL, NULL, NULL,
                                           list(n.trees = 5L)),
                 "are not on the same rows", fixed = TRUE, info = piece)
  }
})

# the test sample is cut the same way: a missing value in one component of
# the test data drops its row from every component under na.omit
local({
  set.seed(32L)
  n <- 60L
  d <- data.frame(g = factor(sample(6L, n, TRUE)), x = rnorm(n), s = runif(n), u = runif(n),
                  v = runif(n))
  d$y <- 2 * d$x + 3 * d$u + (rnorm(6L) + d$s * rnorm(6L))[d$g] + rnorm(n)
  dropped <- c(3L, 9L)
  for (column in c("u", "s", "x")) {
    test <- d[1:20, ]
    test[[column]][dropped] <- NA
    fit <- stan4bart(y ~ x + (1 + s | g) + bart(u + v), d, test = test, na.action = na.omit,
                     cores = 1L, verbose = -1L, chains = 1L, warmup = 3L, iter = 7L, seed = 9L,
                     bart_args = list(n.trees = 5L))
    reference <- stan4bart(y ~ x + (1 + s | g) + bart(u + v), d, test = d[1:20, ][-dropped, ],
                           cores = 1L, verbose = -1L, chains = 1L, warmup = 3L, iter = 7L,
                           seed = 9L, bart_args = list(n.trees = 5L))
    expect_identical(c(nrow(fit$test$X), ncol(fit$test$reTrms$Zt), nrow(fit$bartData@x.test)),
                     rep(18L, 3L), info = column)
    expect_inherits(fit$test$na.action.all, "omit", info = column)
    ev <- tryCatch(suppressWarnings(extract(fit, "ev", sample = "test")), error = function(e) NULL)
    expect_identical(ev, extract(reference, "ev", sample = "test"), info = column)
  }
})
