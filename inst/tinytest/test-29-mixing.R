# end-of-fit mixing check: split rank-normalized R-hat of sigma, each group sd
# and the in-sample fit

split_rhat <- stan4bart:::split_rhat
fit_mixing <- stan4bart:::fit_mixing
format_mixing <- stan4bart:::format_mixing

# the helper against posterior::rhat (1.7.0), computed once on these matrices;
# built from sin and cos so that no random stream is involved
pin <- function(n, m, f) outer(seq_len(n), seq_len(m), f)
cases <- list(
  agree     = list(pin(40, 4, function(i, j) sin(i * 1.7 + j * 2.3) + cos(i * i * 0.37 + j)), 0.9944477630),
  shifted   = list(pin(41, 3, function(i, j) sin(i * 0.9 + j) + (j == 2) * 0.8), 1.1527850310),
  trend     = list(pin(30, 2, function(i, j) sin(i * 1.3 + j * 5) + 0.05 * i), 1.1479243606),
  one_chain = list(pin(24, 1, function(i, j) sin(i * 2.1) + cos(i * i * 0.11)), 0.9770982846),
  scale     = list(pin(20, 4, function(i, j) (1 + 3 * (j == 4)) * sin(i * 1.9 + j * 7.1)), 1.3763920699),
  heavy     = list(pin(21, 3, function(i, j) exp(sin(i * 3.3 + j)) * ifelse(i %% 7 == 0, 40, 1)), 1.0474392455),
  ties      = list(pin(16, 3, function(i, j) round(sin(i * 1.1 + j) * 3)), 0.9710777520))
for (nm in names(cases))
  expect_equal(split_rhat(cases[[nm]][[1L]]), cases[[nm]][[2L]], tolerance = 1e-8, info = nm)

# posterior returns NA for a constant quantity and for draws it cannot split
expect_true(is.na(split_rhat(matrix(2, 10, 3))))
expect_true(is.na(split_rhat(matrix(c(1, 2, 3), 3, 1))))
expect_true(is.na(split_rhat(cbind(c(1:9, NA), 1:10))))
# two-valued columns have no variance left once folded: NA, never NaN
for (two in list(matrix(rep(c(0, 1), each = 8), 16, 1), cbind(rep(0, 8), rep(1, 8)), matrix(rep(c(0, 1), 4), 8, 2)))
  expect_identical(split_rhat(two), NA_real_)

# the summary's lines: a star on a value over 1.1 and a line saying what it
# means, nothing more when none is over; an NA is said to be unchecked
expect_identical(format_mixing(c(sigma = 1.09, "g sd" = 1.1, fit = 1.0), 4L),
                 "R-hat: sigma 1.09, g sd 1.10, fit 1.00")
expect_identical(format_mixing(c(sigma = 3.3, "g sd" = 1.02, fit = 1.01), 4L),
                 c("R-hat: sigma 3.30*, g sd 1.02, fit 1.01", "* over 1.1: the chains disagree"))
expect_identical(format_mixing(c(sigma = 1.5, "g sd" = 2.5, fit = NA), 1L),
                 c("R-hat: sigma 1.50*, g sd 2.50*, fit NA", "* over 1.1: the halves of the chain disagree",
                   "NA: could not be checked"))
expect_identical(format_mixing(structure(c(sigma = 1.5, fit = NA), reasons = c(fit = "out of memory")), 4L),
                 c("R-hat: sigma 1.50*, fit NA", "* over 1.1: the chains disagree",
                   "NA: could not be checked (out of memory)"))

set.seed(1)
n <- 200
d <- data.frame(g = sample(10, n, TRUE), x1 = rnorm(n), xb = runif(n))
d$y <- 1 + 2 * d$x1 + sin(3 * d$xb) + rnorm(10)[d$g] * 0.7 + rnorm(n) * 0.5
d$z <- as.integer(d$y > 1.5)
# without lme4 a fit with random effects warns that the BART response range or
# the binary start is the plain one; that warning is not under test here
fit_stan4bart <- function(...)
  withCallingHandlers(stan4bart(...), warning = function(cond) {
    if (startsWith(conditionMessage(cond), "lme4 is not installed"))
      invokeRestart("muffleWarning")
  })

args <- list(data = d, cores = 2L, verbose = -1L, seed = 3L, iter = 300L, warmup = 150L,
             bart_args = list(n.trees = 25L))

# a healthy fit is silent, carries the three kinds of quantity, and summary()
# prints them on one line, with no flag
fit <- NULL
expect_silent(fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args)))
expect_equal(names(fit$mixing), c("sigma", "g sd", "fit"))
expect_true(all(fit$mixing < 1.1))
out <- capture.output(print(summary(fit)))
line <- grep("^R-hat: ", out, value = TRUE)
expect_equal(length(line), 1L)
expect_equal(line, format_mixing(fit$mixing, 4L)[1L])
expect_false(any(grepl("over 1.1", out)))
expect_false(any(grepl("R-hat", capture.output(print(fit)))))

# a fit with no grouping factor has no group sd; a binary fit has no sigma
fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb)), args))
expect_equal(names(fit$mixing), c("sigma", "fit"))
fit <- do.call(fit_stan4bart, c(list(z ~ x1 + bart(xb) + (1 | g)), args))
expect_equal(names(fit$mixing), c("g sd", "fit"))

# one chain is checked too, by comparing its halves
fit <- NULL
expect_silent(fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args[names(args) != "cores"],
                                         list(chains = 1L))))
expect_true(all(is.finite(fit$mixing)))

# warmup draws are never read: the check is the same with them saved
fit_w <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args[names(args) != "cores"],
                              list(chains = 1L, save_warmup = TRUE)))
expect_equal(fit_w$mixing, fit$mixing)

# a planted disagreement is flagged by summary(): the same fit with one chain's
# sigma draws moved apart
fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args))
sigma_row <- match("aux.1", dimnames(fit$stan)[[1L]])
fit$stan[sigma_row, , 1L] <- fit$stan[sigma_row, , 1L] + 5
fit$mixing <- fit_mixing(fit)
expect_true(fit$mixing[["sigma"]] > 1.1)
expect_true(all(fit$mixing[c("g sd", "fit")] < fit$mixing[["sigma"]]))
out <- capture.output(print(summary(fit)))
line <- grep("^R-hat: ", out, value = TRUE)
expect_true(grepl("sigma [0-9.]+\\*,", line))
expect_false(grepl("g sd [0-9.]+\\*", line))
expect_true("* over 1.1: the chains disagree" %in% out)
# print() of the fit says nothing of it
expect_false(any(grepl("R-hat|over 1.1", capture.output(print(fit)))))

# a real fit too short to mix raises no warning, whatever verbose says; its
# R-hats are stored, and its summary flags them
for (verbose in c(-1L, 0L)) {
  warnings <- character()
  short <- withCallingHandlers(
    do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)),
                         modifyList(args, list(verbose = verbose, seed = 3L, iter = 12L, warmup = 2L,
                                               bart_args = list(n.trees = 5L))))),
    warning = function(cond) {
      warnings <<- c(warnings, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_equal(length(warnings), 0L, info = verbose)
  expect_true(any(short$mixing > 1.1), info = verbose)
}
out <- capture.output(print(summary(short)))
expect_true("* over 1.1: the chains disagree" %in% out)

# the fit quantity is summed in blocks of rows: the per-draw RMSE is the
# unblocked computation on the whole block of expected values to rounding, and
# so is its R-hat, for a stored block and for kept trees, with and without
# parametric parts, and for a binary response (probability scale). Blocked
# sums, and even one block on some builds (a sanitizer build's), differ from
# the unblocked computation in the last bits, which the R-hat's ranks of the
# draws about their median can amplify to a few parts in a thousand where two
# draws sit equally far from it.
unblocked_rmse <- function(fit) {
  ev <- extract(fit, "ev", combine_chains = FALSE)
  d <- dim(ev)
  resid <- ev - fit$y
  dim(resid) <- c(d[1L], d[2L] * d[3L])
  sqrt(colMeans(resid^2))
}
for (store in c("fits", "trees")) {
  for (form in list(y ~ x1 + bart(xb) + (1 | g), y ~ x1 + bart(xb), z ~ x1 + bart(xb) + (1 | g))) {
    fit <- do.call(fit_stan4bart, c(list(form), args[names(args) != "cores"],
                                    list(chains = 2L, store = store)))
    expected <- unblocked_rmse(fit)
    d <- dim(fit$stan)[2L:3L]
    expected_rhat <- split_rhat(matrix(expected, d[1L], d[2L]))
    expect_equal(fit$mixing[["fit"]], expected_rhat, tolerance = 5e-3, info = store)
    # one block, then blocks of seven rows and of one, so that each sum runs over many
    for (cells in c(2e6, 7 * prod(d), prod(d))) {
      expect_equal(stan4bart:::fit_rmse_draws(fit, cells), expected, tolerance = 1e-12, info = store)
      expect_equal(fit_mixing(fit, block_cells = cells)[["fit"]], expected_rhat, tolerance = 5e-3,
                   info = store)
    }
  }
}

# a failure computing the fit quantity leaves it unchecked, NA and said so in
# summary(), and the other quantities stand
fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args))
broken <- fit
broken$reTrms$Zt <- NULL
mixing <- fit_mixing(broken)
expect_equal(names(mixing), c("sigma", "g sd", "fit"))
expect_true(is.na(mixing[["fit"]]))
expect_equal(mixing[c("sigma", "g sd")], fit$mixing[c("sigma", "g sd")])
broken$mixing <- mixing
out <- capture.output(print(summary(broken)))
expect_true(grepl("fit NA", grep("^R-hat: ", out, value = TRUE)))
reason <- attr(mixing, "reasons")[["fit"]]
expect_true(nzchar(reason))
expect_true(paste0("NA: could not be checked (", reason, ")") %in% out)
# the reason is the failure's own message
expect_identical(reason, tryCatch(stan4bart:::fit_rmse_draws(broken, 2e6), error = conditionMessage))

# keep_fits = FALSE keeps no draws: nothing is checked and nothing is stored
fit <- do.call(fit_stan4bart, c(list(y ~ x1 + bart(xb) + (1 | g)), args, list(keep_fits = FALSE)))
expect_null(fit$mixing)
