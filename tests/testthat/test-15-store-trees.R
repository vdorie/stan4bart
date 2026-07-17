context("stan4bart store = 'trees' recompute seam (bart-train-recompute C3)")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

# gaussian and binary tiers, same causal+ranef structure as the other suites.
# g.1/g.2 as factors so test data carries all levels; test = df (whole training
# frame) so the bart_test recompute path is exercised with every level present.
gaussianData <- generateFriedmanData(400, TRUE, TRUE, FALSE)
binaryData   <- generateFriedmanData(400, TRUE, TRUE, TRUE)
rm(generateFriedmanData)

df_g <- with(gaussianData, data.frame(x, g.1 = as.factor(g.1), g.2 = as.factor(g.2), y, z))
df_b <- with(binaryData,   data.frame(x, g.1 = as.factor(g.1), g.2 = as.factor(g.2), y, z))

# A bart()-column factor (categorical splits), unlike df_g/df_b above where
# g.1/g.2 are factors but excluded from bart() and only used for the random
# effect grouping. bartData@x is then a dbartsMixedMatrix, not a plain matrix -
# the recompute seam (recompute_bart_block, generics.R) must densify it the
# same way makeTestModelMatrix already densifies bartData@x.test.
df_factor <- df_g
df_factor$X1 <- cut(df_factor$X1, quantile(df_factor$X1, seq(0, 1, 0.2)), include.lowest = TRUE)
levels(df_factor$X1) <- letters[seq_len(nlevels(df_factor$X1))]

# Formula must be written out at the call site: bart()/bar-term detection walks
# the UNEVALUATED formula from match.call(). store is an ordinary (evaluated)
# argument, so it can be threaded through a variable. Both stores are fit with
# keepTrees = TRUE so the ONLY difference is R-side storage - sampling, and hence
# every parametric draw, is bit-identical.
fit_store <- function(store_mode, df) {
  stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
            data = df, test = df,
            cores = 1L, verbose = -1L, chains = 2L,
            warmup = 5L, iter = 30L, seed = 99L,
            store = store_mode,
            bart_args = list(n.trees = 15L, keepTrees = TRUE))
}

# Fresh-session runner copied from test-14: write a script inheriting this
# session's libpaths, run under Rscript, read the saved result back. deparse()
# wraps long vectors, so collapse the .libPaths line to one expression.
run_in_fresh_session <- function(body_lines) {
  script <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               "suppressMessages(library(stan4bart))",
               body_lines),
             script)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0L, output = out)
}

TOL <- 1e-10   # matches the C0 exactness gate; NEVER identical()-based

test_that("store = 'trees' drops the bart blocks but keeps the parametric draws bit-identical", {
  skip_on_cran()

  fit_fits  <- fit_store("fits",  df_g)
  fit_trees <- fit_store("trees", df_g)

  # opt-in absence / default presence of the stored blocks
  expect_false(is.null(fit_fits$bart_train));  expect_false(is.null(fit_fits$bart_test))
  expect_true(is.null(fit_trees$bart_train));  expect_true(is.null(fit_trees$bart_test))
  expect_true(is.null(fit_trees$warmup$bart_train))
  # the trees are retained so recompute can rebuild the sampler after reload
  expect_false(is.null(fit_trees$state.bart))
  expect_false(is.null(fit_trees$sampler.bart))
  # everything else that is NOT the n x draws block is still stored
  expect_false(is.null(fit_trees$bart_varcount))
  expect_false(is.null(fit_trees$stan))

  # PRECONDITION: storage does not touch sampling, so with the same seed the
  # parametric draws and variable counts are bit-identical across the two stores
  expect_identical(fit_trees$stan,          fit_fits$stan)
  expect_identical(fit_trees$bart_varcount, fit_fits$bart_varcount)
})

test_that("recomputed BART surfaces match the stored path across readers (gaussian)", {
  skip_on_cran()

  fit_fits  <- fit_store("fits",  df_g)
  fit_trees <- fit_store("trees", df_g)

  # extract() materializes the full block via recompute; matches stored to TOL
  for (smp in c("train", "test")) {
    expect_equal(extract(fit_trees, "indiv.bart", sample = smp),
                 extract(fit_fits,  "indiv.bart", sample = smp), tolerance = TOL)
    expect_equal(extract(fit_trees, "ev", sample = smp),
                 extract(fit_fits,  "ev", sample = smp), tolerance = TOL)
  }

  # fitted() streams the posterior mean in row-blocks; matches stored to TOL
  expect_equal(fitted(fit_trees, "indiv.bart"), fitted(fit_fits, "indiv.bart"), tolerance = TOL)
  expect_equal(fitted(fit_trees, "ev"),         fitted(fit_fits, "ev"),         tolerance = TOL)
  expect_equal(fitted(fit_trees, "ev", sample = "test"),
               fitted(fit_fits,  "ev", sample = "test"), tolerance = TOL)

  # predict(no newdata) routes through extract() -> recompute
  expect_equal(predict(fit_trees, type = "ev"), predict(fit_fits, type = "ev"), tolerance = TOL)
  # predict(newdata): both replay the SAME trees over the SAME design -> tighter
  expect_equal(predict(fit_trees, df_g, type = "indiv.bart"),
               predict(fit_fits,  df_g, type = "indiv.bart"), tolerance = 1e-12)

  # record the memory/latency trade at the test scale (plan C3: object-size
  # ratio and one extract latency number)
  size_ratio <- as.numeric(object.size(fit_trees)) / as.numeric(object.size(fit_fits))
  t_extract  <- system.time(invisible(extract(fit_trees, "indiv.bart")))[["elapsed"]]
  message(sprintf("[store=trees] object.size ratio trees/fits = %.3f; extract(indiv.bart) = %.3f s",
                  size_ratio, t_extract))
  expect_true(is.finite(size_ratio))
})

test_that("recomputed BART surfaces match the stored path (binary)", {
  skip_on_cran()

  fit_fits  <- fit_store("fits",  df_b)
  fit_trees <- fit_store("trees", df_b)

  expect_equal(fit_trees$family$family, "binomial")
  expect_identical(fit_trees$stan, fit_fits$stan)

  # ev is pnorm(bart + fixef + ranef): fitted must transform per-draw before
  # averaging, which the streamed path does block-by-block
  expect_equal(extract(fit_trees, "indiv.bart"), extract(fit_fits, "indiv.bart"), tolerance = TOL)
  expect_equal(extract(fit_trees, "ev"),         extract(fit_fits, "ev"),         tolerance = TOL)
  expect_equal(fitted(fit_trees, "ev"),          fitted(fit_fits, "ev"),          tolerance = TOL)
  expect_equal(fitted(fit_trees, "indiv.bart"),  fitted(fit_fits, "indiv.bart"),  tolerance = TOL)
})

test_that("recomputed BART surfaces match the stored path for a factor-bearing bart component", {
  skip_on_cran()

  fit_fits  <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) + X4 + z +
                               (1 + X4 | g.1) + (1 | g.2),
                        data = df_factor, test = df_factor,
                        cores = 1L, verbose = -1L, chains = 2L,
                        warmup = 5L, iter = 30L, seed = 99L,
                        store = "fits",
                        bart_args = list(n.trees = 15L, keepTrees = TRUE))
  fit_trees <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) + X4 + z +
                               (1 + X4 | g.1) + (1 | g.2),
                        data = df_factor, test = df_factor,
                        cores = 1L, verbose = -1L, chains = 2L,
                        warmup = 5L, iter = 30L, seed = 99L,
                        store = "trees",
                        bart_args = list(n.trees = 15L, keepTrees = TRUE))

  expect_true("X1" %in% rownames(fit_trees$bart_varcount))
  expect_identical(fit_trees$stan,          fit_fits$stan)
  expect_identical(fit_trees$bart_varcount, fit_fits$bart_varcount)

  for (smp in c("train", "test")) {
    expect_equal(extract(fit_trees, "indiv.bart", sample = smp),
                 extract(fit_fits,  "indiv.bart", sample = smp), tolerance = TOL)
    expect_equal(extract(fit_trees, "ev", sample = smp),
                 extract(fit_fits,  "ev", sample = smp), tolerance = TOL)
  }

  expect_equal(fitted(fit_trees, "indiv.bart"), fitted(fit_fits, "indiv.bart"), tolerance = TOL)
  expect_equal(fitted(fit_trees, "ev"),         fitted(fit_fits, "ev"),         tolerance = TOL)
})

test_that("a store = 'trees' fit survives saveRDS/readRDS and fitted() works after reload", {
  skip_on_cran()

  fit_trees <- fit_store("trees", df_g)
  fitted_before <- fitted(fit_trees, "ev")

  fit_rds <- tempfile(fileext = ".rds")
  df_rds  <- tempfile(fileext = ".rds")
  out_rds <- tempfile(fileext = ".rds")
  saveRDS(fit_trees, fit_rds)
  saveRDS(df_g, df_rds)

  res <- run_in_fresh_session(c(
    sprintf("fit2 <- readRDS(%s)", deparse(fit_rds)),
    sprintf("df2  <- readRDS(%s)", deparse(df_rds)),
    # the live pointer is dead after reload; recompute must rebuild it lazily
    "dead <- !stan4bart:::bart_pointer_is_live(fit2$sampler.bart)",
    "stopifnot(is.null(fit2$bart_train))",
    "fitted2 <- fitted(fit2, 'ev')",
    "cached  <- stan4bart:::bart_pointer_is_live(fit2$bart_env$ptr)",
    sprintf("saveRDS(list(dead = dead, cached = cached, fitted = fitted2), %s)", deparse(out_rds))))

  expect_true(file.exists(out_rds), info = paste(res$output, collapse = "\n"))
  reloaded <- readRDS(out_rds)
  expect_true(reloaded$dead)
  expect_true(reloaded$cached)
  expect_equal(reloaded$fitted, fitted_before, tolerance = TOL)
})

test_that("store = 'trees' refuses a contradictory keepTrees = FALSE", {
  skip_on_cran()

  expect_error(
    stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
              data = df_g, cores = 1L, verbose = -1L, chains = 1L,
              warmup = 2L, iter = 6L, seed = 99L, store = "trees",
              bart_args = list(n.trees = 10L, keepTrees = FALSE)),
    "requires keepTrees")
})

test_that("a large store = 'fits' fit nudges toward store = 'trees' once per session", {
  # helper-level test: a fit large enough to trip the gigabyte threshold is
  # not buildable in a test suite, so exercise the nudge on the byte count it
  # would receive from package_samples
  rm(list = ls(stan4bart:::.message_env), envir = stan4bart:::.message_env)
  expect_message(stan4bart:::nudge_if_bart_store_large(2^31), "store = \"trees\"")
  expect_message(stan4bart:::nudge_if_bart_store_large(2^31), NA)
  rm(list = ls(stan4bart:::.message_env), envir = stan4bart:::.message_env)
  expect_message(stan4bart:::nudge_if_bart_store_large(2^20), NA)
})
