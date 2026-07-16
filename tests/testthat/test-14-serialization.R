context("stan4bart bart-state serialization and lazy pointer rebuild")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(120, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

# Run an expression in a FRESH R session (no callr dependency): write a script
# that inherits this session's library paths, run it with Rscript, and read
# back its saved result. Mirrors the standing reload bug's reproduction - the
# fit's live externalptr dies across the process boundary.
run_in_fresh_session <- function(body_lines) {
  script <- tempfile(fileext = ".R")
  # deparse() wraps long vectors across elements; collapse so the generated
  # line stays one valid expression
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               "suppressMessages(library(stan4bart))",
               body_lines),
             script)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0L, output = out)
}

test_that("keepTrees fit retains serializable bart state and rebuilds the pointer after reload", {
  skip_on_cran()

  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                   cores = 1L, verbose = -1L, chains = 2L,
                   warmup = 2L, iter = 8L,
                   bart_args = list(n.trees = 10L, keepTrees = TRUE))

  # C2: the serializable inputs are retained and a session cache cell exists
  expect_false(is.null(fit$state.bart))
  expect_true(is.environment(fit$bart_env))

  # in-session path is UNCHANGED: getBartSampler returns the original live
  # pointer with no rebuild (the cache cell stays empty)
  expect_identical(stan4bart:::getBartSampler(fit), fit$sampler.bart)
  expect_null(fit$bart_env$ptr)

  # in-session reference values to reproduce after reload
  pred_before  <- predict(fit, df, type = "indiv.bart", combine_chains = TRUE)
  trees_before <- extract(fit, "trees")

  fit_rds   <- tempfile(fileext = ".rds")
  df_rds    <- tempfile(fileext = ".rds")
  pred_rds  <- tempfile(fileext = ".rds")
  trees_rds <- tempfile(fileext = ".rds")
  out_rds   <- tempfile(fileext = ".rds")
  saveRDS(fit, fit_rds)
  saveRDS(df,  df_rds)

  res <- run_in_fresh_session(c(
    sprintf("fit2 <- readRDS(%s)", deparse(fit_rds)),
    sprintf("df2  <- readRDS(%s)", deparse(df_rds)),
    # the live pointer arrives DEAD after reload (this is the standing bug)
    "dead <- !stan4bart:::bart_pointer_is_live(fit2$sampler.bart)",
    # predict/extract for the bart component succeed via the lazy rebuild
    "pred2  <- predict(fit2, df2, type = 'indiv.bart', combine_chains = TRUE)",
    "trees2 <- extract(fit2, 'trees')",
    # rebuilt once and now held in the session cache cell
    "cached <- stan4bart:::bart_pointer_is_live(fit2$bart_env$ptr)",
    sprintf("saveRDS(list(dead = dead, cached = cached), %s)", deparse(out_rds)),
    sprintf("saveRDS(pred2,  %s)", deparse(pred_rds)),
    sprintf("saveRDS(trees2, %s)", deparse(trees_rds))))

  expect_true(file.exists(out_rds), info = paste(res$output, collapse = "\n"))
  flags <- readRDS(out_rds)

  # the reload genuinely killed the live pointer, and the rebuild recovered it
  expect_true(flags$dead)
  expect_true(flags$cached)

  # reloaded predict/extract match the in-session values to tight tolerance
  expect_equal(readRDS(pred_rds),  pred_before, tolerance = 1e-12)
  expect_equal(readRDS(trees_rds), trees_before)
})

test_that("keepTrees = FALSE fit retains no bart state (no object-size change)", {
  fit0 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                    cores = 1L, verbose = -1L, chains = 1L,
                    warmup = 2L, iter = 8L,
                    bart_args = list(n.trees = 10L))

  # nothing bart-state related is attached, so the object is byte-for-byte the
  # size it was before C2 (the retained state lands ONLY under keepTrees)
  expect_null(fit0$sampler.bart)
  expect_null(fit0$state.bart)
  expect_null(fit0$bart_env)
})
