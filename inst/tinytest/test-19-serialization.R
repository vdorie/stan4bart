# stan4bart bart-state serialization and lazy pointer rebuild

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(120, TRUE, TRUE, FALSE)
binaryData <- generateFriedmanData(120, TRUE, TRUE, TRUE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
df_b <- with(binaryData, data.frame(x, g.1, g.2, y, z))

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

# keepTrees fit retains serializable bart state and rebuilds the pointer after reload
if (at_home()) local({

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

# keepTrees = FALSE fit retains no bart state (no object-size change)
local({
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

# A serialized external pointer reads back dead in the same session, so the
# rebuild after a reload is reached without a fresh one.
reload <- function(fit) {
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(fit, file)
  readRDS(file)
}
fit_kept <- function(df) {
  stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
            cores = 1L, verbose = -1L, chains = 3L, warmup = 7L, iter = 13L, seed = 5L,
            bart_args = list(n.trees = 10L, keepTrees = TRUE))
}

# a continuous fit's kept trees replay to its stored fits after a reload, each
# chain through a sampler of its own at the response mapping it ended warmup on
local({
  fit <- fit_kept(df)

  # PRECONDITION: no two chains ended warmup on the same mapping, so no one
  # sampler's mapping reads every chain's trees
  mappings <- vapply(fit$state.bart$state, function(state) state[[1L]]$fit.scale, double(2L))
  expect_identical(anyDuplicated(mappings, MARGIN = 2L), 0L)

  stored <- extract(fit, "indiv.bart", combine_chains = FALSE)
  expect_equal(predict(fit, df, type = "indiv.bart", combine_chains = FALSE), stored,
               tolerance = 1e-10, check.attributes = FALSE)

  fit2 <- reload(fit)
  expect_false(stan4bart:::bart_pointer_is_live(fit2$sampler.bart))
  expect_equal(predict(fit2, df, type = "indiv.bart", combine_chains = FALSE), stored,
               tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(predict(fit2, df, type = "ev"), extract(fit, "ev"), tolerance = 1e-10)
  expect_identical(length(stan4bart:::getBartSampler(fit2)), 3L)
  expect_identical(extract(fit2, "trees"), extract(fit, "trees"))

  # each restored sampler reports the mapping its chain recorded, as the
  # width and the midpoint dbarts's getLeafPrior reads off the sampler
  for (restored in list(fit, fit2)) for (i in seq_len(3L)) {
    leafPrior <- stan4bart:::getBartSampler(restored)[[i]]$getLeafPrior()
    expect_identical(leafPrior$response.scale, mappings[2L, i] - mappings[1L, i])
    expect_equal(leafPrior$response.shift, mean(mappings[, i]), tolerance = 1e-12)
  }

  # the re-anchor lands on the recorded mapping only with no offset in force
  restore <- fit$state.bart
  restore$data@offset <- double(nrow(df))
  expect_error(stan4bart:::restoreBartSampler(restore$control, restore$model, restore$data,
                                              restore$state, restore$active, TRUE),
               "carrying an offset", fixed = TRUE)

  # the reloaded trees are tied to the fits the run stored: a leaf's value is
  # its share of the fit on its chain's mapping, so over a draw's leaves the
  # values weighted by their row counts sum to that draw's stored fits
  trees <- extract(fit2, "trees")
  leaves <- trees[trees$var == -1L, ]
  leafSums <- tapply(leaves$n * leaves$value, list(leaves$sample, leaves$chain), sum)
  storedSums <- sweep(sweep(leafSums, 2L, mappings[2L, ] - mappings[1L, ], "*"), 2L,
                      nrow(df) * colMeans(mappings), "+")
  expect_equal(storedSums, apply(stored, c(2L, 3L), sum), tolerance = 1e-10, check.attributes = FALSE)
})

# a binary fit's mapping is fixed, and its kept trees replay to its stored fits
# the same way
local({
  fit <- fit_kept(df_b)

  mappings <- vapply(fit$state.bart$state, function(state) state[[1L]]$fit.scale, double(2L))
  expect_identical(anyDuplicated(mappings, MARGIN = 2L), 2L)

  stored <- extract(fit, "indiv.bart", combine_chains = FALSE)
  expect_equal(predict(fit, df_b, type = "indiv.bart", combine_chains = FALSE), stored,
               tolerance = 1e-10, check.attributes = FALSE)

  fit2 <- reload(fit)
  expect_false(stan4bart:::bart_pointer_is_live(fit2$sampler.bart))
  expect_equal(predict(fit2, df_b, type = "indiv.bart", combine_chains = FALSE), stored,
               tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(predict(fit2, df_b, type = "ev"), extract(fit, "ev"), tolerance = 1e-10)
})
