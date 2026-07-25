# stan4bart print/summary and the mis-tuning diagnostic

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, FALSE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y))
rm(testData)

# testthat's expect_message(expr, NA) - "emits no message at all" - has no
# tinytest counterpart: collect the messages an expression emits and assert
# on that. Stdout is swallowed, since print() legitimately writes there.
messagesOf <- function(expr) {
  msgs <- character()
  withCallingHandlers(
    invisible(capture.output(expr)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  msgs
}

# default warmup ratio (iter %/% 2), comfortably above the documented floor,
# so this fit's mean_leapfrog should sit well under the warning threshold;
# seeded for a reproducible (non-flaky) diagnostic value
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4) + X4 + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = -1L, chains = 2, iter = 300,
                 bart_args = list(n.trees = 10), seed = 7)

# a well-warmed fit's print/summary do not warn
expect_true(max(fit$adaptation$mean_leapfrog) < 30)
expect_equal(length(messagesOf(print(fit))), 0L)
expect_equal(length(messagesOf(summary(fit))), 0L)

# print and summary warn when mean_leapfrog is inflated
local({
  bad_fit <- fit
  bad_fit$adaptation$mean_leapfrog[1] <- 50
  msg <- "parametric sampler tuning looks poor"
  expect_true(any(grepl(msg, messagesOf(print(bad_fit)))))
  expect_true(any(grepl(msg, messagesOf(summary(bad_fit)))))
})

# a fit with no adaptation record does not error or warn
local({
  no_adapt_fit <- fit
  no_adapt_fit$adaptation <- NULL
  expect_equal(length(messagesOf(print(no_adapt_fit))), 0L)
  expect_equal(length(messagesOf(summary(no_adapt_fit))), 0L)
})

# summary returns a summary.stan4bartFit carrying the per-chain diagnostics
local({
  s <- summary(fit)
  expect_true(inherits(s, "summary.stan4bartFit"))
  expect_equal(s$adaptation$mean_leapfrog, fit$adaptation$mean_leapfrog)
  expect_equal(s$adaptation$mean_leapfrog_warmup, fit$adaptation$mean_leapfrog_warmup)
})

# fitted(type = 'ppd') nudges toward 'ev' once per session
rm(list = ls(stan4bart:::.message_env), envir = stan4bart:::.message_env)
expect_true(any(grepl("computes it exactly and faster",
                      messagesOf(invisible(fitted(fit, "ppd"))))))
expect_equal(length(messagesOf(invisible(fitted(fit, "ppd")))), 0L)
