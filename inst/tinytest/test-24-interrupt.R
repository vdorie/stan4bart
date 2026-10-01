# An interrupt during a fit, and the slow-count warning of a monotone BART
# component. dbarts' run polls for R's interrupt and turns it into an error
# that unwinds through stan4bart's sweep loop; the loop's draws land in R
# storage, so nothing is stranded and a later fit is unaffected. The interrupt
# is injected through dbarts' internal count hooks, which report an interrupt
# on the Nth poll without touching R's signal state. They are an internal
# entry, so this file skips when the installed dbarts' entry takes another
# number of arguments than the one it was written against.

hooks <- tryCatch(
  get("C_dbarts_bartcore_setMonotoneCountHooks", asNamespace("dbarts")),
  error = function(e) NULL
)
if (is.null(hooks) || !identical(hooks$numParameters, 3L)) {
  exit_file("dbarts' count hooks are not the ones this file drives")
}
# returns the slow-count threshold it replaced, so a caller can put it back
countHooks <- function(slowSeconds = NA_real_, interruptAfterPolls = NA_integer_) {
  invisible(.Call(hooks, as.double(slowSeconds), FALSE,
                  as.integer(interruptAfterPolls)))
}

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)
testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)
df <- with(testData, data.frame(x, g.1, g.2, y, z))

fitOnce <- function(...) {
  stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
            cores = 1, verbose = -1, chains = 1, warmup = 10, iter = 20, seed = 2,
            ...)
}

before <- fitOnce()

# the fifth poll falls inside the warmup loop: one at creation's first draw,
# then one per sweep. The hook is process-wide, so it is disarmed whatever
# the fit did.
interrupted <- local({
  on.exit(countHooks(interruptAfterPolls = 0L))
  countHooks(interruptAfterPolls = 5L)
  tryCatch(fitOnce(), error = conditionMessage)
})
expect_true(is.character(interrupted) &&
            grepl("sampler run interrupted", interrupted))

# a fresh fit afterwards is the fit it would have been
after <- fitOnce()
expect_equal(after$bart_train, before$bart_train)
expect_equal(after$stan, before$stan)

# a monotone component whose leaf-order counts are slow warns once per fit:
# warmup and sampling run on one sampler, and dbarts warns once per sampler
slowWarnings <- local({
  previous <- countHooks(slowSeconds = -1)
  on.exit(countHooks(slowSeconds = previous))
  count <- 0L
  withCallingHandlers(
    fitOnce(bart_args = list(
      n.trees = 4,
      monotone = dbarts:::monotone(c(X1 = "increasing"), prior = "leaf"))),
    warning = function(w) {
      if (inherits(w, "dbartsSlowCountWarning")) count <<- count + 1L
      invokeRestart("muffleWarning")
    })
  count
})
expect_equal(slowWarnings, 1L)
