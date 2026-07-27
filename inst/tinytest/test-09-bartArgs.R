# extract for trees

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

# bart args are set correctly
output <- character()
messages <- character()
outputConnection <- textConnection("output", "w", local = TRUE)
messagesConnection <- textConnection("messages", "w", local = TRUE)
sink(outputConnection)
sink(messagesConnection, type = "message")

fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = 2, chains = 1, warmup = 0, iter = 1,
                 bart_args = list(n.trees = 2, power = 2.5, base = 0.9, split.probs = c(X3 = 2, .default = 1)),
                 treatment = z)

sink(type = "message")
sink()

powerAndBaseLine <- grep("\tpower and base for tree prior:", output)
expect_true(length(powerAndBaseLine) > 0L)

powerAndBase <-
  as.numeric(strsplit(sub("^[^0-9]+([0-9.]+ [0-9.]+)$", "\\1", output[powerAndBaseLine]), " ")[[1L]])

expect_equal(powerAndBase[1L], 2.5)
expect_equal(powerAndBase[2L], 0.9)

splitProbsLine <- grep("\ttree split probabilities:", output)
expect_true(length(splitProbsLine) > 0L)
splitProbs <-
  as.numeric(strsplit(sub("^[^0-9]+((?:[0-9.]+, )*[0-9.]+)$", "\\1", output[splitProbsLine], perl = TRUE), ", ")[[1L]])
expect_equal(splitProbs, c(0.1, 0.1, 0.2, 0.1, 0.1))


# dbartsSpec's own arguments ride through bart_args ---------------------------

fitWith <- function(bart_args) {
  set.seed(21)
  stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
            cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
            bart_args = bart_args)
}

# The priors themselves are reachable, and the k/power/base shorthands are
# exactly the node.prior/tree.prior they write into. dbarts's prior vocabulary
# is resolved by dbarts, so a prior has to be spelled inline in the stan4bart
# call - the same requirement 'chi' has always had.
set.seed(21)
fit.k <- stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
                   cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                   bart_args = list(n.trees = 3, k = 3))
set.seed(21)
fit.node <- stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
                      cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                      bart_args = list(n.trees = 3, node.prior = normal(k = 3)))
expect_equal(fit.k$bart_train, fit.node$bart_train)

set.seed(21)
fit.pb <- stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
                    cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                    bart_args = list(n.trees = 3, power = 2.5, base = 0.9))
set.seed(21)
fit.tree <- stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
                      cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
                      bart_args = list(n.trees = 3, tree.prior = cgm(power = 2.5, base = 0.9)))
expect_equal(fit.pb$bart_train, fit.tree$bart_train)

# a shorthand and the prior it writes into cannot both be given
expect_error(
  stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
            cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
            bart_args = list(k = 3, node.prior = normal(k = 2))),
  "both 'k' and 'node.prior'")
expect_error(
  stan4bart(y ~ bart(X1 + X2 + X3) + X4 + z + (1 | g.1), df,
            cores = 1, verbose = -1L, chains = 1, warmup = 3, iter = 6,
            bart_args = list(power = 2.5, tree.prior = cgm())),
  "both 'tree.prior' and 'power'")

# model-level constraints reach dbarts, which is what validates them
expect_true(inherits(fitWith(list(n.trees = 3, monotone = c(X1 = 1))), "stan4bartFit"))
expect_error(fitWith(list(n.trees = 3, monotone = c(X1 = 1),
                          proposal.probs = c(birth_death = 0.6, swap = 0.1,
                                             change = 0.3, birth = 0.5))),
             "birth/death-only proposals")
expect_error(fitWith(list(n.trees = 3, monotone = c(nosuchvariable = 1))),
             "unrecognized variable name")
expect_error(fitWith(list(n.trees = 3, interactions = dbarts::interactions(max.order = 0L))),
             "max.order")

# the residual model belongs to the parametric component, not the forest
for (reserved in c("sigma", "resid.prior", "resid.dist", "variance"))
  expect_error(fitWith(structure(list(3, 1), names = c("n.trees", reserved))),
               paste0("cannot set '", reserved, "'"))

# the response family follows the response, not bart_args
expect_error(fitWith(list(n.trees = 3, family = "logistic")),
             "must be \"gaussian\"")
