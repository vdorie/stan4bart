# extract for trees

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

# stan4bart extracts trees correctly
n.chains <- 1L
n.trees <- 3L
n.samples <- 4L
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = -1L,
                 chains = n.chains,
                 warmup = 0L,
                 iter = n.samples,
                 bart_args = list(n.trees = n.trees, keepTrees = TRUE))

allTrees <- extract(fit, "trees")

expect_true(all(c("sample", "tree") %in% colnames(allTrees)))
expect_true(!("chain" %in% colnames(allTrees)))

combinations <- data.frame(
  sample = rep(seq_len(n.samples), each = n.trees),
  tree   = rep(seq_len(n.trees), times = n.samples)
)
expect_true(all(interaction(combinations$sample, combinations$tree) %in%
                interaction(allTrees$sample, allTrees$tree)))

individualSamples <-
  lapply(seq_len(n.samples), function(i) extract(fit, "trees", sampleNums = i))
individualSamples <- Reduce(rbind, individualSamples)
row.names(individualSamples) <- as.character(seq_len(nrow(individualSamples)))

expect_equal(allTrees, individualSamples)

# several chains: each chain's trees carry its own label, in the order asked for
n.chains <- 3L
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = -1L, seed = 5L,
                 chains = n.chains,
                 warmup = 3L,
                 iter = 3L + n.samples,
                 bart_args = list(n.trees = n.trees, keepTrees = TRUE))

allTrees <- extract(fit, "trees")
expect_identical(colnames(allTrees)[1L], "chain")

# the splits of each chain and sample are the ones the sampler counted while
# it ran; the chains differ, so a mislabeled chain cannot match
varcount <- extract(fit, "varcount", combine_chains = FALSE)
expect_false(identical(varcount[, , 1L], varcount[, , 2L]))
expect_false(identical(varcount[, , 1L], varcount[, , 3L]))
expect_false(identical(varcount[, , 2L], varcount[, , 3L]))
splits <- allTrees[allTrees$var > 0L, ]
splitCounts <- table(factor(splits$var, seq_len(dim(varcount)[1L])),
                     factor(splits$sample, seq_len(n.samples)),
                     factor(splits$chain, seq_len(n.chains)))
expect_identical(as.vector(splitCounts), as.vector(varcount))

# chains come back in the order named, a repeated one repeated
chainNums <- c(3L, 1L, 3L)
someTrees <- do.call(rbind, lapply(chainNums, function(i) allTrees[allTrees$chain == i, ]))
row.names(someTrees) <- as.character(seq_len(nrow(someTrees)))
expect_identical(extract(fit, "trees", chainNums = chainNums), someTrees)

noTrees <- extract(fit, "trees", chainNums = integer(0L))
expect_identical(nrow(noTrees), 0L)
expect_identical(colnames(noTrees), colnames(allTrees))

for (chainNum in c(4L, 0L, -1L))
  expect_error(extract(fit, "trees", chainNums = chainNum), "'chainNums' must be in [1, 3]", fixed = TRUE)
expect_error(extract(fit, "trees", chainNums = NA), "'chainNums' contains missing values", fixed = TRUE)

