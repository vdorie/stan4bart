context("extract for trees")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

test_that("stan4bart extracts trees correctly", {
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
})

