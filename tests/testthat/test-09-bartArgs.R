context("extract for trees")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

test_that("bart args are set correctly", {
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
})

