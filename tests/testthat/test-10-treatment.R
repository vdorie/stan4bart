context("stan4bart treatment argument")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))

test_that("works with logical treatment", {
  df$z <- df$z > 0

  fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   treatment = z)
  expect_is(fit, "stan4bartFit")

  expect_is(fitted(fit, sample = "test"), "numeric")
})

