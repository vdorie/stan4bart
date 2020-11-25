
context("mstan4bart continuous response")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriendmanData(100, TRUE, TRUE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

test_that("extract matches fitted in causal setting model", {
  df.train <- df[seq_len(floor(0.8 * nrow(df))),]
  df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]
  
  fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df.train,
                    cores = 2, verbose = 0, chains = 3, warmup = 7, iter = 13,
                    bart_args = list(n.trees = 11),
                    treatment = z)
  
  samples.pred <- predict(fit, df.train)
  
  samples.ev.train <- extract(fit)
  samples.ev.test  <- extract(fit, sample = "test")
  
  samples.icate <- (samples.ev.train - samples.ev.test) * (2 * testData$z - 1)
  samples.cate <- apply(samples.icate, 2, mean)
  cate <- mean(samples.cate)
  
  samples.ppd.test <- extract(fit, type = "ppd", sample = "test")
  
  samples.ite <- (df.train$y - samples.ppd.test) * (2 * df.train$z - 1)
  samples.sate <- apply(samples.ite, 2, mean)
  sate <- mean(samples.sate)
  
  samples.ppd.test <- extract(fit, type = "ppd", sample = "train")
  samples.pate <- apply((samples.ppd.test - samples.ppd.test) * (2 * df.train$z - 1), 2, mean)
  pate <- mean(samples.pate)
  
  fitted.ev.train <- apply(samples.ev.train, 1, mean)
  fitted.ev.test  <- apply(samples.ev.test,  1, mean)
  
  expect_equal(fitted.ev.train, fitted(fit, "ev", "train"))
  expect_equal(fitted.ev.train, fitted(fit, "ev", "test"))
})
