context("mstan4bart factor levels")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

df$X1 <- cut(df$X1, quantile(df$X1, seq(0.0, 1.0, 0.2)), include.lowest = TRUE)
levels(df$X1) <- letters[seq_len(nlevels(df$X1))]

df.train <- df[seq_len(floor(0.8 * nrow(df))),]
df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]

levels(df.train$X1) <- c(levels(df.train$X1), letters[1L:2L + nlevels(df.train$X1)])
levels(df.test$X1) <- c(levels(df.test$X1), letters[c(1L, 3L) + nlevels(df.test$X1)])

test_that("model fits with empty factor levels", {
  fit <- mstan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) + 
                             X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df.train, test = df.test,
                  cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11))
  
  types <- c("ev", "indiv.fixef", "indiv.ranef", "indiv.bart")
  for (type in types) {
    value <- fitted(fit, type = type)
    
    expect_true(length(value) == nrow(df.train))
    expect_true(!anyNA(value))
    
    value <- fitted(fit, type = type, sample = "test")
    
    expect_true(length(value) == nrow(df.test))
    expect_true(!anyNA(value))
  }
})

