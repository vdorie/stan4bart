context("stan4bart na.action")

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

df.train <- df[seq_len(floor(0.8 * nrow(df))),]
df.test  <- df[seq.int(floor(0.8 * nrow(df)) + 1L, nrow(df)),]

df.train$X1[1L:3L] <- NA
df.train$X4[4L] <- NA
df.train$g.1[5L:6L] <- NA

df.test$X1[1L:3L] <- NA
df.test$X4[4L] <- NA
df.test$g.1[5L:6L] <- NA


test_that("na.action = na.omit works correctly", {
  fit <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) + 
                         X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   df.train, test = df.test,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   na.action = na.omit)
  
  types <- c("ev", "indiv.fixef", "indiv.ranef", "indiv.bart")
  for (type in types) {
    value <- fitted(fit, type = type)
    
    expect_true(length(value) == nrow(df.train) - 6L)
    expect_true(!anyNA(value))
    
    value <- fitted(fit, type = type, sample = "test")
    
    expect_true(length(value) == nrow(df.test) - 6L)
    expect_true(!anyNA(value))
  }
})

test_that("na.action = na.exclude works correctly", {
  fit <- stan4bart(y ~ bart(X1 + X2 + X3 + X5 + X6 + X7 + X8 + X9 + X10) + 
                         X4 + z + (1 + X4 | g.1) + (1 | g.2),
                   df.train, test = df.test,
                   cores = 1, verbose = -1L, chains = 3, warmup = 7, iter = 13,
                   bart_args = list(n.trees = 11),
                   na.action = na.exclude)

  value <- fitted(fit, type = "indiv.ranef")
  expect_true(length(value) == nrow(df.train))
  # 3 NA: two missing group and one missing group predictor
  expect_true(sum(is.na(value)) == 3L) 
  
  # This should complain as we had to fit with data minus 6 rows and 
  # adding the three ones for which predictions can be made to the
  # test data is too much work.
  expect_warning(value <- fitted(fit, type = "indiv.bart"))
  expect_true(length(value) == nrow(df.train))
  expect_true(sum(is.na(value)) == 6L)
  
  value <- fitted(fit, type = "indiv.fixef")
  expect_true(length(value) == nrow(df.train))
  expect_true(sum(is.na(value)) == 1L) 
  
  value <- fitted(fit)
  expect_true(length(value) == nrow(df.train))
  expect_true(sum(is.na(value)) == 6L)
 
  
  value <- fitted(fit, type = "indiv.ranef", sample = "test")
  expect_true(length(value) == nrow(df.test))
  expect_true(sum(is.na(value)) == 3L) 
  
  value <- fitted(fit, type = "indiv.bart", sample = "test")
  expect_true(length(value) == nrow(df.test))
  # With the test data, we can give it a full data frame containing
  # all possible rows, whether or not we would have to drop them
  # elsewhere.
  expect_true(sum(is.na(value)) == 3L)
  
  value <- fitted(fit, type = "indiv.fixef", sample = "test")
  expect_true(length(value) == nrow(df.test))
  expect_true(sum(is.na(value)) == 1L) 
  
  value <- fitted(fit, sample = "test")
  expect_true(length(value) == nrow(df.test))
  expect_true(sum(is.na(value)) == 6L)
})


