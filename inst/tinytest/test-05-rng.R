# stan4bart continuous response

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(100, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

# stan4bart with global seed and one thread is reproducible
set.seed(12345L)
fit1 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  chains = 2, cores = 1)


fit2 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  seed = 12345L,
                  chains = 2, cores = 1)
expect_equal(fit1$bart_train, fit2$bart_train)

# stan4bart with fixed seed is reproducible
fit1 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  seed = 12345L,
                  chains = 2, cores = 1)

fit2 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  seed = 12345L,
                  chains = 2, cores = 1)

expect_equal(fit1$bart_train, fit2$bart_train)

fit3 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  seed = 12345L,
                  chains = 2, cores = 2)
fit4 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                  df,
                  verbose = -1L, warmup = 7, iter = 13,
                  bart_args = list(n.trees = 11),
                  seed = 12345L,
                  chains = 2, cores = 2)


expect_equal(fit3$bart_train, fit4$bart_train)

expect_true(any(fit1$bart_train != fit3$bart_train))

# skip's parametric element takes that many transitions per sweep during
# sampling, so it moves the draws; warmup keeps one transition per sweep at any
# skip, so the frozen tuning is untouched.
fitSkip1 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                      df,
                      verbose = -1L, warmup = 7, iter = 13,
                      bart_args = list(n.trees = 11),
                      seed = 12345L, skip = c(bart = 1L, stan = 1L),
                      chains = 2, cores = 1)
fitSkip3 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2),
                      df,
                      verbose = -1L, warmup = 7, iter = 13,
                      bart_args = list(n.trees = 11),
                      seed = 12345L, skip = c(bart = 1L, stan = 3L),
                      chains = 2, cores = 1)

expect_equal(fitSkip1$bart_train, fit1$bart_train)
expect_true(any(fitSkip3$bart_train != fitSkip1$bart_train))
expect_equal(fitSkip3$adaptation$step_size, fitSkip1$adaptation$step_size)
