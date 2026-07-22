# as.array / as.matrix shape and value parity

source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

testData <- generateFriedmanData(60, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))
rm(testData)

n.chains <- 2L
n.warmup <- 3L
n.iter   <- 8L

fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 cores = 1, verbose = -1L, chains = n.chains, warmup = n.warmup, iter = n.iter,
                 bart_args = list(n.trees = 5))

arr <- as.array(fit)
mat <- as.matrix(fit)

# as.array is (iterations, chain, parameters); as.matrix collapses the first
# two into a single (iterations x chains) dimension, one column per parameter
expect_equal(dim(arr), c(n.iter - n.warmup, n.chains, dim(arr)[3L]))
expect_equal(dim(mat), c((n.iter - n.warmup) * n.chains, dim(arr)[3L]))

expect_equal(names(dimnames(arr)), c("iterations", "chain", "parameters"))
expect_equal(names(dimnames(mat)), c("iterations:chains", "parameters"))

expect_equal(colnames(mat), dimnames(arr)[[3L]])
expect_true("sigma" %in% colnames(mat))

# as.matrix is just as.array with the leading two dimensions flattened,
# preserving the per-parameter values exactly
for (par in dimnames(arr)[[3L]])
  expect_equal(as.vector(mat[,par]), as.vector(arr[,,par]))

# include_warmup is passed through to as.array()
arr_warmup_only <- as.array(fit, include_warmup = "only")
arr_both <- as.array(fit, include_warmup = TRUE)
mat_both <- as.matrix(fit, include_warmup = TRUE)

expect_equal(dim(arr_warmup_only), c(n.warmup, n.chains, dim(arr)[3L]))
expect_equal(dim(arr_both), c(n.iter, n.chains, dim(arr)[3L]))
expect_equal(dim(mat_both), c(n.iter * n.chains, dim(arr)[3L]))

for (par in dimnames(arr)[[3L]]) {
  expect_equal(unname(rbind(arr_warmup_only[,,par], arr[,,par])), unname(arr_both[,,par]))
  expect_equal(as.vector(mat_both[,par]), as.vector(arr_both[,,par]))
}

# invalid include_warmup values are rejected the same way by both generics
expect_error(as.array(fit, include_warmup = "bogus"))
expect_error(as.array(fit, include_warmup = c(TRUE, FALSE)))
expect_error(as.array(fit, include_warmup = 1L))
expect_error(as.matrix(fit, include_warmup = "bogus"))

# single chain still yields a well-formed 2-D as.array and as.matrix
fit1 <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                  cores = 1, verbose = -1L, chains = 1, warmup = n.warmup, iter = n.iter,
                  bart_args = list(n.trees = 5))
arr1 <- as.array(fit1)
mat1 <- as.matrix(fit1)
expect_equal(dim(arr1), c(n.iter - n.warmup, 1L, dim(arr1)[3L]))
expect_equal(dim(mat1), c(n.iter - n.warmup, dim(arr1)[3L]))
expect_equal(as.vector(mat1[,"sigma"]), as.vector(arr1[,,"sigma"]))
