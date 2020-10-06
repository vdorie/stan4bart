This mess is an initial stab at a C++ version of sampler that optionally uses BART for non-parametric mean components and Stan for multilevel/parametric ones.

To install:
  1. Install all prerequisites (see DESCRIPTION file)
  2. Try install from this root directory

Windows not tested yet.

Here's some test code to get started with the package. It accepts an lme4 syntax and should be pretty flexible in that regard, but the results are not well packaged at present. You will have to manually dig into the samples, but they should at least have a familar structure.

```R
require(stan4bart)

generateFriedmanData <- function(n = 100) {
  f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
  
  set.seed(99)
  sigma <- 1.0
  
  x <- matrix(runif(n * 10), n, 10)
  y <- rnorm(n, f(x), sigma)
  
  list(x = x, y = y, f.x = f(x))
}
testData <- generateFriedmanData(1000)
rm(generateFriedmanData)


set.seed(0)

n.g.1 <- 5L
g.1 <- sample(n.g.1, length(testData$y), replace = TRUE)

Sigma.b <- matrix(c(1.5^2, .2, .2, 1^2), 2)
R.b <- chol(Sigma.b)
b.1 <- matrix(rnorm(2 * n.g.1), n.g.1) %*% R.b

n.g.2 <- 8L
g.2 <- sample(n.g.2, length(testData$y), replace = TRUE)

sigma.b <- 1.2
b.2 <- rnorm(n.g.2, 0, sigma.b)

testData$y <- testData$y + b.1[g.1,1] + testData$x[,4] * b.1[g.1,2] + b.2[g.2]
testData$g.1 <- g.1
testData$b.1 <- b.1
testData$g.2 <- g.2
testData$b.2 <- b.2
rm(Sigma.b, R.b, g.1, b.1, g.2, b.2, sigma.b)

df <- with(testData, data.frame(x, g.1, g.2, y))

ranef_true <- with(testData, b.1[g.1,1] + x[,4] * b.1[g.1,2] + b.2[g.2])
fixef_true <- testData$f.x
mean_true <- fixef_true + ranef_true

fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4) + X4 + (1 + X4 | g.1) + (1 | g.2), df,
                  verbose = 2, chains = 1)

# samples from the posterior of the expected value of the response
mean((fitted(fit) - mean_true)^2)

mean((fitted(fit, "mu.fixef") + fitted(fit, "bart") - fixef_true)^2)
mean((fitted(fit, "mu.ranef") - ranef_true)^2)

lmer_fit <- lme4::lmer(y ~ . - g.1 - g.2 - (1 + X4 | g.1) + (1 | g.2), df)
mean((fitted(lmer_fit) - mean_true)^2)
```
