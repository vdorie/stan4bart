This package is a early implementation of a C++ sampler that uses BART for non-parametric mean components and Stan for multilevel/parametric ones.

To install:
  1. Install the developer tools for your platform ([Mac OS](https://cran.r-project.org/bin/macosx/tools), [Windows](https://cran.r-project.org/bin/windows/Rtools/)). Mac OS users will need the (linked) gfortan for their respective platforms.
  2. Execute:

```
install.packages("remotes")
remotes::install_github("vdorie/dbarts")
remotes::install_github("vdorie/stan4bart")
```

Windows compilation has been minimally tested and may need more attention.

Here's some test code to get started with the package. It accepts an lme4 syntax and should be pretty flexible in that regard. See the package documentation `?mstan4bart` and `?stan4bart::mstan4bart-generics` for more information.

```R
require(stan4bart)

generateFriedmanData <- function(n, ranef = FALSE, causal = FALSE, binary = FALSE) {
  f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 5 * x[,5]
  
  set.seed(99)
  sigma <- 1.0
  
  x <- matrix(runif(n * 10), n, 10)
  mu <- f(x)
  
  result <- list(x = x, sigma = sigma)
  
  if (ranef) {
    n.g.1 <- 5L
    n.g.2 <- 8L
    
    result <- within(result, {
      g.1 <- sample(n.g.1, n, replace = TRUE)
      
      Sigma.b.1 <- matrix(c(1.5^2, .2, .2, 1^2), 2)
      R.b <- chol(Sigma.b.1)
      b.1 <- matrix(rnorm(2 * n.g.1), n.g.1) %*% R.b
      rm(R.b)
      
      g.2 <- sample(n.g.2, n, replace = TRUE)
      
      Sigma.b.2 <- as.matrix(1.2)
      b.2 <- rnorm(n.g.2, 0, sqrt(Sigma.b.2))
      
      mu.fixef <- x[,4] * 10
      mu.bart <- mu - mu.fixef
      mu.ranef <- b.1[g.1,1] + x[,4] * b.1[g.1,2] + b.2[g.2]
      mu <- mu + mu.ranef
    })
  } else {
    mu.fixef <- x[,4] * 10
    mu.bart <- mu - mu.fixef 
  }
  
  if (causal) {
    result <- within(result, {
      tau <- 5
      z <- rbinom(n, 1, 0.2)
    })
    
    if (ranef) {
      result <- within(result, {
        mu.fixef.0 <- mu.fixef
        mu.fixef.1 <- mu.fixef.0 + tau
        mu.bart.0  <- mu.bart.1  <- mu.bart
        mu.ranef.0 <- mu.ranef.1 <- mu.ranef
        
        mu.0 <- mu.bart.0 + mu.fixef.0 + mu.ranef.0
        mu.1 <- mu.bart.1 + mu.fixef.1 + mu.ranef.1
      })
    } else {
      result <- within(result, {
        mu.fixef.0 <- mu.fixef
        mu.fixef.1 <- mu.fixef.0 + tau
        mu.bart.0  <- mu.bart.1  <- mu.bart
        
        mu.0 <- mu.bart.0 + mu.fixef.0
        mu.1 <- mu.bart.1 + mu.fixef.1
      })
    }
    
    if (binary) {
      result <- within(result, {
        loc   <- mean(c(mu.0, mu.1))
        scale <- sd(c(mu.0, mu.1)) / qnorm(0.15)
        mu.0 <- (mu.0 - loc) / scale
        mu.1 <- (mu.1 - loc) / scale
        
        mu.fixef.0 <- (mu.fixef.0 - loc) / scale
        mu.fixef.1 <- (mu.fixef.1 - loc) / scale
        mu.bart.0 <- mu.bart.0 / scale
        mu.bart.1 <- mu.bart.1 / scale
        if (ranef) {
          mu.ranef.0 <- mu.ranef.0 / scale
          mu.ranef.1 <- mu.ranef.1 / scale
        }
        
        rm(loc, scale)
        
        y.0 <- rbinom(n, 1L, pnorm(mu.0))
        y.1 <- rbinom(n, 1L, pnorm(mu.1))
        y <- y.1 * z + y.0 * (1 - z)
      })
    } else {
      result <- within(result, {
        y.0 <- mu.0 + rnorm(n, 0, sigma)
        y.1 <- mu.1 + rnorm(n, 0, sigma)
        y <- y.1 * z + y.0 * (1 - z)
      })
    }
    
    
    result$mu <- NULL
    result$mu.fixef <- NULL
    result$mu.ranef <- NULL
  } else {
    if (binary) {
      result <- within(result, {
        loc <- mean(mu)
        scale <- sd(mu) / qnorm(0.15)
        mu <- (mu - loc) / scale
        
        mu.fixef <- (mu.fixef - loc) / scale
        mu.bart <- mu.bart / scale
        if (ranef)
          mu.ranef <- mu.ranef / scale
        
        rm(loc, scale)
        
        y <- rbinom(n, 1L, pnorm(mu))
      })
    } else {
      result <- within(result, {
        y <- mu + rnorm(n, 0, sigma)
      })
    }
  }
  
  result
}
testData <- generateFriedmanData(1000, TRUE, TRUE, FALSE)
rm(generateFriedmanData)

df <- with(testData, data.frame(x, g.1, g.2, y, z))


# Causal inference example
fit <- mstan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                  cores = 1, verbose = 1,
                  treatment = z)

samples.mu.train <- extract(fit)
samples.mu.test  <- extract(fit, sample = "test")

# Individual conditional treatment effects
samples.icate <- (samples.mu.train - samples.mu.test) * (2 * testData$z - 1)
# Conditional average treatment effect
samples.cate <- apply(samples.icate, 2, mean)
cate <- mean(samples.cate)

samples.ppd.test <- extract(fit, type = "ppd", sample = "test")

# Individual sample treatment effects
samples.ite <- (testData$y - samples.ppd.test) * (2 * testData$z - 1)
# Sample average treatment effect
samples.sate <- apply(samples.ite, 2, mean)
sate <- mean(samples.sate)

# Population average treatment effect
samples.ppd.test <- extract(fit, type = "ppd", sample = "train")
samples.pate <- apply((samples.ppd.test - samples.ppd.test) * (2 * testData$z - 1), 2, mean)
pate <- mean(samples.pate)


fitted.mu.train <- fitted(fit)
# equal to: apply(samples.mu.train, 1, mean)
fitted.mu.test  <- fitted(fit, sample = "test")
# equal to: apply(samples.mu.test,  1, mean)

# Observed and conterfactual MSE
mse.train <- with(testData, mean((fitted.mu.train - mu.1 * z - mu.0 * (1 - z))^2))
mse.test  <- with(testData, mean((fitted.mu.test  - mu.1 * (1 - z) - mu.0 * z)^2))
```

