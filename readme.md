This package is an implementation of a C++ sampler that uses [BART](https://cran.r-project.org/package=dbarts) for non-parametric mean components and [Stan](https://mc-stan.org) for multilevel/parametric ones.

## Installation

  1. Install the developer tools for your platform ([Mac OS](https://cran.r-project.org/bin/macosx/tools), [Windows](https://cran.r-project.org/bin/windows/Rtools/)). Mac OS users will need the (linked) gfortan for their respective platforms.
  2. Execute:

```
if (length(find.package("remotes", quiet = TRUE)) == 0L)
  install.packages("remotes")
remotes::install_github("vdorie/dbarts")
remotes::install_github("vdorie/stan4bart")
```

## Use

The package utilizes the flexible, expressive lme4 syntax for specifying group-level structures. See the package documentation `?stan4bart` and `?stan4bart::stan4bart-generics` for more information.

### Formulas

* `y ~ bart(x_nuisance_1 + x_nuisance_1) + x_inference` - wrap the formula for any variables that you want to be fit non-parametrically in a call to `bart()`; `+`s within `bart()` are interpretted figuratively as indicating the names of the variables to be included
* `y ~ bart(x_fixef) + (1 + x_ranef | g)` - supports all [`lme4`](https://www.rdocumentation.org/packages/lme4/versions/1.1-26/topics/lmer) varying intercept and slope formula constructs
* `y ~ x_1 + (1 + x_2 | g)` - if there's no BART component, use [`rstanarm`](https://cran.r-project.org/package=rstanarm) or [`lme4`](https://cran.r-project.org/package=lme4) instead

### Main Event

The main function is `stan4bart`.

### Results

Results are retrieved using the `extract`, `fitted`, and `predict` generics. See `?"stan4bart-generics"` for more information.

## Known Issues

The name and definition of `extract` conflict with `rstan`. The `rstan` package is not needed to use `stan4bart` and does not need to be loaded. If a name-collision occurs, the `stan4bart` extract can be referenced as in:

```R
stan4bart:::extract.stan4bartFit(stan4bart_fit)
```

## Example Code

```R
library(stan4bart)

# Load a test-data function
source(system.file("common", "friedmanData.R", package = "stan4bart"), local = TRUE)

# Relatively low n for illustrative purposes
testData <- generateFriedmanData(n = 100, ranef = TRUE, causal = TRUE, binary = FALSE)

# First level model is:
#   y ~ f(x_1, x_2) + a * x_3^2 + b * x_4 + c * x_5 + z
# Random intercepts are added for g.1 and g.2, and a random slope is placed on x4
# x_6 through x_10 are pure noise
df <- with(testData, data.frame(x, g.1, g.2, y, z))

# Causal inference example
fit <- stan4bart(y ~ bart(. - g.1 - g.2 - X4 - z) + X4 + z + (1 + X4 | g.1) + (1 | g.2), df,
                 treatment = z,
                 cores = 1, seed = 0,
                 verbose = 1)

samples.mu.train <- extract(fit)
samples.mu.test  <- extract(fit, sample = "test")

# Individual conditional treatment effects
samples.icate <- (samples.mu.train - samples.mu.test) * (2 * testData$z - 1)
# Conditional average treatment effect
samples.cate <- apply(samples.icate, 2, mean)
cate <- mean(samples.cate)
cate.int <- c(cate - 1.96 * sd(samples.cate), cate + 1.96 * sd(samples.cate))

# Samples of the posterior predictive distribution are used in calculating
# the counterfactuals for SATE and for calculating the response under
# the observed treatment condition when estimating PATE.
samples.ppd.test <- extract(fit, type = "ppd", sample = "test")

# Individual sample treatment effects
samples.ite <- (testData$y - samples.ppd.test) * (2 * testData$z - 1)
# Sample average treatment effect
samples.sate <- apply(samples.ite, 2, mean)
sate <- mean(samples.sate)
sate.int <- c(sate - 1.96 * sd(samples.sate), sate + 1.96 * sd(samples.sate))

# Population average treatment effect
samples.ppd.train <- extract(fit, type = "ppd", sample = "train")
samples.pate <- apply((samples.ppd.train - samples.ppd.test) * (2 * testData$z - 1), 2, mean)
pate <- mean(samples.pate)
pate.int <- c(pate - 1.96 * sd(samples.pate), pate + 1.96 * sd(samples.pate))


fitted.mu.train <- fitted(fit)
# equal to: apply(samples.mu.train, 1, mean)
fitted.mu.test  <- fitted(fit, sample = "test")
# equal to: apply(samples.mu.test,  1, mean)

# Observed and conterfactual MSE
mse.train <- with(testData, mean((fitted.mu.train - mu.1 * z - mu.0 * (1 - z))^2))
mse.test  <- with(testData, mean((fitted.mu.test  - mu.1 * (1 - z) - mu.0 * z)^2))
```
