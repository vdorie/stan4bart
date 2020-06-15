# change to the path where you've checked out the repo and subdir to ihdp
setwd("~/Repositories/stan4bart/ihdp")

require(bartCause) # requires the development version from github
require(stan4bart)

source("util.R")
source("data.R")
source("sim.R")

ihdp <- loadIHDPData()

n.iters <- 1000L

for (iter in seq_len(n.iters)) {
  resp <- generateResponseForIter(ihdp, iter, momage.is.group = TRUE)
  # use response surface C, noted below by suffixes
  # mu.xx are the mean structure for trt/control
  # b.xx are the varying intercepts and slopes for each observation
  y.0 <- with(resp, mu.0.c + b.0 + eps.0.c)
  y.1 <- with(resp, mu.1.c + b.1 + eps.1.c)
  y   <- with(ihdp, y.0 * (1 - z) + y.1 * z)
  
  df <- as.data.frame(ihdp$x.z)
  df$g1 <- ihdp$g1
  df$z  <- ihdp$z
  df$y  <- y
  
  # see R/mstan4bart.R for a description of what this function returns
  fit.s4b <- mstan4bart(y ~ . - g1 - (1 + z | g1), df, treatment = z, verbose = 0, chains = 10)
  rows.b <- grepl("^b\\.", rownames(fit.s4b[[1L]]$sample$stan$raw))
  n.obs <- nrow(df)
  n.samp <- ncol(fit.s4b[[1L]]$sample$bart$test)
  Zt.cf <- attr(fit.s4b, "Zt.cf")
 
  
  # by way of example, this computes the individual treatment effect estimates,
  #   (y(obs) - y(cf)) * sign flip for controls
  # can be modified to compute sub-group estimates or confidence intervals
  ite <- apply(sapply(fit.s4b, function(fit.chain) {
    sample.sigma <- fit.chain$sample$stan$raw["aux",]
    y.hat.cf <- fit.chain$sample$bart$test + Matrix::crossprod(Zt.cf, fit.chain$sample$stan$raw[rows.b,,drop = FALSE]) +
      rnorm(n.obs * n.samp, 0, rep(sample.sigma, rep(n.obs, n.samp)))
    
    sample.ite <- (y - y.hat.cf) * (2 * ihdp$z - 1) # num obs x num samples
    ite <- apply(sample.ite, 1, mean)
  }), 1L, mean)
  # mean(ite) is the sate
  
  
  # bart with fixed effects
  confounders <- paste0(colnames(ihdp$x.z), collapse = " + ")
  fit.bart.fixef <- bartc(y, z, confounders, df, estimand = "ate", group.by = g1, use.ranef = FALSE, n.samples = 2000, n.burn = 1000, verbose = FALSE)
  sum.bart.fixef <- summary(fit.bart.fixef, target = "sate")
  # sum.bart.fixef$estimates$estimate
  
  # bart with varying intercepts
  fit.bart.varint <- bartc(y, z, confounders, df, estimand = "ate", group.by = g1, use.ranef = TRUE, n.samples = 2000, n.burn = 1000, verbose = FALSE)
  sum.bart.varint <- summary(fit.bart.varint, target = "sate")
  # sum.bart.varint$estimates$estimate
  
  
  break
}
