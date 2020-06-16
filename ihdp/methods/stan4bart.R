# bart with unmodeled group parameters
getFit <- function(df, cache) {
  
  # see R/mstan4bart.R for a description of what this function returns
  fit.s4b <- mstan4bart(y ~ . - g1 - (1 + z | g1), df, treatment = z, verbose = 0, chains = 10)
  rows.b <- grepl("^b\\.", rownames(fit.s4b[[1L]]$sample$stan$raw))
  n.obs <- nrow(df)
  n.samp <- ncol(fit.s4b[[1L]]$sample$bart$test)
  n.chain <- length(fit.s4b)
  Zt.obs <- attr(fit.s4b, "Zt.obs")
  Zt.cf  <- attr(fit.s4b, "Zt.cf")
 
  
  # by way of example, this computes the individual treatment effect estimates,
  #   (mu(obs) - mu(cf)) * sign flip for controls
  # can be modified to compute sub-group estimates or confidence intervals
    icate.samples <- t(matrix(sapply(fit.s4b, function(fit.chain) {
    # this commented out bit computes individual sample effect estimates
    #sample.sigma <- fit.chain$sample$stan$raw["aux",]
    #y.hat.cf <- fit.chain$sample$bart$test + Matrix::crossprod(Zt.cf, fit.chain$sample$stan$raw[rows.b,,drop = FALSE]) +
    #  rnorm(n.obs * n.samp, 0, rep(sample.sigma, rep(n.obs, n.samp)))
    
    mu.hat.obs <- fit.chain$sample$bart$train + Matrix::crossprod(Zt.obs, fit.chain$sample$stan$raw[rows.b,,drop = FALSE])
    mu.hat.cf  <- fit.chain$sample$bart$test  + Matrix::crossprod(Zt.cf,  fit.chain$sample$stan$raw[rows.b,,drop = FALSE])
    
    as.vector(mu.hat.obs - mu.hat.cf * (2 * ihdp$z - 1)) # num obs x num samples
  }), n.obs, n.samp * n.chain))
  gcatt <- sapply(g.sel, function(sel) mean(apply(icate.samples[,sel,drop = FALSE], 1L, mean)))
  icatt.samples <- icate.samples[,df$z == 1]
  icatt <- apply(icatt.samples, 2L, mean)
  catt.samples <- apply(icatt.samples, 1L, mean)
  catt <- mean(catt.samples)
  catt.se <- sd(catt.samples)
    
  
  list(catt = catt,
       catt.lower = catt - qnorm(0.975) * catt.se,
       catt.upper = catt + qnorm(0.975) * catt.se,
       icatt = icatt,
       gcatt = gcatt)
}

getCache <- function(ihdp)
{
  list()
}

init <- function() {
  if (!require(stan4bart, quietly = TRUE)) {
    remotes::install_github("vdorie/stan4bart")
    if (!require(stan4bart, quietly = TRUE))
      stop("stan4bart must be installed")
  }
  
  invisible(NULL)
}
