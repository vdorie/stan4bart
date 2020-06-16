# bart with varying intercepts but no slopes
getFit <- function(df, cache) {
  confounders <- paste0(colnames(df)[colnames(df) %not_in% c("z", "y")], collapse = " + ")
  fit <- bartc(y, z, confounders, df, estimand = "att", group.by = g1, use.ranef = TRUE,
               n.samples = 2000, n.burn = 1000, verbose = FALSE)
  icate.samples <- extract(fit, "icate", sample = "all")
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
  if (!require(bartCause, quietly = TRUE) || packageVersion("bartCause") < "1.0.3") {
    remotes::install_github("vdorie/bartCause")
    if (!require(bartCause, quietly = TRUE))
      stop("bartCause must be installed")
  }
  if (packageVersion("bartCause") < "1.0.3") 
    stop("bartCause needs to be updated, execute remotes::install_github('vdorie/bartCause')")
  
  invisible(NULL)
}
