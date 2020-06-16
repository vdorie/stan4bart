
# uses post-stratification to compute catt
# TODO: compute standard error
getFit <- function(df, cache) {
  
  p.groups <- cache$p.groups  

  form <- as.formula(paste0("y ~ ", paste0(colnames(df)[colnames(df) %not_in% c("y", "g1")], collapse = " + "), " + (1 + z | g1)"))
  
  ignored <- capture.output(fit <- suppressWarnings(lmer(form, df)))
  
  icatt <- fitted(fit)[df$z == 1] - predict(fit, within(df, z <- 1 - z)[df$z == 1,])
  
  # not being clever about this
  new.data <- df[sapply(seq_len(n.groups), function(j) which.max(df$g1 == levels(df$g1)[j])),]
  new.data$z <- 1
  mu.1 <- predict(fit, new.data)
  new.data$z <- 0
  mu.0  <- predict(fit, new.data)
  
  gcatt <- mu.1 - mu.0
  catt <- sum(gcatt * p.groups)
  
  se.lmer <- NA_real_
  
  catt <- sum(gcatt * p.groups)

  
  list(catt = catt,
       catt.lower = catt - qnorm(0.975) * se.lmer,
       catt.upper = catt + qnorm(0.975) * se.lmer,
       icatt = icatt,
       gcatt = gcatt)
}

getCache <- function(ihdp)
{
  nlist(p.groups = unname(table(ihdp$g1[ihdp$z == 1])) / sum(ihdp$z == 1))
}

init <- function() {
  if (!require(lme4, quietly = TRUE)) {
    install.packages("lme4")
    if (!require(lme4, quietly = TRUE))
      stop("lme4 must be installed")
  }
  invisible(NULL)
}
