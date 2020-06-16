# uses post-stratification to compute catt
getFit <- function(df, cache) {
  p.groups <- cache$p.groups  

  fit <- lm(y ~ . + g1 * z, df)
  
  icatt <- fitted(fit)[df$z == 1] - predict(fit, within(df, z <- 1 - z)[df$z == 1,])
  
  # not being clever about this
  new.data <- df[sapply(seq_len(n.groups), function(j) which.max(df$g1 == levels(df$g1)[j])),]
  new.data$z <- 1
  v <- model.matrix(fit$terms, new.data)
  new.data$z <- 1 - new.data$z
  v <- v - model.matrix(fit$terms, new.data)
  # est is prop' v beta
  gcatt <- as.vector(v %*% fit$coef)
  catt <- sum(gcatt * p.groups)
  
  # se is prop'v beta beta' v' prop
  se.lm <- sqrt(crossprod(p.groups, v %*% tcrossprod(vcov(fit), v)) %*% p.groups)[1L]
  
  catt <- sum(gcatt * p.groups)

  
  list(catt = catt,
       catt.lower = catt - qnorm(0.975) * se.lm,
       catt.upper = catt + qnorm(0.975) * se.lm,
       icatt = icatt,
       gcatt = gcatt)
}

getCache <- function(ihdp)
{
  nlist(p.groups = unname(table(ihdp$g1[ihdp$z == 1])) / sum(ihdp$z == 1))
}

init <- function() {
  invisible(NULL)
}
