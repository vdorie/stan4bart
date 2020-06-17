
# could use a better way of dividing the different simulation settings
#   if momage is a grouping variable, it shouldn't be used in the mean
#   function; however, doing so removes consistency with historical
#   sims
#
# if we don't care about matching old sims, best to break into 
#   a function that only generates a single response surface at a time
#   e.g., generateResponseForIter(ihdp, iter, mean.type = "a", group.type = "b")
#
# if we do care, then have to generate all three mean-response surfaces
#   sequentially after setting the seed


# NOTE: !!if group.var is not momage, references to g1 will need to be
# changed to g2 elsewhere!!
generateResponseForIter <- function(ihdp, iter,
                                    grouping.var = c("momage", "site.num"),
                                    Sigma.b = diag(c(1, 0.5)))
{
  # only warning is that the sample.kind is old, which is intended to
  # match the original dgp
  if (iter <= 500)
    suppressWarnings(set.seed(565 + iter * 5, sample.kind = "Rounding"))
  else
    suppressWarnings(set.seed(7565 + iter * 5, sample.kind = "Rounding"))
  
  grouping.var <- match.arg(grouping.var)
  
  if (grouping.var == "momage") {
    ihdp$x.z  <- ihdp$x.z[,colnames(ihdp$x.z) != "momage"]
    ihdp$x.o2 <- ihdp$x.o2[,!grepl("momage", colnames(ihdp$x.o2))]
  }
  
  with(ihdp, {
    n <- nrow(x.z)
    p <- ncol(x.z)
    
    sig.y <- 1
    tau <- 4
    
    beta.a <- sample(0:4, p + 1L, replace = TRUE, prob = c(.5, .2, .15, .1, .05))
    
    mu.a <- cbind(1, x.z) %*% beta.a
    eps.0.a <- rnorm(n, 0, sig.y)
    eps.1.a <- rnorm(n, 0, sig.y)
     
    beta.b <- sample(c(0.0, 0.1, 0.2, 0.3, 0.4), p + 1L, replace = TRUE, prob = c(.6, .1, .1, .1, .1))
    
    eps.0.b <- rnorm(n, 0, sig.y)
    eps.1.b <- rnorm(n, 0, sig.y)
    
    # main effects coefficients
    beta.m.0.c <- sample(c(0, 1, 2), p + 1L, replace = TRUE, prob = c(.6, .3, .1))
    beta.m.1.c <- sample(c(0, 1, 2), p + 1L, replace = TRUE, prob = c(.6, .3, .1))
    # quadratic coefficients
    #these we make pretty rare since they really represent 3-way interactions
    beta.q.0.c <- sample(c(0, .5, 1), ncol(x.o2), replace = TRUE, prob = c(.8, .15, .05))
    beta.q.1.c <- sample(c(0, .5, 1), ncol(x.o2), replace = TRUE, prob = c(.8, .15, .05))
    eps.0.c <- rnorm(n, 0, sig.y)
    eps.1.c <- rnorm(n, 0, sig.y)
    
    ##################
    #  group stuff
    ##################
    g <- if (grouping.var == "momage") g1 else g2
    
    n.g <- length(unique(g))
    b <- matrix(rnorm(n.g * 2, 0, 1), n.g) %*% chol(Sigma.b)
    
    b.0 <- b[as.integer(g1),1L]
    b.1 <- b.0 + b[as.integer(g2),2L]
        
    mu.0.b <- exp(cbind(1, x.z + 0.5) %*% beta.b)
    mu.1.b <-     cbind(1, x.z + 0.5) %*% beta.b
    # offset.b <- mean(mu.1.b[z == 1] - mu.0.b[z == 1]) - 4
    offset.b <- mean(mu.1.b[z == 1] + b.1[z == 1] - mu.0.b[z == 1] - b.0[z == 1]) - 4
    mu.1.b <- mu.1.b - offset.b
    
    mu.0.c <- cbind(1, x.z) %*% beta.m.0.c + x.o2 %*% beta.q.0.c
    mu.1.c <- cbind(1, x.z) %*% beta.m.1.c + x.o2 %*% beta.q.1.c
    #offset.c <- mean(mu.1.c[z == 1] - mu.0.c[z == 1]) - 4
    offset.c <- mean(mu.1.c[z == 1] + b.1[z == 1] - mu.0.c[z == 1] - b.0[z == 1]) - 4
    mu.1.c <- mu.1.c - offset.c
    
    nlist(mu.0.a = mu.a, mu.1.a = mu.a + tau, eps.0.a, eps.1.a,
          mu.0.b, mu.1.b, eps.0.b, eps.1.b,
          mu.0.c, mu.1.c, eps.0.c, eps.1.c,
          b, b.0, b.1)
  })
}
