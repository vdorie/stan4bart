loadIHDPData <- function(root.path = ".", file.name = "sim.data.gz")
{
  ihdp <- read.table(gzfile(file.path(root.path, file.name)), header = TRUE, sep = "\t")
  ihdp <- subset(ihdp, treat != 1 | momwhite != 0)
  
  
  covs.cont <- c("bw", "b.head", "preterm", "birth.o", "nnhealth", "momage")
  covs.cat  <- c("sex", "twin", "b.marr", "mom.lths", "mom.hs", "mom.scoll",
                 "cig", "first", "booze", "drugs", "work.dur", "prenatal",
                 "ark", "ein", "har", "mia", "pen", "tex", "was")
  
  z <- ihdp$treat
  x <- as.matrix(ihdp[,c(covs.cont, covs.cat)])
  ps <- fitted(glm(z ~ ., data.frame(z, x), family = binomial))
  
  g1 <- ihdp$momage
  g1[g1 < 16] <- 15
  g1[g1 > 39] <- 40
  g1 <- as.factor(g1)
  g2 <- as.factor(ihdp$site.num)

  
  p <- length(c(covs.cont, covs.cat))
  
  trans.x <- list(m = rep_len(0, p), s = rep_len(1, p))
  trans.x <- within(trans.x, {
    m[seq_along(covs.cont)] <- apply(x[,covs.cont], 2L, mean)
    s[seq_along(covs.cont)] <- apply(x[,covs.cont], 2L, sd)
  })
  
  x.z <- with(trans.x, t(t(x) - m) %*% diag(1 / s))
  colnames(x.z) <- colnames(x)
  
  ### create matrix of all interactions etc for third response surface
  formula <- as.formula(paste0("y.tmp ~ -1 + (", paste0(c(covs.cont, covs.cat), collapse = " + "), ")^2 + ", paste0(paste0("I(", covs.cont, ")^2"), collapse = " + ")))
  y.tmp <- rnorm(nrow(x.z))
  temp <- glm(formula, as.data.frame(x.z), x = TRUE, family = gaussian)
  x.o2 <- temp$x[,!is.na(temp$coef)]
  x.o2 <- x.o2[,colnames(x.o2) %not_in% colnames(x.z)]
  
  # returns:
  #  x       - ihdp data, original form
  #  x.z     - ihdp data, standardized
  #  x.o2    - second order terms (interactions and quadratics) of standardized matrix
  #  z       - treament
  #  g1      - first grouping variable
  #  trans.x - mean and scale used to transform x into x.z
  nlist(x, x.z, x.o2, z, ps, g1, g2, trans.x)
}

