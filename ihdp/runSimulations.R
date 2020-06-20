# change to the path where you've checked out the repo and subdir to ihdp
setwd("~/Repositories/stan4bart/ihdp")

require(bartCause) # requires the development version from github
require(stan4bart)

source("util.R")
source("data.R")
source("sim.R")

ihdp <- loadIHDPData()

# methods <- c("lm", "lmer", "bart", "bart_vi", "stan4bart")
method.files <- list.files("methods")
n.methods <- length(method.files)
methods <- gsub("\\.R$", "", method.files)

# load methods
fit.functions <- vector("list", length(n.methods))
fit.caches    <- vector("list", length(n.methods))
for (i in seq_along(method.files)) {
  source.env <- new.env()
  source(file.path("methods", method.files[[i]]), source.env)
  source.env$init()
  fit.functions[[i]] <- source.env$getFit
  fit.caches[[i]]    <- source.env$getCache(ihdp)
}
names(fit.functions) <- names(fit.caches) <- methods

n.iters <- 1000L
n.groups  <- nlevels(ihdp$g1)
p.groups  <- unname(table(ihdp$g1[ihdp$z == 1])) / sum(ihdp$z == 1)

# vector of treated indices corresponding to each group
g.sel <- lapply(seq_len(n.groups), function(j) which(ihdp$g1 == levels(ihdp$g1)[j] & ihdp$z == 1))

results <- array(NA_real_,
                 c(n.iters, n.methods, 5L),
                  dimnames = list(NULL, methods, c("bias", "cover", "ci_len", "pehe", "pegste")))

truthIsFour <- TRUE

set.seed(1)
# sample a covariance matrix for random effects by first picking
# a high correlation
# then a ratio of variances, while fixing the det of Sigma.b to 1
rho <- rbeta(1, 16, 4)
r.var <- rf(1, 5, 7)
# r.var = sigma1.sq / sigma2.sq
#  | Sigma.b | = sigma1.sq * sigma2.sq - sigma1.sq * sigma2.sq * rho^2
#              = sigma1.sq * sigma2.sq * (1 - rho^2)
# 1 / (1 - rho) = sigma1.sq^2 / r.var
sigma1.sq <- sqrt(r.var / (1 - rho^2))
sigma2.sq <- sigma1.sq / r.var
Sigma.b <- matrix(c(sigma1.sq, rho * sqrt(sigma1.sq) * sqrt(sigma2.sq), rho * sqrt(sigma1.sq) * sqrt(sigma2.sq), sigma2.sq),
                  2, 2)
ranef.dist <- function(n) rnorm(n, 0, 1)
#######################
# Simulation settings #
#######################
sim.setting <- 1
# 1 - correlated covariance, rho ~= 0.86, |Sigma.b| = 1, group coefficients normally distributed
# 2 - |Sigma.b| = 2
# 3 - |Sigma.b| = 1, group coefficients t-7 distributed
if (sim.setting == 2)
  Sigma.b <- Sigma.b * 2
if (sim.setting == 3)
  ranef.dist <- function(n) rt(n, 7, 0)

rm(rho, r.var, sigma1.sq, sigma2.sq)


startIter <- which.max(apply(results[,,"bias"], 1L, anyNA))
startTime <- proc.time()

for (iter in seq.int(startIter, n.iters)) {
  if (!anyNA(results[iter,,"bias"])) next
  
  if (iter %% 10L == 0L) {
    timeDiff <- proc.time() - startTime
    
    cat("fitting iter: ", iter,
        ", elapsed time: ", format.time(timeDiff[["elapsed"]]), 
        ", time/iter: ", format.time(timeDiff[["elapsed"]] / (iter - startIter + 1L)),
        "\n", sep = "")
  }
  resp <- generateResponseForIter(ihdp, iter, grouping.var = "momage", Sigma.b = Sigma.b, ranef.dist = ranef.dist)
  
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
  
  icate.truth <- with(resp, (mu.1.c + b.1 - mu.0.c - b.0))
  gcatt.truth <- sapply(g.sel, function(sel) mean(icate.truth[sel]))
  icatt.truth <- icate.truth[df$z == 1]
  truth <- if (truthIsFour) 4 else mean(icatt.truth)
  
  sd.y <- sd(y)
  
  for (method in methods) {
    fit <- fit.functions[[method]](df, fit.caches[[method]])
    
    results[iter,method,"bias"]  <- (fit$catt - truth) / sd.y
    results[iter,method,"cover"] <- as.double(fit$catt.lower <= truth && truth <= fit$catt.upper)
    results[iter,method,"ci_len"] <- (fit$catt.upper - fit$catt.lower) / sd.y
    results[iter,method,"pehe"]  <- sqrt(mean((fit$icatt - icatt.truth)^2)) / sd.y
    results[iter,method,"pegste"] <- sqrt(mean((fit$gcatt - gcatt.truth)^2)) / sd.y
  }
}

maxIter <- if (!anyNA(results[,,"bias"])) n.iters else which.max(apply(results[,,"bias"], 1L, anyNA)) - 1L

source("plotDualAxis.R")

pdf(paste0("results_sim_", sim.setting, ".pdf"), 6, 6)
par(mfrow = c(2, 2))
plotBiasRMSE(results[seq_len(maxIter),,"bias"], methods, main = "Bias/RMSE")
#plotCovLenDualAxis(results[seq_len(maxIter),,"cover"],
                   #results[seq_len(maxIter),,"ci_len"],
                   #methods, main = "Coverage/Interval Length")
coverage <- apply(results[seq_len(maxIter),,"cover"], 2L, mean)
ci_len   <- apply(results[seq_len(maxIter),,"ci_len"], 2L, mean)
par(mar = c(2.2, 2.2, 1.5, 2.5),
      mgp = c(1, 0.2, 0.0),
      tcl = -0.3,
      cex.main = 1.2, cex.axis = 0.7)
plot(NULL, type = "n", xlim = c(0.9, 1.02), ylim = range(ci_len, na.rm = TRUE),
     xlab = "coverage", ylab = "interval length", main = "Coverage/Interval Length")
abline(v = 0.95, col = "gray")
text(coverage, ci_len, methods)
par(mar = c(1.5, 2.5, 1.5, 2.5),
      mgp = c(1, 0.2, 0.0),
      tcl = -0.3,
      cex.main = 1.2, cex.axis = 0.7)
boxplot(results[seq_len(maxIter),,"pehe"], main = "PEHE", ylab = "pehe")
boxplot(results[seq_len(maxIter),,"pegste"], main = "PEGSTE", ylab = "pegste")
dev.off()

# old plots
if (FALSE) {
pdf(paste0("results_sim_", sim.setting, ".pdf"), 6, 6)
par(mfrow = c(2, 2), cex.axis = 0.8)
boxplot(results[seq_len(maxIter),,"bias"], main = "bias")

boxplot(results[seq_len(maxIter),,"pehe"], main = "pehe")
boxplot(results[seq_len(maxIter),,"pegste"], main = "pegste")
coverage <- apply(results[seq_len(maxIter),,"cover"], 2L, mean)
rmse     <- apply(results[seq_len(maxIter),,"bias"], 2L, function(x) sqrt(mean(x^2)))

plot(NULL, type = "n", xlim = range(coverage, na.rm = TRUE), ylim = range(rmse),
     xlab = "coverage", ylab = "rmse")
text(coverage, rmse, methods)
dev.off()
}

save(results, file = paste0("results_sim_", sim.setting, ".RData"))


