
sim.setting <- 1L
plotIndividually <- TRUE

resultsFile <- paste0("results_sim_", sim.setting, ".rds")

if (!file.exists(resultsFile))
  stop("cannot find results file '", resultsFile, "'")

results <- readRDS(resultsFile)

n.iters <- dim(results)[1L]
methods <- dimnames(results)[[2L]]
maxIter <- if (!anyNA(results[,,"bias"])) n.iters else which.max(apply(results[,,"bias"], 1L, anyNA)) - 1L


if (!plotIndividually) {
  pdf(paste0("all_sim_", sim.setting, ".pdf"), 6, 6)
  par(mfrow = c(2, 2))
}

if (plotIndividually)
  pdf(paste0("bias_rmse_", sim.setting, ".pdf"), 3.5, 3.5)

#source("plotDualAxis.R")
#plotBiasRMSE(results[seq_len(maxIter),,"bias"], methods, main = "Bias/RMSE")

rmse <- apply(results[seq_len(maxIter),,"bias"], 2L, function(x) sqrt(mean(x^2)))

par(mar = c(1.5, 2.5, 1.5, 2.5),
      mgp = c(1, 0.2, 0.0),
      tcl = -0.3,
      cex.main = 1.2, cex.axis = 0.7)
y.range <- range(results[seq_len(maxIter),,"bias"])
y.range <- c(-1, 1) * 0.5 * 1.1 * diff(y.range) + mean(y.range)
x.range <- c(1, length(methods))
x.range <- c(-1, 1) * 0.5 * 1.2 * diff(x.range) + mean(x.range)
plot(NULL, type = "n", xlim = x.range, ylim = y.range,
     main = "Bias/RMSE", ylab = "bias", xlab = "", xaxt = "n")
abline(h = 0, col = "gray")
boxplot(results[seq_len(maxIter),,"bias"],
        ylim = y.range, staplewex = 0, add = TRUE)
plotRegion <- par("usr")
cex.text <- 0.7
lineHeight <- diff(grconvertY(0:1, "inches", "user")) * par("cin")[2L] * par("lheight") * par("cex") * cex.text
text(1, plotRegion[4] - 1.5 * lineHeight, labels = "rmse", adj = c(0.5, 0.5), cex = cex.text)
x.vals <- seq_len(length(methods))
y.vals <- plotRegion[4] - 0.5 * lineHeight
y.vals <- rep_len(y.vals, length(methods))
text(x.vals, y.vals, labels = round(rmse, 3),
     cex = cex.text,
     adj = c(0.5, 0.5))

if (plotIndividually) {
  dev.off()
  pdf(paste0("coverage_length", sim.setting, ".pdf"), 3.5, 3.5)
}

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

if (plotIndividually) {
  dev.off()
  pdf(paste0("pehe", sim.setting, ".pdf"), 3.5, 3.5)
}

boxplot(results[seq_len(maxIter),,"pehe"],
        main = "PEHE", ylab = "pehe",
        staplewex = 0)

if (plotIndividually) {
  dev.off()
  pdf(paste0("pegste", sim.setting, ".pdf"), 3.5, 3.5)
}

boxplot(results[seq_len(maxIter),,"pegste"],
        main = "PEGSTE", ylab = "pegste",
        staplewex = 0)
dev.off()

