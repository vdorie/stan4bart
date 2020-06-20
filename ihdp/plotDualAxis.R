plotBiasRMSE <- function(bias, methods, main, fileName = NULL)
{
  #collapse.bias <- function(x) if (NCOL(x) > 1L) apply(x, 1L, function(row) mean(row, na.rm = TRUE)) else x
  #collapse.r  <- function(x) if (NCOL(x) > 1L) apply(x, 1L, function(row) sqrt(mean(row^2, na.rm = TRUE))) else abs(x)
  collapse.bias <- function(x) apply(x, 2L, mean, na.rm = TRUE)
  collapse.r    <- function(x) apply(x, 2L, function(x) sqrt(mean(x^2, na.rm = TRUE)))
  
  #getIntervals.bias <- function(x) {
  #   if (NCOL(x) > 1L)
  #     list(lower = apply(x, 1L, function(row) quantile(row, 0.25, na.rm = TRUE)),
  #          upper = apply(x, 1L, function(row) quantile(row, 0.75, na.rm = TRUE)))
  #   else list(lower = x, upper = x)
  #}
  getIntervals.bias <- function(x) {
    res <- apply(x, 2L, quantile, na.rm = TRUE, probs = c(0.25, 0.75))
    list(lower = res[1L,], upper = res[2L,])
  }
  
  yRange <- c(-0.03, 0.04)
  yMinGrowth <- 0
  yMaxGrowth <- 0.25
  xGrowth <- 0.2
  
  ylab <- "bias"
  digits <- 2L
  abline <- quote(abline(h = 0, col = "gray"))
  labelPos <- "upperleft"
  
  
  if (!is.null(fileName)) {
    suffix <- sub(".*\\.([^.]+)$", "\\1", fileName, perl = TRUE)
    if (identical(suffix, "pdf")) {
      pdf(fileName, width = 3.5, height = 3.5)
    } else if (identical(suffix, "svg")) {
      svg(fileName, width = 3.5, height = 3.5, onefile = TRUE)
    } else if (identical(suffix, "png")) {
      png(fileName, width = 3.5, height = 3.5, units = "in", res = 144)
    }
  }
  
  parArgs <- list(mar = c(0.1, 2.5, 1.5, 2.5),
      mgp = c(1.25, 0.3, 0.0),
      tcl = -0.3,
      cex.main = 1.2)
  
  do.call("par", parArgs)
  
  
  yRange.orig <- yRange
  yMax <- yRange[2L] + 0.05 * (yRange[2L] - yRange[1L])
  yRange[1L] <- yRange[1L] - yMinGrowth * (yRange.orig[2L] - yRange.orig[1L])
  yRange[2L] <- yRange[2L] + yMaxGrowth * (yRange.orig[2L] - yRange.orig[1L])
  multiplier <- 10^digits
  
  xDelta.r <- 0.4
  
  xValues <- seq_along(methods)
  
  xRange <- range(xValues)
  xRange[1L] <- xRange[1L] - xDelta.r / 2 ## include space for moving 2nd set of plot symbols over
  xRange[2L] <- xRange[2L] + xDelta.r / 2
  xRange.orig <- xRange
  xRange[1L] <- xRange[1L] - xGrowth * diff(xRange.orig)
  xRange[2L] <- xRange[2L] + xGrowth * diff(xRange.orig)
  
  pointTypes <- c(1L, 20L, 2L, 17L)
  pointSizes <- c(0.875, 1.225, 0.75, 0.8)
  cex.axis <- 0.8
  
    x.i <- bias
    
    pts <- collapse.bias(x.i)
    if (!is.null(getIntervals.bias)) {
      intervals <- getIntervals.bias(x.i)
      lower <- intervals$lower
      upper <- intervals$upper
    }
    
    plot(NULL, type = "n", bty = "n",
         xlim = xRange,
         ylim = rev(yRange),
         main = main, xaxs = "i",
         yaxt = "n", xaxt = "n", ylab = ylab, xlab = "")
    
    minTick <- floor(multiplier * yRange.orig[1L]) / multiplier
    maxTick <- ceiling(multiplier * yRange.orig[2L]) / multiplier
    tickMarks <- round(axTicks(2, axp = c(minTick, maxTick, par("yaxp")[3L])), digits)
    if (!is.null(getIntervals.bias)) {
      axis(2, at = tickMarks, usr = range(lower, upper, na.rm = TRUE), cex.axis = cex.axis)
    } else {
      axis(2, at = tickMarks, usr = range(yRange.orig[1L], yRange.orig[2L], na.rm = TRUE), cex.axis = cex.axis)
    }
    
    if (!is.null(abline)) eval(abline)
    
    if (!is.null(getIntervals.bias)) {
      lines(rbind(xValues - xDelta.r / 2, xValues - xDelta.r / 2, rep_len(NA, length(pts))),
            rbind(lower, upper, rep_len(NA, length(pts))),
            lwd = 1)
    }
    
    points(xValues - xDelta.r / 2, pts, pch = pointTypes[1L], cex = pointSizes[1L])
    for (j in seq_along(methods))
      text(xValues[j], yMax, methods[j], srt = 90, cex = 0.85, adj = c(1, 0.5), font = 1L)
    
    #plotRegion <- par("usr")
    leftAxis <- c(maxTick, minTick)
    yValues.r <- collapse.r(x.i)
    yRange.r <- range(yValues.r)
    
    tickRange.r <- 0.85 * round((minTick - maxTick) * max(yValues.r) / min(pts), 2)
    maxTick.r <- tickRange.r * minTick / (minTick - maxTick)
    minTick.r <- tickRange.r * maxTick / (minTick - maxTick)
    
    points(xValues + xDelta.r / 2,
           diff(leftAxis) * (yValues.r - minTick.r) / (maxTick.r - minTick.r) + leftAxis[1L],
           pch = pointTypes[3L], cex = pointSizes[3L])
    
    tickMarks.r <- round(axTicks(2, axp = c(minTick.r, maxTick.r, par("yaxp")[3L])), digits)
    labels.r <- as.character(tickMarks.r)
    if (any(tickMarks.r < 0)) labels.r[tickMarks.r < 0] <- ""
    axis(4, at = diff(leftAxis) * (tickMarks.r - minTick.r) / (maxTick.r - minTick.r) + leftAxis[1L], labels = labels.r,
         cex.axis = cex.axis)
        
    plotRegion <- par("usr")
    #text(plotRegion[2L] + 1.6 * diff(grconvertX(0:1, "inches", "user")) * par("cin")[2L] * par("cex") * par("lheight"),
    #     mean(plotRegion[3L:4L]), "mean absolute deviation", cex = 0.9, srt = 90, xpd = TRUE)
    text(plotRegion[2L] + 1.6 * diff(grconvertX(0:1, "inches", "user")) * par("cin")[2L] * par("cex") * par("lheight"),
         mean(plotRegion[3L:4L]), "root mean squared error", cex = 0.9, srt = 90, xpd = TRUE)
    
    if (!is.null(labelPos)) {
      if (labelPos == "upperleft") {
        xRange <- par("usr")[c(1L, 2L)]
        x.label <- xRange[1L] + 0.05 * (xRange[2L] - xRange[1L])
        y.label <- yRange.orig[2L] - 0.05 * (yRange.orig[2L] - yRange.orig[1L])
        adj.label <- c(0, 0.5)
      } else if (labelPos == "upperright") {
        xRange <- par("usr")[c(1L, 2L)]
        x.label <- xRange[2L] - 0.05 * (xRange[2L] - xRange[1L])
        y.label <- yRange.orig[2L] - 0.05 * (yRange.orig[2L] - yRange.orig[1L])
        adj.label <- c(1, 0)
      }
    }
    
    #legend("topleft", legend = c("bias", "mad", "new"), pch = pointTypes[c(1L, 3L, 2L)], pt.cex = pointSizes[c(1L, 3L, 2L)], inset = 0.025, bty = "n")
    legend("topleft", legend = c("bias", "rmse"), cex = 0.65,
           pch = pointTypes[c(1L, 3L)], pt.cex = pointSizes[c(1L, 3L)],
           inset = 0.025, bty = "n")
  
  if (!is.null(fileName)) dev.off()
  
  invisible(NULL)
}


plotCovLenDualAxis <- function(cvrg, clen, methods, main, fileName = NULL)
{
  yRange.l <- c(0.75, 1)
  yRange.r <- c(.1, 0.2)
  
  yMinGrowth <- 0.25
  yMaxGrowth <- 0
  xGrowth <- 0.25
  
  x.l <- apply(cvrg, 2L, mean, na.rm = TRUE)
  x.r <- apply(clen, 2L, mean, na.rm = TRUE)
  ylab.l <- "coverage"
  ylab.r <- "ave interval length"
  
  tcks.l <- 5
  tcks.r <- 4
  
  digits <- 2L
  #abline <- quote(abline(h = 0, col = "gray"))
  #labelPos <- "bottomright"
  labelPos <- NULL
  
  if (!is.null(fileName)) {
    suffix <- sub(".*\\.([^.]+)$", "\\1", fileName, perl = TRUE)
    if (identical(suffix, "pdf")) {
      pdf(fileName, width = 3.5, height = 3.5)
    } else if (identical(suffix, "svg")) {
      svg(fileName, width = 3.5, height = 3.5, onefile = TRUE)
    } else if (identical(suffix, "png")) {
      png(fileName, width = 3.5, height = 3.5, units = "in", res = 144)
    }
  }
  
  parArgs <- list(mar = c(0.1, 2.5, 1.5, 2.5),
      mgp = c(1.25, 0.3, 0.0),
      tcl = -0.3,
      cex.main = 1.2)
  do.call("par", parArgs)
  
  
  yRange.orig <- yRange.l
  yMin.l <- min(yRange.l) - 0.05 * (max(yRange.l) - min(yRange.l))
  yMax.l <- max(yRange.l) + 0.05 * (max(yRange.l) - min(yRange.l))
  if (yRange.l[1L] > yRange.l[2L]) { temp <- yMax.l; yMax.l <- yMin.l; yMin.l <- temp }
  
  yRange.l[which.min(yRange.l)] <- yRange.l[which.min(yRange.l)] - yMinGrowth * (max(yRange.orig) - min(yRange.orig))
  yRange.l[which.max(yRange.l)] <- yRange.l[which.max(yRange.l)] + yMaxGrowth * (max(yRange.orig) - min(yRange.orig))
  
  multiplier <- 10^digits
  xDelta.r <- 0
  
  xValues <- seq_along(methods)
  xRange <- range(xValues)
  xRange[1L] <- xRange[1L] - xDelta.r / 2 ## include space for moving 2nd set of plot symbols over
  xRange[2L] <- xRange[2L] + xDelta.r / 2
  xRange.orig <- xRange
  xRange[1L] <- xRange[1L] - xGrowth * diff(xRange.orig)
  xRange[2L] <- xRange[2L] + xGrowth * diff(xRange.orig)
  
  pointTypes <- c(1L, 20L, 2L, 17L)
  pointSizes <- c(0.875, 1.225, 0.75, 0.8)
  cex.axis <- 0.8
  
    pts.l <- x.l
    plot(NULL, type = "n", bty = "n",
         xlim = xRange,
         ylim = yRange.l,
         main = main, xaxs = "i",
         yaxt = "n", xaxt = "n", ylab = ylab.l, xlab = "")
    
    minTick.l <- floor(multiplier * min(yRange.orig)) / multiplier
    maxTick.l <- ceiling(multiplier * max(yRange.orig)) / multiplier
    tickMarks <- round(axTicks(2, axp = c(minTick.l, maxTick.l, tcks.l)), digits)
    axis(2, at = tickMarks, usr = range(yRange.orig[1L], yRange.orig[2L], na.rm = TRUE), cex.axis = cex.axis)
    
    # abline(h = 0.95, col = "gray")
    truncateToRange <- function(x, r) {
      u <- x > max(r); l <- x < min(r)
      if (any(u, na.rm = TRUE)) x[u] <- max(r)
      if (any(l, na.rm = TRUE)) x[l] <- min(r)
      list(x = x, cut = u | l)
    }
    
    y.l <- pts.l
    temp <- truncateToRange(y.l, yRange.orig)
    y.l.cut <- temp$x
    cut.l <- temp$cut
    for (j in seq_along(methods)) {
      text(xValues[j], yMin.l, methods[j], srt = 90, cex = 0.85, adj = c(1, 0.5), font = 1L)
      #lines(rbind(xValues, xValues, NA), rep(c(min(yRange.orig), max(yRange.orig), NA), length(methods)), lwd = 0.25, col = rgb(0.8, 0.8, 0.8))
    }
    points(xValues - xDelta.r / 2, y.l.cut, pch = pointTypes[1L], cex = pointSizes[1L],
           col = ifelse(cut.l, "gray", "black"))
    
    axis.l <- c(minTick.l, maxTick.l)
    yValues.r <- x.r
    
    minTick.r <- min(yRange.r)
    maxTick.r <- max(yRange.r)
    axis.r <- c(minTick.r, maxTick.r)
    
    y.r <- yValues.r
    temp <- truncateToRange(y.r, yRange.r)
    y.r.cut <- temp$x
    cut.r <- temp$cut
    
    coordConvert <- function(x)
      abs(diff(axis.l)) * (x - min(axis.r)) / (max(axis.r) - min(axis.r)) + min(axis.l)
    
    y.r.cut <- coordConvert(y.r.cut)
    #axis.r <- abs(diff(axis.l)) * (axis.r - min(axis.r)) / abs(diff(axis.r)) + min(axis.l)
    
    points(xValues + xDelta.r / 2, y.r.cut,
           pch = pointTypes[3L], cex = pointSizes[3L],
           col = ifelse(cut.r, "gray", "black"))
    
    #tickMarks.r <- round(axTicks(2, axp = c(maxTick.r, minTick.r, par("yaxp")[3L])), digits)
    tickMarks.r <- round(axTicks(2, axp = c(axis.r[1L], axis.r[2L], tcks.r)), digits)
    labels.r <- as.character(tickMarks.r)
    if (any(tickMarks.r < 0)) labels.r[tickMarks.r < 0] <- ""
    axis(4, at = coordConvert(tickMarks.r), labels = labels.r, cex.axis = cex.axis)
    
    #abline(h = coordConvert(0.95), col = "gray")
    abline(h = 0.95, col = "gray")
    
    plotRegion <- par("usr")
    #text(plotRegion[2L] + 1.6 * diff(grconvertX(0:1, "inches", "user")) * par("cin")[2L] * par("cex") * par("lheight"),
    #     mean(plotRegion[3L:4L]), "mean absolute deviation", cex = 0.9, srt = 90, xpd = TRUE)
    text(plotRegion[2L] + 1.6 * diff(grconvertX(0:1, "inches", "user")) * par("cin")[2L] * par("cex") * par("lheight"),
         mean(plotRegion[3L:4L]), ylab.r, cex = 0.9, srt = 90, xpd = TRUE)
    
    if (!is.null(labelPos)) {
      xRange <- par("usr")[c(1L, 2L)]
      if (labelPos == "upperleft") {
        x.label <- xRange[1L] + 0.05 * (xRange[2L] - xRange[1L])
        y.label <- yRange.orig[2L] - 0.05 * (yRange.orig[2L] - yRange.orig[1L])
        adj.label <- c(0, 0.5)
      } else if (labelPos == "upperright") {
        x.label <- xRange[2L] - 0.05 * (xRange[2L] - xRange[1L])
        y.label <- yRange.orig[2L] - 0.05 * (yRange.orig[2L] - yRange.orig[1L])
        adj.label <- c(1, 0)
      } else if (labelPos == "bottomright") {
        x.label <- xRange[2L] - 0.05 * (xRange[2L] - xRange[1L])
        y.label <- yRange.orig[1L] + 0.05 * (yRange.orig[2L] - yRange.orig[1L])
        adj.label <- c(1, 0.5)
      }
      text(x.label, y.label, paste0("k = ", sum(subsets[[i]]), ", r = 100"),
           adj = adj.label, cex = 1.1)
    }
    
    y.legend <- 0.025 * (max(yRange.orig) - min(yRange.orig)) + min(yRange.orig)
    # x.legend <- max(xRange) - 0.025 * (max(xRange) - min(xRange))
     x.legend <- min(xRange) + 0.025 * (max(xRange) - min(xRange))
    legend(x.legend, y.legend, xjust = 0, yjust = 0, legend = c("coverage", "int len"), cex = 0.75, pch = pointTypes[c(1L, 3L)], pt.cex = pointSizes[c(1L, 3L)], inset = 0.025, bty = "o", bg = "white", box.lwd = NA)
  
  if (!is.null(fileName)) dev.off()
  
  invisible(NULL)
}

