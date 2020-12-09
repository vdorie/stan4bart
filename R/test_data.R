getTestDataFrames <- function(object, newdata, na.action = na.pass,
                              type = c("all", "random", "fixed", "bart"))
{
  type <- match.arg(type)
  
  result <- list()
  
  fixed.na.action <- NULL
  terms.fixed <- terms(object, type = "fixed")
  if (!identical(na.action, na.pass)) {
    mfnew <- model.frame(delete.response(terms.fixed), newdata, na.action = na.action)
    fixed.na.action <- attr(mfnew, "na.action")
  }
  nobs <- nrow(newdata)
  
  if (type %in% c("all", "fixed") && !is.null(object$X)) {
    X <- object$X
    X.col.dropped <- attr(X, "col.dropped")
    RHS <- formula(substitute(~R, list(R = RHSForm(formula(object, type = "fixed")))))
    
    orig.fixed.levs <- get.orig.levs(object, type = "fixed")
    
    mfnew <- suppressWarnings(
      model.frame(delete.response(terms.fixed), newdata,
                  na.action = fixed.na.action, xlev = orig.fixed.levs)
    )
    X <- model.matrix(RHS, data = mfnew, contrasts.arg = attr(X, "contrasts"))
    if (is.numeric(X.col.dropped) && length(X.col.dropped) > 0)
      X <- X[, -X.col.dropped, drop = FALSE]
    result$X <- X
  }
  if (type %in% c("all", "bart")) {
    orig.bart.levs <- get.orig.levs(object, type = "bart")
    
    mf.bart <- suppressWarnings(
      model.frame.default(delete.response(terms(object, type = "bart")), newdata,
                  na.action = fixed.na.action, xlev = orig.bart.levs)
    )

    result$X.bart <- dbarts::makeTestModelMatrix(object$bartData, mf.bart)
  }
  if (type %in% c("all", "random") && !is.null(object$reTrms)) {
    terms.random <- terms(object, type = "random")
    form.random <- formula(object, type = "random")
    
    newdata.NA <- newdata
    if (!is.null(fixed.na.action)) {
      newdata.NA <- newdata.NA[-fixed.na.action,]
    }
    
    tt <- delete.response(terms.random)
    frame.random <- model.frame(object, type = "random")
    orig.random.levs <- get.orig.levs(object, newdata = newdata.NA, type = "random")
    sparse <- !is.null(orig.random.levs) && max(lengths(orig.random.levs)) > 100
    orig.random.cntr <- get.orig.levs(object, FUN = contrasts, sparse = sparse, type = "random")
    
    re.form <- reOnly(object$formula)
    
    pv <- attr(tt, "predvars")
    for (i in 2:(length(pv))) {
      missvars <- setdiff(all.vars(pv[[i]]), all.vars(re.form))
      for (mv in missvars)
        newdata.NA[[mv]] <- NA
    }
    
    rfd <- suppressWarnings(model.frame(tt, newdata.NA, na.action = na.pass, 
                                        xlev = orig.random.levs))
    termvars <- unique(unlist(lapply(findbars(form.random), function(x) all.vars(x[[2]]))))
    for (fn in Reduce(intersect, list(names(orig.random.cntr), termvars, names(rfd)))) {
      if (!is.factor(rfd[[fn]])) 
        rfd[[fn]] <- factor(rfd[[fn]])
      contrasts(rfd[[fn]]) <- orig.random.cntr[[fn]]
    }
    if (!is.null(fixed.na.action)) 
      attr(rfd, "na.action") <- fixed.na.action
    
        
    ReTrms.test <- mkReTrms(findbars(re.form[[2]]), rfd)
    
    attr(ReTrms.test$Zt, "na.action") <- fixed.na.action
    
    result$reTrms <- list(Zt      = ReTrms.test$Zt,
                          Lambdat = ReTrms.test$Lambdat,
                          flist   = ReTrms.test$flist,
                          cnms    = ReTrms.test$cnms)
  }
  result
}


