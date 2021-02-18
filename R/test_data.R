getTestDataFrames <- function(object, newdata, na.action = na.pass,
                              type = c("all", "random", "fixed", "bart"))
{
  type <- match.arg(type)
  
  # Create a single frame that has data for the whole formula
  formula <- subbart(subbars(formula(object)))
  formula[[2L]] <- NULL
  environment(formula) <- environment()
  
  mf_call <- quote(stats::model.frame(formula = formula, data = newdata, na.action = "na.pass"))
  
  result <- list(frame = eval(mf_call))
  
  # define the sub-model frames as applicable
  if (type %in% c("all", "fixed") && !is.null(object$X)) {
    orig.fixed.levs <- get.orig.levs(object, type = "fixed")
    
    mf.fixed <- suppressWarnings(
      model.frame(delete.response(terms(object, type = "fixed")), newdata,
                  na.action = na.action, xlev = orig.fixed.levs)
    )
        
    if (!identical(na.action, na.pass))
      result$na.action.fixed <- attr(mf.fixed, "na.action")
  }
  
  if (type %in% c("all", "bart")) {
    orig.bart.levs <- attr(terms(object), "levels.bart")
    
    mf.bart <- suppressWarnings(
      model.frame(delete.response(terms(object, type = "bart")), newdata,
                  na.action = na.action, xlev = orig.bart.levs)
    )
    
    if (!identical(na.action, na.pass))
      result$na.action.bart <- attr(mf.bart, "na.action")
  }
  
  if (type %in% c("all", "random") && !is.null(object$reTrms)) {
    terms.random <- terms(object, type = "random")
    form.random <- formula(object, type = "random")
    
    tt <- delete.response(terms.random)
    frame.random <- model.frame(object, type = "random")
    orig.random.levs <- get.orig.levs(object, newdata = newdata, type = "random")
    sparse <- !is.null(orig.random.levs) && max(lengths(orig.random.levs)) > 100
    orig.random.cntr <- get.orig.levs(object, FUN = contrasts, sparse = sparse, type = "random")
    
    re.form <- reOnly(object$formula)
    
    newdata.random <- newdata
    
    pv <- attr(tt, "predvars")
    if (length(pv) > 1L) for (i in seq.int(2L, length(pv))) {
      missvars <- setdiff(all.vars(pv[[i]]), all.vars(re.form))
      for (mv in missvars)
        newdata.random[[mv]] <- NA
    }
    
    mf.random <- suppressWarnings(model.frame(tt, newdata.random, na.action = na.action, 
                                              xlev = orig.random.levs))
    termvars <- unique(unlist(lapply(findbars(form.random), function(x) all.vars(x[[2]]))))
    for (fn in Reduce(intersect, list(names(orig.random.cntr), termvars, names(mf.random)))) {
      if (!is.factor(mf.random[[fn]])) 
        mf.random[[fn]] <- factor(mf.random[[fn]])
      contrasts(mf.random[[fn]]) <- orig.random.cntr[[fn]]
    }
    if (!identical(na.action, na.pass))
      result$na.action.random <- attr(mf.random, "na.action") 
  }
  
  na.action.all <- c(result$na.action.fixed, result$na.action.random, result$na.action.bart)
  if (length(na.action.all) > 0L) {
    na.action.all <- sort(na.action.all[!duplicated(na.action.all)])
    if (!is.null(class(result$na.action.fixed)))
      class(na.action.all) <- class(result$na.action.fixed)
    else if (!is.null(class(result$na.action.bart)))
      class(na.action.all) <- class(result$na.action.bart)
    else
      class(na.action.all) <- class(result$na.action.random)
    
    result$na.action.all <- na.action.all
    
    all_rows <- seq_len(nrow(result$frame)) %not_in% result$na.action.all
  }
  
  # If na is omit, subset the model frames down to just their shared 
  # complete cases.
  if (exists("mf.fixed")) {
    if (inherits(na.action.all, "omit")) {
      fixed_rows <- seq_len(nrow(result$frame)) %not_in% (result$na.action.fixed %ORifNULL% integer(0L))
      mf.fixed <- mf.fixed[all_rows[fixed_rows],,drop = FALSE]
    }
    
    rhs.fixed <- formula(substitute(~R, list(R = RHSForm(formula(object, type = "fixed")))))
    X.col.dropped <- attr(object$X, "col.dropped")
    X <- model.matrix(rhs.fixed, data = mf.fixed, contrasts.arg = attr(object$X, "contrasts"))
    if (is.numeric(X.col.dropped) && length(X.col.dropped) > 0)
      X <- X[, -X.col.dropped, drop = FALSE]
    result$X <- X

  }
  if (exists("mf.bart")) {
    if (inherits(na.action.all, "omit")) {
      bart_rows <- seq_len(nrow(result$frame)) %not_in% (result$na.action.bart %ORifNULL% integer(0L))
      mf.bart <- mf.bart[all_rows[bart_rows],,drop = FALSE]
    }
    
    result$X.bart <- dbarts::makeTestModelMatrix(object$bartData, mf.bart)
  }
  
  
  if (exists("mf.random")) {
    if (inherits(na.action.all, "omit")) {
      random_rows <- seq_len(nrow(result$frame)) %not_in% (result$na.action.random %ORifNULL% integer(0L))
      mf.random <- mf.random[all_rows[random_rows],,drop = FALSE]
    }
    
    ReTrms.test <- mkReTrms(findbars(re.form[[2]]), mf.random)
    
    if (inherits(result$na.action.random, "omit")) {
      attr(ReTrms.test$Zt, "na.action") <- result$na.action.all
    } else {
      attr(ReTrms.test$Zt, "na.action") <- attr(mf.random, "na.action")
    }
    
    result$reTrms <- list(Zt      = ReTrms.test$Zt,
                          Lambdat = ReTrms.test$Lambdat,
                          flist   = ReTrms.test$flist,
                          cnms    = ReTrms.test$cnms)
  }
  
  result
}


