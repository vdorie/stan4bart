# Part of the lme4 package for estimating generalized linear
# mixed effect models. Copyright (C) 2003-2021 Douglas Bates,
# Martin Maechler, Ben Bolker, and Steven Walker.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


# Part of the reformulas package for specifying formulas for lme4 
# models.

glFormula <- function(formula, data = NULL, subset, weights, 
    na.action, offset, contrasts = NULL,
    start, mustart, etastart, 
    control = glmerControl(), verbose = FALSE, ...) 
{
    control <- control$checkControl
    mf <- mc <- match.call()
    
    ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer", 
        "control", "nAGQ")
    l... <- list(...)
    l... <- l...[names(l...) %not_in% ignoreArgs]
    do.call(checkArgs, c(list("glmer"), l...))
    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr, control[[cstr]])
    denv <- checkFormulaData(formula, data,
                             checkLHS = control$check.formula.LHS == "stop")
    mc$formula <- formula <- as.formula(formula, env = denv)
    m <- match(c("data", "subset", "weights", "na.action", "offset", 
                 "mustart", "etastart"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf$na.action <- quote(stats::na.pass)
    
    # fr and fr.form contains everything, build specialized
    # frames below
    fr.form <- subbart(subbars(formula))
    environment(fr.form) <- environment(formula)
    for (i in c("weights", "offset")) {
      if (!eval(bquote(missing(x = .(i))))) 
        assign(i, get(i, parent.frame()), environment(fr.form))
    }
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    fr <- factorize(fr.form, fr, char.only = TRUE)
    attr(fr, "formula") <- formula
    attr(fr, "offset") <- mf$offset
    if (!missing(start) && is.list(start)) {
      attr(fr, "start") <- start$fixef
    }
    
    
    mf$na.action <- na.action

    getVariableNames <- function(terms) {
      # variables contains everything, even if there is a . - XX construction
      # Subtracted out terms don't appear in the term labels.
      
      factors <- attr(terms, "factors")
      if (NROW(factors) > 0L) result <- rownames(factors)[rowSums(factors) > 0L]
      else result <- attr(terms, "term.labels")

      variables <- as.character(attr(terms, "variables"))[-1L]
      response <- attr(terms, "response")
       
      if (length(response) > 0L) {
        response <- variables[response]
        setdiff(result, response)
      } else {
        result
      }
    }
    
    fixedform <- formula
    RHSForm(fixedform) <- nobart(nobars(RHSForm(fixedform)))
    RHSForm(mf$formula) <- RHSForm(fixedform)
    fixedfr <- eval(mf, parent.frame())
    fixedterms <- attr(fixedfr, "terms")
    attr(attr(fr, "terms"), "predvars.fixed") <- attr(fixedterms, "predvars")
    attr(attr(fr, "terms"), "varnames.fixed") <- getVariableNames(fixedterms)
    attr(fr, "na.action.fixed") <- attr(fixedfr, "na.action")
    
    
    bartform <- formula
    RHSForm(bartform) <- allbart(nobars(RHSForm(bartform)))
    RHSForm(mf$formula) <- RHSForm(bartform)
    mf$drop.unused.levels <- FALSE
    bartfr <- eval(mf, parent.frame())
    bartterms <- attr(bartfr, "terms")
    attr(attr(fr, "terms"), "predvars.bart") <- attr(bartterms, "predvars")
    attr(attr(fr, "terms"), "varnames.bart") <- getVariableNames(bartterms)
    bartlevels <-
      lapply(colnames(attr(bartterms, "factors")), function(n) levels(bartfr[[n]]))
    names(bartlevels) <- colnames(attr(bartterms, "factors"))
    attr(attr(fr, "terms"), "levels.bart") <- bartlevels
    attr(fr, "na.action.bart") <- attr(bartfr, "na.action")
    
    ranform <- formula
    RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))
    mf$formula <- ranform
    mf$drop.unused.levels <- TRUE
    ranfr <- eval(mf, parent.frame())
    ranterms <- attr(ranfr, "terms")
    attr(attr(fr, "terms"), "predvars.random") <- attr(ranterms, "predvars")
    attr(attr(fr, "terms"), "varnames.random") <- getVariableNames(ranterms)
    attr(fr, "na.action.random") <- attr(ranfr, "na.action")
    
    # We create a total frame with all of the data, and then subset the individual frames
    # based on the total na situation
    na.action.all <- c(attr(fixedfr, "na.action"), attr(bartfr, "na.action"), attr(ranfr, "na.action"))
    if (length(na.action.all) > 0L) {
      na.action.all <- sort(na.action.all[!duplicated(na.action.all)])
      if (!is.null(class(attr(fr, "na.action.fixed"))))
        class(na.action.all) <- class(attr(fr, "na.action.fixed"))
      else if (!is.null(class(attr(fr, "na.action.bart"))))
        class(na.action.all) <- class(attr(fr, "na.action.bart"))
      else
        class(na.action.all) <- class(attr(fr, "na.action.random"))
      
      attr(fr, "na.action.all") <- na.action.all
      
      all_rows <- seq_len(nrow(fr)) %not_in% na.action.all
      
      if (!is.null(attr(fixedfr, "na.action"))) {
        fixed_rows <- seq_len(nrow(fr)) %not_in% attr(fixedfr, "na.action")
        fixedfr <- fixedfr[all_rows[fixed_rows],,drop = FALSE]
      }
      
      if (!is.null(attr(bartfr, "na.action"))) {
        bart_rows <- seq_len(nrow(fr)) %not_in% attr(bartfr, "na.action")
        bartfr <- bartfr[all_rows[bart_rows],,drop = FALSE]
      }
      if (!is.null(attr(ranfr, "na.action"))) {
        random_rows <- seq_len(nrow(fr)) %not_in% attr(ranfr, "na.action")
        ranfr <- ranfr[all_rows[random_rows],,drop = FALSE]
      }
    }
    
    if (length(attr(ranterms, "term.labels")) > 0L) {
      reTrms <- mkReTrms(findbars(RHSForm(formula)), ranfr)
      
      wmsgNlev <- checkNlevels(reTrms$flist, n = nrow(ranfr), control, allow.n = TRUE)
      wmsgZdims <- checkZdims(reTrms$Ztlist, n = nrow(ranfr), control, allow.n = TRUE)
      wmsgZrank <- checkZrank(reTrms$Zt, n = nrow(ranfr), control, nonSmall = 1e+06, 
          allow.n = TRUE)
    } else {
      reTrms = list()
    }
    
    # dbartsData can't handle columns in the model frame not in the formula terms, so we
    # pull out  "(offset)" if it exists
    bartprednames <- as.character(attr(bartterms, "predvars"))[-1L]
    if (any(names(bartfr) %not_in% bartprednames)) {
      keepcols <- names(bartfr) %in% bartprednames
      bartfr <- bartfr[,keepcols, drop = FALSE]
      attr(bartfr, "terms") <- bartterms
      attr(attr(bartfr, "terms"), "dataClasses") <- attr(bartterms, "dataClasses")[keepcols]
      bartterms <- attr(bartfr, "terms")
      attr(attr(fr, "terms"), "varnames.bart") <- getVariableNames(bartterms)
    }
    
    bartData <- dbarts::dbartsData(bartform, bartfr)
    if (ncol(bartData@x) == 0L)
      stop("no bart component detected in formula; consider using rstanarm package instead")
    
    X <- model.matrix(fixedform, fixedfr, contrasts)
    # drop intercept, offset
    keepcols <- colnames(X) %not_in% c("(Intercept)", "`(offset)`")
    X <- X[,keepcols,drop=FALSE]
    if (is.null(rankX.chk <- control[["check.rankX"]])) 
      rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
    if (is.null(scaleX.chk <- control[["check.scaleX"]])) 
      scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind = scaleX.chk)
    
    terms <- attr(fr, "terms")
    
    varnames.random <- attr(terms, "varnames.random")
    varnames.fixed  <- attr(terms, "varnames.fixed")
    varnames.bart   <- attr(terms, "varnames.bart")
    
    varnames.mixed <- varnames.bart %in% varnames.random | varnames.bart %in% varnames.fixed
    if (any(varnames.mixed) && as.integer(verbose) > -1L) {
      warning("variable(s) '", paste0(varnames.bart[varnames.mixed], collapse = "', '"),
              "' that appear in both parametric and nonparametric are not identifiable; ",
              "model will fit but some results may be uninterpretable")
    }
    
    if (length(varnames.random) == 0L && length(varnames.fixed) == 0L)
      stop("no parametric component detected in formula; consider using dbarts package instead")
    
    result <- list(fr = fr, X = X, bartData = bartData, reTrms = reTrms, formula = formula, 
                   terms = terms)
    if (length(reTrms) > 0L)
      result$wmsgs <- c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)
    
    result
}

`%ORifNotInLme4%` <- tryCatch(
  {
    lme4ns <- getNamespace("lme4")
    function(a, b) {
      mc <- match.call()
      if (is.symbol(mc[["a"]])) a <- deparse(mc[["a"]])
      if (!is.null(lme4ns[[a]])) lme4ns[[a]] else b
    }
  },
  error = function(e) function(a, b) b
)

`%ORifNotInReformulas%` <- tryCatch(
  {
    reformulasns <- getNamespace("reformulas")
    function(a, b) {
      mc <- match.call()
      if (is.symbol(mc[["a"]])) a <- deparse(mc[["a"]])
      if (!is.null(reformulasns[[a]])) reformulasns[[a]] else b
    }
  },
  error = function(e) function(a, b) b
)

getMerControls <- function() {
  merControl <-
    function(optimizer="nloptwrap",
             restart_edge=TRUE,
             boundary.tol=1e-5,
             calc.derivs=NULL,
             use.last.params=FALSE,
             sparseX=FALSE,
             standardize.X=FALSE,
             autoscale=NULL,
             check.nobs.vs.rankZ="ignore",
             check.nobs.vs.nlev="stop",
             check.nlev.gtreq.5="ignore",
             check.nlev.gtr.1="stop",
             check.nobs.vs.nRE="stop",
             check.rankX = c("message+drop.cols",
                             "silent.drop.cols", "warn+drop.cols",
                             "stop.deficient", "ignore"),
             check.scaleX = c("warning","stop","silent.rescale",
                              "message+rescale","warn+rescale","ignore"),
             check.formula.LHS = "stop",
             check.conv.nobsmax = 1e4,
             check.conv.nparmax = 10,
             check.conv.grad     = .makeCC("warning", tol = 2e-3, relTol = NULL),
             check.conv.singular = .makeCC(action = "message", tol = getSingTol()),
             check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
             optCtrl = list(),
             mod.type="lmer"
             ) {

        stopifnot(is.list(optCtrl))
        if (mod.type=="glmer" && length(optimizer)==1) {
            optimizer <- replicate(2,optimizer)
        }
        c.opts <- paste0(mod.type,"Control")
        merOpts <- getOption(c.opts)
        if (!is.null(merOpts)) {
            nn <- names(merOpts)
            nn.ok <- .get.checkingOpts(names(merOpts))
            if (length(nn.ignored <- setdiff(nn,nn.ok))>0) {
                warning("some options in ",shQuote(sprintf("getOption('%s')",c.opts)),
                        " ignored : ",paste(nn.ignored,collapse=", "))
            }
            for (arg in nn.ok) {
                if (do.call(missing,list(arg)))
                    assign(arg,merOpts[[arg]])
            }
        }
        check.rankX <- match.arg(check.rankX)
        check.scaleX <- match.arg(check.scaleX)
        me <- sys.function()
        chk.cconv(check.conv.grad,     me)
        chk.cconv(check.conv.singular, me)
        chk.cconv(check.conv.hess    , me)

        if (mod.type=="glmer" && use.last.params && calc.derivs) {
            warning("using ",shQuote("use.last.params"),"=TRUE and ",
                    shQuote("calc.derivs"),"=TRUE with ",shQuote("glmer"),
                    " will not give backward-compatible results")
        }

        ret <- namedList(optimizer,
                         restart_edge,
                         boundary.tol,
                         calc.derivs,
                         use.last.params,
                         checkControl =
                             namedList(autoscale,
                                       check.nobs.vs.rankZ,
                                       check.nobs.vs.nlev,
                                       check.nlev.gtreq.5,
                                       check.nlev.gtr.1,
                                       check.nobs.vs.nRE,
                                       check.rankX,
                                       check.scaleX,
                                       check.formula.LHS),
                         checkConv=
                             namedList(check.conv.nobsmax,
                                       check.conv.nparmax,
                                       check.conv.grad,
                                       check.conv.singular,
                                       check.conv.hess),
                         optCtrl=optCtrl)
        if (mod.type=="glmer") {
            ret <- c(ret, namedList(tolPwrss,
                                    compDev,
                                    nAGQ0initStep))
            ret$checkControl <- c(ret$checkControl,
                                 namedList(check.response.not.const))
        }
        class(ret) <- c(c.opts, "merControl")
        ret
    }

  lmerControl <- merControl
  glmerControl <- merControl
  formals(glmerControl)[["check.conv.nparmax"]] <- 20
  formals(glmerControl)[["optimizer"]] <- c("bobyqa","Nelder_Mead")
  formals(glmerControl)[["mod.type"]] <- "glmer"
  formals(glmerControl)[["restart_edge"]] <- FALSE
  formals(glmerControl) <- c(formals(glmerControl),
                             list(tolPwrss=1e-7,
                                  compDev = TRUE,
                                  nAGQ0initStep = TRUE,
                                  check.response.not.const="stop")
                             )
  list(lmer = lmerControl, glmer = glmerControl)
}
merControls <- getMerControls()
rm(getMerControls)

glmerControl <- glmerControl %ORifNotInLme4% merControls$glmer

tolPwrss <- tolPwrss %ORifNotInLme4% NULL
compDev <- compDev %ORifNotInLme4% NULL
nAGQ0initStep <- nAGQ0initStep %ORifNotInLme4% NULL
check.response.not.const <- check.response.not.const %ORifNotInLme4% NULL

lmerControl <- merControls$lmer
rm(merControls)

checkArgs <- checkArgs %ORifNotInLme4% function (type, ...) 
{
  l... <- list(...)
  if (isTRUE(l...[["sparseX"]])) warning("sparseX = TRUE has no effect at present",call.=FALSE)
  if(length(l... <- list(...))) {
    if (!is.null(l...[["family"]])) {
      warning("calling lmer with family() is deprecated: please use glmer() instead",call.=FALSE)
      type <- "glmer"
    }
    if (!is.null(l...[["method"]])) {
      msg <- paste("Argument", sQuote("method"), "is deprecated.")
      if (type == "lmer")
        msg <- paste(msg, "Use the REML argument to specify ML or REML estimation.")
      else if (type == "glmer")
        msg <- paste(msg, "Use the nAGQ argument to specify Laplace (nAGQ=1) or adaptive",
                     "Gauss-Hermite quadrature (nAGQ>1).  PQL is no longer available.")
      warning(msg,call.=FALSE)
      l... <- l...[names(l...) != "method"]
    }
    if(length(l...)) {
      warning("extra argument(s) ",
              paste(sQuote(names(l...)), collapse=", "),
              " disregarded",call.=FALSE)
    }
  }
}

checkCtrlLevels <- checkCtrlLevels %ORifNotInLme4% function (cstr, val, smallOK=FALSE)
{
  bvals <- c("message","warning","stop","ignore")
  if (smallOK) bvals <- outer(bvals, c("","Small"), paste0)
  if (!is.null(val) && !val %in% bvals)
    stop("invalid control level ",sQuote(val)," in ",cstr,": valid options are {",
         paste(sapply(bvals,sQuote),collapse=","),"}")
  invisible(NULL)
}

checkFormulaData <- checkFormulaData %ORifNotInLme4% function (formula, data, checkLHS = TRUE, checkData = TRUE, debug = FALSE) 
{
    wd <- tryCatch(force(data), error = identity)
    if (bad.data <- inherits(wd,"error")) {
        bad.data.msg <- wd$message
    }
    if (bad.data || debug) {
        varex <- function(v, env) exists(v, envir=env, inherits=FALSE)
        allvars <- all.vars(as.formula(formula))
        allvarex <- function(env, vvec=allvars) all(vapply(vvec, varex, NA, env))
    }
    if (bad.data) {
        if (allvarex(environment(formula))) {
            stop("bad 'data', but variables found in environment of formula: ",
                 "try specifying 'formula' as a formula rather ",
                 "than a string in the original model",call.=FALSE)
        } else {
            stop("bad 'data': ", bad.data.msg, call. = FALSE)
        }
    } else {
        denv <-
            if (is.null(data)) {
                if (!is.null(ee <- environment(formula))) {
                    ee
            } else {
                parent.frame(2L)
            }
        } else
            list2env(data)
    }
    if (debug) {
        cat("Debugging parent frames in checkFormulaData:\n")
        glEnv <- 1L
        while (!identical(parent.frame(glEnv),.GlobalEnv)) {
            glEnv <- glEnv+1L
        }
        for (i in 1:glEnv) {
            OK <- allvarex(parent.frame(i))
            cat("vars exist in parent frame ", i)
            if (i == glEnv) cat(" (global)")
            cat(" ",OK, "\n")
        }
        cat("vars exist in env of formula ", allvarex(denv), "\n")
    }

    stopifnot(!checkLHS || length(as.formula(formula,env=denv)) == 3)
    return(denv)

}

mkReTrms <- mkReTrms %ORifNotInReformulas% function (bars, fr, drop.unused.levels=TRUE,
                     reorder.terms=TRUE,
                     reorder.vars=FALSE,
                     calc.lambdat = TRUE,
                     sparse = NULL)
{
  if (!length(bars))
    stop("No random effects terms specified in formula",call.=FALSE)
  stopifnot(is.list(bars), vapply(bars, is.language, NA),
            inherits(fr, "data.frame"))
  names(bars) <- barnames(bars)
  term.names <- vapply(bars, deparse1, "")
  blist <- lapply(bars, mkBlist, fr, drop.unused.levels,
                  reorder.vars = reorder.vars, sparse = sparse)
  nl <- vapply(blist, `[[`, 0L, "nl")
  
  ord <- seq_along(nl)
  if (reorder.terms) {
      if (any(diff(nl) > 0)) {
          ord <- rev(order(nl))
          blist      <- blist     [ord]
          nl         <- nl        [ord]
          term.names <- term.names[ord]
      }
  }
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(rbind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)

  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)

  if (calc.lambdat)
      nth <- as.integer((nc * (nc+1))/2)
      nb <- nc * nl

  if (sum(nb) != q) {
      stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)",
                   sum(nb),q))
  }
  boff <- cumsum(c(0L, nb))
  if (calc.lambdat) thoff <- cumsum(c(0L, nth))

  if (calc.lambdat) {
      mk_b <-function(i) {
          mm <- matrix(seq_len(nb[i]), ncol = nc[i],
                       byrow = TRUE)
          dd <- diag(nc[i])
          ltri <- lower.tri(dd, diag = TRUE)
          ii <- row(dd)[ltri]
          jj <- col(dd)[ltri]

          data.frame(i = as.vector(mm[, ii]) + boff[i],
                     j = as.vector(mm[, jj]) + boff[i],
                     x = as.double(rep.int(seq_along(ii),
                                           rep.int(nl[i], length(ii))) +
                                   thoff[i]))
      }
      Lambdat <- t(do.call(Matrix::sparseMatrix,
                           do.call(rbind,
                                   lapply(seq_along(blist), mk_b))))
      Lind <- as.integer(Lambdat@x)
  } else {
      Lambdat <- Lind <- NULL
  }
  thet <- NULL
  if (calc.lambdat)  thet <- numeric(sum(nth))
  ll <- list(Zt = drop0(Zt), theta = thet, Lind = Lind,
             Gp = unname(c(0L, cumsum(nb))))

  ll$lower <- -Inf * (thet + 1)
  if (calc.lambdat) {
      ll$lower[unique(diag(Lambdat))] <- 0
      Lambdat@x[] <- ll$theta[ll$Lind]
  }      
  ll$theta[] <- is.finite(ll$lower)
  ll$Lambdat <- Lambdat

  fl <- lapply(blist, `[[`, "ff")

  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else asgn <- seq_along(fl)
  names(fl) <- ufn

  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  ll$nl <- nl
  ll$ord <- ord
  ll
}

makeOp <- makeOp %ORifNotInReformulas% function(x, y, op=NULL)
{
  if (is.null(op) || missing(y)) {
    if (is.null(op)) {
      substitute(OP(X),list(X=x,OP=y))
    } else {
      substitute(OP(X),list(X=x,OP=op))
    }
  } else substitute(OP(X,Y), list(X=x,OP=op,Y=y))
}

expandDoubleVert <- expandDoubleVert %ORifNotInReformulas% function(term) {
  frml <- formula(substitute(~x,list(x=term[[2]])))
  tt <- terms(frml)
  newtrms <- lapply(attr(tt, "term.labels"),
                    function(t) {
                        sumTerms(list(0, toLang(t)))
                    })
  if(attr(tt, "intercept") != 0) {
    newtrms <- c(1, newtrms)
  }
  res <- lapply(newtrms,
         function(t) {
             makeOp(
                 makeOp(t, term[[3]], quote(`|`)),
                 quote(`(`)
             )
         })
  return(res)
}

sumTerms <- sumTerms %ORifNotInReformulas% function(termList) {
    Reduce(function(x,y) makeOp(x,y,op=quote(`+`)),termList)
}

anySpecial <- anySpecial %ORifNotInReformulas% function(term, specials=findReTrmClasses(), fast = FALSE)
{
  if (fast) return(any(specials %in% all.names(term)))
  has_s <- FALSE
  as2 <- function(expr) {
    if (length(expr) == 1) return(NULL)
    for (ss in specials) {
      if (identical(expr[[1]], as.name(ss))) {
        assign("has_s", TRUE, environment(as2))
        break
      }
    }
    if (has_s) return(NULL)
    lapply(expr[-1], as2)
  }
  as2(term)
  return(has_s)
}

inForm <- inForm%ORifNotInReformulas% function(form, value) {
    if (any(sapply(form,identical,value))) return(TRUE)
    if (all(sapply(form,length)==1)) return(FALSE)
    return(any(vapply(form,inForm,value,FUN.VALUE=logical(1))))
}


esfun <- esfun %ORifNotInReformulas% function(x) {
  if (length(x)==1 || !anySpecial(x, "|")) return(x)
  if (length(x)==2) {
    return(lapply(esfun(x[[2]]),  makeOp, y=x[[1]]))
  }
  if (length(x)==3) {
    if (x[[1]]==quote(`|`)) {
      return(lapply(expandGrpVar(x[[3]]),
                    makeOp, x=x[[2]], op=quote(`|`)))
    } else {
      return(x)
    }
  }
}

.valid_covstruct <- .valid_covstruct %ORifNotInReformulas% c(
  diag = 0,
  us   = 1,
  cs   = 2,
  ar1  = 3,
  ou   = 4,
  exp = 5,
  gau = 6,
  mat = 7,
  toep = 8,
  rr = 9,
  homdiag = 10
)

findReTrmClasses <- findReTrmClasses %ORifNotInReformulas% function() {
    c(names(.valid_covstruct), "s")
}

toLang <- toLang %ORifNotInReformulas% function(x) parse(text=x)[[1]]

expandGrpVar <- expandGrpVar %ORifNotInReformulas% function(f) {
  form <- as.formula(makeOp(f,quote(`~`)))
  mm <- terms(form)
  tl <- attr(mm,"term.labels")
  switch_order <- function(x) paste(rev(unlist(strsplit(x, ":"))), collapse = ":")
  if (inForm(f, quote(`/`))) {
    tl <- unname(vapply(tl, switch_order, character(1)))
    tl <- rev(tl)
  }
  res <- lapply(tl, toLang)
  return(res)
}


expandAllGrpVar <- expandAllGrpVar %ORifNotInReformulas% function(bb) {
  if (!is.list(bb))
    expandAllGrpVar(list(bb))
  else {
    for (i in seq_along(bb)) {
      return(unlist(lapply(bb,esfun)))
    }
  }
}

findbars_x <- findbars_x %ORifNotInReformulas% function(term,
                debug=FALSE,
                specials=character(0),
                default.special="us",
                target = '|',
                expand_doublevert_method = c("diag_special", "split"))
{
  expand_doublevert_method <- match.arg(expand_doublevert_method)

  if (length(term) == 3 && identical(term[[1]], quote(`~`))) {
    term <- RHSForm(term, as.form = TRUE)
  }
  ds <- if (is.null(default.special)) {
    NULL
  } else {
    eval(substitute(as.name(foo),list(foo=default.special)))
  }
  
  fbx <- function(term) {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (list(term[[1]]) %in% lapply(specials,as.name)) {
      if (debug) cat("special: ",deparse(term),"\n")
      return(term)
    }
    if (term[[1L]] == as.name(target)) {
      if (debug) {
        tt <- if (target == '|') "bar" else sprintf('"%s"', target)
        cat(sprintf("%s term: %s\n", tt, deparse(term)))
      }
      if (is.null(ds)) return(term)
      return(makeOp(term, ds))
    }
    if (term[[1L]] == as.name("||")) {
      if (expand_doublevert_method == "diag_special") {
          return(makeOp(makeOp(term[[2]], term[[3]],
                               op = quote(`|`)),
                        as.name("diag")))
      }
      return(lapply(expandDoubleVert(term), fbx))
    }
    if (term[[1L]] == as.name("(")) {
      if (debug) cat("paren term:",deparse(term),"\n")
      return(fbx(term[[2]]))
    }
    stopifnot(is.call(term))
    if (length(term) == 2) {
        ## unary operator, decompose argument
        if (debug) cat("unary operator:",deparse(term[[2]]),"\n")
        return(fbx(term[[2]]))
    }
    ## binary operator, decompose both arguments
    f2 <- fbx(term[[2]])
    f3 <- fbx(term[[3]])

    if (debug) { cat("binary operator:",deparse(term[[2]]),",",
                     deparse(term[[3]]),"\n")
                     cat("term 2: ", deparse(f2), "\n")
                     cat("term 3: ", deparse(f3), "\n")
    }
    c(f2, f3)
  }

  fbx_term <- fbx(term)
  if (debug) cat("fbx(term): ", deparse(fbx_term))
  expandAllGrpVar(fbx_term)
}

findbars <- findbars %ORifNotInReformulas% function (term) 
{
  findbars_x(term,
             default.special=NULL,
             expand_doublevert_method = "split")
}

RHSForm <- RHSForm %ORifNotInLme4% function (form, as.form = FALSE) 
{
  if (!as.form) return(form[[length(form)]])
  if (length(form)==2) return(form)
  form[[2]] <- NULL
  if (length(vars <- attr(form,"variables"))>0) {
    attr(form,"variables") <- vars[-2]
  }
  if (is.null(attr(form,"response"))) {
    attr(form,"response") <- 0
  }
  if (length(facs <- attr(form,"factors"))>0) {
    attr(form,"factors") <- facs[-1,]
  }
  return(form)
}

makeFac <- makeFac %ORifNotInReformulas% function(x,char.only=FALSE) {
    if (!is.factor(x) && (!char.only || is.character(x))) factor(x) else x
}

factorize <- factorize %ORifNotInReformulas% function (x, frloc, char.only = FALSE) 
{
  for (i in all.vars(RHSForm(x))) {
    if (!is.null(curf <- frloc[[i]]))
      frloc[[i]] <- makeFac(curf,char.only)
  }
  return(frloc)
}

checkNlevels <- checkNLevels %ORifNotInLme4% function (flist, n, ctrl, allow.n = FALSE) 
{
  stopifnot(is.list(ctrl), is.numeric(n))
  nlevelVec <- vapply(flist, function(x) nlevels(factor(x, exclude=NA)), 1)

  cstr <- "check.nlev.gtr.1"
  checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (doCheck(cc) && any(nlevelVec < 2)) {
      wstr <- "grouping factors must have > 1 sampled level"
      switch(cc,
             "warning" = warning(wstr,call.=FALSE),
             "stop"    =    stop(wstr,call.=FALSE),
             stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
  } else wstr <- character()

  cstr <- "check.nobs.vs.nlev"
  checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (doCheck(cc) && any(if(allow.n) nlevelVec > n else nlevelVec >= n)) {
      w <- if (allow.n) which(nlevelVec>n) else which(nlevelVec>=n)
      bad_facs <- names(nlevelVec)[w]
      wst2 <- gettextf(
          "number of levels of each grouping factor must be %s number of observations",
          if(allow.n) "<=" else "<")
      wst2 <- paste0(wst2," (problems: ",paste(bad_facs,collapse=", "),")")
      switch(cc,
             "warning" = warning(wst2, call.=FALSE),
             "stop"    =    stop(wst2, call.=FALSE)
             )
  } else wst2 <- character()

  cstr <- "check.nlev.gtreq.5"
  checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  if (doCheck(cc) && any(nlevelVec < 5)) {
      wst3 <- "grouping factors with < 5 sampled levels may give unreliable estimates"
      switch(cc,
             "warning" = warning(wst3, call.=FALSE),
             "stop"    = stop   (wst3, call.=FALSE),
             stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
  } else wst3 <- character()
  c(wstr, wst2, wst3)
}

checkZdims <- checkZdims %ORifNotInLme4% function (Ztlist, n, ctrl, allow.n = FALSE) 
{
  stopifnot(is.list(Ztlist), is.numeric(n))
  cstr <- "check.nobs.vs.nRE"
  checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
  term.names <- names(Ztlist)
  rows <- vapply(Ztlist, nrow, 1L)
  cols <- vapply(Ztlist, ncol, 1L)
  stopifnot(all(cols == n))
  if (doCheck(cc)) {
      unique(unlist(lapply(seq_along(Ztlist), function(i) {
          ww <- wmsg(cols[i], rows[i], allow.n, "number of observations",
                     "number of random effects",
                     sprintf(" for term (%s)", term.names[i]))
          if(ww$unident) {
              switch(cc,
                     "warning" = warning(ww$wstr, call.=FALSE),
                     "stop"    = stop   (ww$wstr, call.=FALSE),
                     stop(gettextf("unknown check level for '%s'", cstr), domain=NA))
              ww$wstr
          } else character()
      })))
  } else character()
}

checkZrank <- checkZrank %ORifNotInLme4% function(Zt, n, ctrl, nonSmall = 1e6, allow.n=FALSE)
{
  stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
  cstr <- "check.nobs.vs.rankZ"
  if (doCheck(cc <- ctrl[[cstr]])) {
    checkCtrlLevels(cstr, cc, smallOK=TRUE)
    d <- dim(Zt)
    doTr <- d[1L] < d[2L]
    if(!(grepl("Small",cc) && prod(d) > nonSmall)) {
      rankZ <- rankMatrix(if(doTr) t(Zt) else Zt, method="qr")
      ww <- wmsg(n, rankZ, allow.n, "number of observations", "rank(Z)")
      if(is.na(rankZ)) {
        cc <- "stop"
        ww <-
            list(unident = TRUE,
                 wstr = sub("^.*;",
                            "rank(Z) is NA: invalid random effect factors?",
                            ww$wstr))
      }
      if (ww$unident) {
        switch(cc,
               "warningSmall" =, "warning" = warning(ww$wstr,call.=FALSE),
               "stopSmall" =, "stop" = stop(ww$wstr,call.=FALSE),
               stop(gettextf("unknown check level for '%s'", cstr),
                    domain=NA))
        ww$wstr
      } else character()
    } else character()
  } else character()
}
 
subbart <- function(term)
{
  if (is.name(term) || !is.language(term)) return(term)
  
  if (length(term) == 2) {
    if (term[[1]] == as.name("bart"))
      term[[1]] <- as.name("(")
    term[[2]] <- subbart(term[[2]])
    return(term)
  }
  
  for (j in 2:length(term)) term[[j]] <- subbart(term[[j]])
  term
}

sub_specials <- sub_specials %ORifNotInReformulas% function (term,
                          specials = c("|", "||", "s"),
                        keep_args = c(2L, 2L, NA_integer_)) {
  if (is.name(term) || !is.language(term)) 
    return(term)
  for (i in seq_along(specials)) {
    if (is.call(term) && term[[1]] == as.name(specials[i])) {
      if (is.na(keep_args[i])) {
        if (!is.null(names(term))) {
          term <- term[names(term)==""]
        }
      } else {
        term <- term[1:(1+keep_args[i])]
      }
      term[[1]] <- as.name("+")
    }
  }
  for (j in 2:length(term)) {
    term[[j]] <- sub_specials(term[[j]],
                              specials = specials,
                              keep_args = keep_args)
  }
  term
}

subbars <- subbars %ORifNotInReformulas% function (term) 
  sub_specials(term, specials = c("|", "||"), keep_args = c(2L, 2L))

chkRank.drop.cols <- chkRank.drop.cols %ORifNotInLme4% function(X, kind, tol = 1e-7, method = "qr")
{
  stopifnot(is.matrix(X))
  kinds <- eval(formals(lmerControl)[["check.rankX"]])
  if (!kind %in% kinds) stop(gettextf("undefined option for 'kind': %s", kind))

  if(kind == "ignore") return(X)
  p <- ncol(X)
  if (kind == "stop.deficient") {
    if ((rX <- Matrix::rankMatrix(X, tol=tol, method=method)) < p)
      stop(gettextf(sub("\n +", "\n",
            "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);
             the fixed effects will be jointly unidentifiable"),
                    rX, p), call. = FALSE)
  } else {
    qr.X <- qr(X, tol = tol, LAPACK = FALSE)
    rnkX <- qr.X$rank
    if (rnkX == p)
        return(X)

    msg <- sprintf(ngettext(p - rnkX,
            "fixed-effect model matrix is rank deficient so dropping %d column / coefficient",
            "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"),
                   p - rnkX)
    if (kind != "silent.drop.cols")
      (if(kind == "warn+drop.cols") warning else message)(msg, domain = NA)
    contr <- attr(X, "contrasts")
    asgn <- attr(X, "assign")

    keep <- qr.X$pivot[seq_len(rnkX)]
    dropped.names <- colnames(X[,-keep,drop=FALSE])
    X <- X[, keep, drop = FALSE]
    if (Matrix::rankMatrix(X, tol=tol, method=method) < ncol(X))
        stop(gettextf("Dropping columns failed to produce full column rank design matrix"),
             call. = FALSE)

    if(!is.null(contr)) attr(X, "contrasts") <- contr
    if(!is.null(asgn))  attr(X, "assign")    <- asgn[keep]
    attr(X, "msgRankdrop") <- msg
    attr(X, "col.dropped") <- setNames(qr.X$pivot[(rnkX+1L):p],
                                       dropped.names)
  }
  X
}
checkScaleX <- checkScaleX %ORifNotInLme4% function(X,  kind="warning", tol=1e3) {
  kinds <- eval(formals(lmerControl)[["check.scaleX"]])
  if (!kind %in% kinds) stop(gettextf("unknown check-scale option: %s",kind))
  if (is.null(kind) || kind == "ignore") return(X)

  cont.cols <- apply(X,2,function(z) !all(z %in% c(0,1)))
  col.sd <- apply(X[,cont.cols, drop=FALSE], 2L, sd)
  sdcomp <- outer(col.sd,col.sd,"/")
  logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
  logsd <- abs(log(col.sd))
  if (any(c(logcomp,logsd) > log(tol))) {
    wmsg <- "Some predictor variables are on very different scales:"
    if (kind %in% c("warning","stop")) {
        msg2 <- "\nYou may also use (g)lmerControl(autoscale = TRUE) to improve numerical stability."
        wmsg <- paste(wmsg, "consider rescaling.", msg2)
        switch(kind,
               "warning" = warning(wmsg, call.=FALSE),
               "stop" = stop(wmsg, call.=FALSE))
    } else {
        wmsg <- paste(wmsg, "auto-rescaled (results NOT adjusted)")
        X[,cont.cols] <- sweep(X[,cont.cols,drop=FALSE],2,col.sd,"/")
        attr(X,"scaled:scale") <- setNames(col.sd,colnames(X)[cont.cols])
        if (kind == "warn+rescale") warning(wmsg, call.=FALSE)
    }
  } else
      wmsg <- character()
  structure(X, msgScaleX = wmsg)
}

if (getRversion() < "4.0.0") {
    deparse1 <- function (expr, collapse = " ", width.cutoff = 500L, ...) {
        paste(deparse(expr, width.cutoff, ...), collapse = collapse)
    }
}


barnames <- barnames %ORifNotInLme4% function (bars) 
  vapply(bars, function(x) deparse1(x[[3]]), "")

mkBlist <- mkBlist %ORifNotInReformulas% function(x,frloc, drop.unused.levels=TRUE,
                    reorder.vars=FALSE, sparse = NULL)
{
  frloc <- factorize(x,frloc)
  ff0 <- replaceTerm(x[[3]], quote(`:`), quote(`%i%`))
  ff <- try(eval(substitute(makeFac(fac),
                            list(fac = ff0)),
                 frloc), silent = TRUE)
  if (inherits(ff, "try-error")) {
    stop("couldn't evaluate grouping factor ",
         deparse1(x[[3]])," within model frame:",
         "error =",
         c(ff),
         " Try adding grouping factor to data ",
         "frame explicitly if possible",call.=FALSE)
  }
  if (all(is.na(ff)))
    stop("Invalid grouping factor specification, ",
         deparse1(x[[3]]),call.=FALSE)
  if (drop.unused.levels) ff <- factor(ff, exclude=NA)
  nl <- length(levels(ff))
  has.sparse.contrasts <- function(x) {
    cc <- attr(x, "contrasts")
    !is.null(cc) && inherits(cc, "sparseMatrix")
  }
  any.sparse.contrasts <- any(vapply(frloc, has.sparse.contrasts, FUN.VALUE = logical(1)))
  mMatrix <- if (!isTRUE(sparse) && !any.sparse.contrasts) model.matrix else Matrix::sparse.model.matrix
  mm <- mMatrix(eval(substitute( ~ foo, list(foo = x[[2]]))), frloc)
  if (reorder.vars) {
    mm <- mm[colSort(colnames(mm)),]
  }
  sm <- fac2sparse(ff, to = "d",
                   drop.unused.levels = drop.unused.levels)
  sm <- KhatriRao(sm, t(mm))
  dimnames(sm) <- list(
      rep(levels(ff),each=ncol(mm)),
      rownames(mm))
  list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

`%i%` <- `%i%` %ORifNotInReformulas% function(f1, f2, fix.order = TRUE) {
  if (!is.factor(f1) || !is.factor(f2)) stop("both inputs must be factors")
  f12 <- paste(f1, f2, sep = ":")
  u <- which(!duplicated(f12))
  if (!fix.order) return(factor(f12, levels = f12[u]))
  levs_rank <- length(levels(f2))*as.numeric(f1[u])+as.numeric(f2[u])
  return(factor(f12, levels = (f12[u])[order(levs_rank)]))
}

replaceTerm <- replaceTerm %ORifNotInReformulas% function(term,target,repl) {
  if (identical(term,target)) return(repl)
  if (!inForm(term,target)) return(term)
  if (length(term) == 2) {
    return(substitute(OP(x),list(OP=replaceTerm(term[[1]],target,repl),
                                 x=replaceTerm(term[[2]],target,repl))))
  }
  return(substitute(OP(x,y),list(OP=replaceTerm(term[[1]],target,repl),
                                 x=replaceTerm(term[[2]],target,repl),
                                 y=replaceTerm(term[[3]],target,repl))))
}

doCheck <- doCheck %ORifNotInLme4% function (x) 
{
  is.character(x) && !any(x == "ignore")
}

nobars_ <- nobars_ %ORifNotInReformulas% function (term) 
{
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
      nb <- nobars_(term[[2]])
      if(is.null(nb)) return(NULL)
      term[[2]] <- nb
      return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}


nobars <- nobars %ORifNotInReformulas% function (term) 
{
  e <- environment(term)
  nb <- nobars_(term)
  if (inherits(term,"formula") && length(term)==3 && is.symbol(nb)) {
    nb <- reformulate("1", response=deparse(nb))
  }
  if (is.null(nb)) {
    nb <- if (inherits(term,"formula")) ~1 else 1
  }
  environment(nb) <- e
  nb
}

anyBars <- anyBars  %ORifNotInLme4% function (term) 
{
  any(c("|", "||") %in% all.names(term))
}

isBar <- isBar %ORifNotInLme4% function(term)
{
  if (is.call(term)) {
    if ((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}

nobart_ <- function(term) 
{
  if (length(term) == 1L)
    return(if (term == as.name("bart")) NULL else term)
  
  if (length(term) == 2L) {
    if (term[[1]] == as.name("bart")) return(NULL)
    nb <- nobart_(term[[2]])
    if (is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  
  nb2 <- nobart_(term[[2]])
  nb3 <- nobart_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

nobart <- function(term)
{
  nb <- nobart_(term)
  if (inherits(term, "formula") && length(term) == 3 && is.symbol(nb)) {
    nb <- reformulate("1", response = deparse(nb))
  }
  if (is.null(nb)) {
    nb <- if (inherits(term, "formula")) 
      ~1
    else 1
  }
  nb
}

allbart_ <- function(term, inbart)
{
  if (length(term) == 1L)
    return(if (term == as.name("bart") || !inbart) NULL else term)
  
  if (length(term) == 2L) {
    if (term[[1]] == as.name("bart")) return(allbart_(term[[2L]], TRUE))
    ab <- allbart_(term[[2]], inbart)
    if (is.null(ab)) return(NULL)
    term[[2L]] <- ab
    return(term)
  }
  
  # should be a + b, or some other binary operator
  ab2 <- allbart_(term[[2L]], inbart)
  ab3 <- allbart_(term[[3L]], inbart)
  if (is.null(ab2)) return(ab3)
  if (is.null(ab3)) return(ab2)
  
  term[[2L]] <- ab2
  term[[3L]] <- ab3
  term
}

allbart <- function(term)
{
  ab <- allbart_(term, FALSE)
  if (inherits(term, "formula") && length(term) == 3 && is.symbol(ab)) {
    ab <- reformulate("1", response = deparse(ab))
  }
  if (is.null(ab)) {
    ab <- if (inherits(term, "formula")) 
      ~1
    else 1
  }
  ab
}

`RHSForm<-` <- `RHSForm<-` %ORifNotInReformulas% function (formula, value) 
{
    formula[[length(formula)]] <- value
    formula
}

reOnly <- function (f, response = FALSE) 
{
  found_bars <- findbars(f)
  if (is.null(found_bars)) return(~0)
  
  reformulate(paste0("(", vapply(findbars(f), deparse1, 
      ""), ")"), response = if (response && length(f) == 3L) 
      f[[2]])
}

mkVarCorr <- mkVarCorr %ORifNotInLme4% function (sc, cnms, nc, theta, nms) 
{
  ncseq <- seq_along(nc)
  thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
  if(!all(nms == names(cnms)))
      warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
              "Please report!", immediate. = TRUE)
  ans <- lapply(ncseq, function(i)
          {
              Li <- diag(nrow = nc[i])
              Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
              rownames(Li) <- cnms[[i]]
              val <- tcrossprod(sc * Li)
              stddev <- sqrt(diag(val))
              corr <- t(val / stddev)/stddev
              diag(corr) <- 1
              structure(val, stddev = stddev, correlation = corr)
          })
  if(is.character(nms)) {
    if (anyDuplicated(nms))
        nms <- make.names(nms, unique = TRUE)
    names(ans) <- nms
  }
  structure(ans, sc = sc)
    
}

wmsg <- wmsg %ORifNotInLme4% function (n, cmp.val, allow.n, msg1 = "", msg2 = "", msg3 = "") 
{
  if (allow.n) {
    unident <- n < cmp.val
    cmp <- "<"
    rstr <- ""
  }
  else {
    unident <- n <= cmp.val
    cmp <- "<="
    rstr <- " and the residual variance (or scale parameter)"
  }
  wstr <- sprintf("%s (=%d) %s %s (=%d)%s; the random-effects parameters%s are probably unidentifiable", 
      msg1, n, cmp, msg2, cmp.val, msg3, rstr)
  list(unident = unident, wstr = wstr)
}

.makeCC <- .makeCC %ORifNotInLme4% function (action, tol, relTol, ...) 
{
    stopifnot(is.character(action), length(action) == 1)
    if (!is.numeric(tol)) 
        stop("must have a numeric 'tol' component")
    if (length(tol) != 1 || tol < 0) 
        stop("'tol' must be number >= 0")
    mis.rt <- missing(relTol)
    if (!mis.rt && !is.null(relTol)) 
        stopifnot(is.numeric(relTol), length(relTol) == 1, relTol >= 
            0)
    c(list(action = action, tol = tol), if (!mis.rt) list(relTol = relTol), 
        list(...))
}

getSingTol <- getSingTol %ORifNotInLme4% function() {
  getOption("lme4.singular.tolerance", 1e-4)
}

isSingular <- isSingular %ORifNotInLme4% function (x, tol = getSingTol()) 
{
    lwr <- getME(x, "lower")
    theta <- getME(x, "theta")
    any(theta[lwr == 0] < tol)
}

.get.checkingOpts <- .get.checkingOpts %ORifNotInLme4% function (nms) 
  nms[grepl("^check\\.(?!conv|rankX|scaleX)", nms, perl = TRUE)]

chk.cconv <- chk.cconv %ORifNotInLme4% function (copt, callingFn) 
{
    cnm <- deparse(substitute(copt))
    if (is.character(copt)) {
        def <- eval(formals(callingFn)[[cnm]])
        def$action <- copt
        assign(cnm, def, envir = sys.frame(sys.parent()))
    }
    else chk.convOpt(copt)
}

namedList <- namedList %ORifNotInLme4% function (...) 
{
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L))) 
        nm <- snm
    if (any(nonames <- nm == "")) 
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

colSort <- colSort %ORifNotInLme4% function (x) 
{
    termlev <- vapply(strsplit(x, ":"), length, integer(1))
    iterms <- split(x, termlev)
    iterms <- sapply(iterms, sort, simplify = FALSE)
    ilab <- "(Intercept)"
    if (ilab %in% iterms[[1]]) {
        iterms[[1]] <- c(ilab, setdiff(iterms[[1]], ilab))
    }
    unlist(iterms)
}

isAnyArgBar <- isAnyArgBar %ORifNotInReformulas% function (term) 
{
    if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
        for (i in seq_along(term)) {
            if (isBar(term[[i]])) 
                return(TRUE)
        }
    }
    FALSE
}

chk.convOpt <- chk.convOpt %ORifNotInLme4% function (opt) 
{
    cnm <- deparse(nm <- substitute(opt))[[1]]
    if (!is.list(opt)) 
        stop("check.conv* option ", cnm, " must be a list")
    if (!is.character(opt$action)) 
        stop("check.conv* option ", cnm, " has no 'action' string")
    if (!is.numeric(tol <- opt$tol)) 
        stop("check.conv* option ", cnm, " must have a numeric 'tol' component")
    if (length(tol) != 1 || tol < 0) 
        stop("check.conv* option ", cnm, "$tol must be number >= 0")
    if (!is.null(relTol <- opt$relTol)) 
        stopifnot(is.numeric(relTol), length(relTol) == 1, relTol >= 
            0)
    invisible()
}

getME <- getME %ORifNotInLme4% function (object, name, ...) 
UseMethod("getME")

levelfun <- function (x, nl.n, sample_new_levels, Sigma) 
{
  old_levels <- dimnames(x)[["group"]]
  if (!all(nl.n %in% old_levels)) {
    nl.n.comb <- union(old_levels, nl.n)
    d <- dim(x)
    dn <- dimnames(x)
    dnn <- names(dn)
    # newx: num_predictors x num_levels x num_samples x num_chains
    newx <- array(0, c(d[1L], length(nl.n.comb), d[3L:4L]),
                  dimnames = list(dn[[1L]], nl.n.comb, dn[[3L]], dn[[4L]]))
    newx[,old_levels,,] <- x
    
    if (sample_new_levels) {
      new_levels <- nl.n.comb[!(nl.n.comb %in% old_levels)]
      L <- apply(Sigma, c(3L, 4L), function(x) t(base::chol(x)))
      n_predictors <- d[1L]
      n_groups  <- length(new_levels)
      n_samples <- d[3L]
      n_chains  <- d[4L]
      L <- array(Sigma, c(n_predictors, n_predictors, n_samples * n_chains))
      L <- bdiag(lapply(seq_len(n_samples * n_chains), function(i) t(base::chol(L[,,i,drop = FALSE]))))
      
      # L: block diagonal where each block is a sample of the covariance
      #    matrix for the random effects at that level
      #    (p x p) x (n_samp x n_chain)
      u <- matrix(rnorm(n_predictors * n_predictors * n_groups * n_samples * n_chains),
                  n_predictors * n_predictors * n_samples * n_chains,
                  n_groups)
      
      # L %*% u: (n_predictors x n_samp x n_chain) x n_groups
      newx[,new_levels,,] <- aperm(array(as.vector(Matrix::crossprod(L, u)), c(n_predictors, n_samples, n_chains, n_groups)), c(1L, 4L, 2L, 3L))
    }
    x <- newx
    names(dimnames(x)) <- dnn
    old_levels <- dimnames(x)[["group"]]
  }
  if (!all(r.inn <- old_levels %in% nl.n)) {
    x <- x[,r.inn,,,drop = FALSE]
  }
  x
}

# can't use lme4 version, as we need to use stan4bartFit terms and model.frame
get.orig.levs <- function (object, FUN = levels, newdata = NULL, sparse = FALSE, ...) 
{
  terms <- terms(object, ...)
  frame <- model.frame(object, ...)
  
  isFac <- vapply(frame, is.factor, FUN.VALUE = TRUE)
  isFac[attr(terms, "response")] <- FALSE
  frame <- frame[isFac]
  hasSparse <- any(grepl("sparse", names(formals(FUN))))
  orig_levs <- if (any(isFac) && hasSparse) lapply(frame, FUN, sparse = sparse) else if(any(isFac) && !hasSparse) lapply(frame, FUN)
  
  if (!is.null(newdata)) {      
    for (n in names(frame)) {
      orig_levs[[n]] <- c(orig_levs[[n]],
                          setdiff(unique(as.character(newdata[[n]])), orig_levs[[n]]))
    }
    
  }
  if (!is.null(orig_levs)) attr(orig_levs, "isFac") <- isFac
  orig_levs
}

formula.stan4bartFit <- function(x, ...)
{
  dots_list <- list(...)
  type <- match.arg(dots_list$type, c("all", "fixed", "random", "bart"))

  if (is.null(formula <- x$formula)) {
    if (!grepl("stan4bart$", deparse(x$call[[1]]))) 
      stop("can't find formula stored in model frame or call")
    form <- as.formula(formula(x$call))
  }
  
  if (type == "fixed") {
    RHSForm(formula) <- nobart(nobars(RHSForm(formula)))
  } else if (type == "random") {
    formula <- reOnly(formula, response = TRUE)
  } else if (type == "bart") {
    RHSForm(formula) <- allbart(nobars(RHSForm(formula)))
  }
  
  formula
}

terms.stan4bartFit <- function(x, ...)
{
  dots_list <- list(...)
  type <- match.arg(dots_list$type, c("all", "fixed", "random", "bart"))
  
  terms <- attr(x$frame, "terms")
  if (type == "all")
    return(terms)
  
  if (type == "fixed") {
    tt <- terms.formula(formula(x, type = "fixed"), data = x$frame)
    attr(tt, "predvars") <- attr(terms, "predvars.fixed")
  } else if (type == "random") {
    tt <- terms.formula(subbars(formula(x, type = "random")), data = x$frame)
    attr(tt, "predvars") <- attr(terms, "predvars.random")
  } else if (type == "bart") {
    tt <- terms.formula(formula(x, type = "bart"), data = x$frame)
    pred_vars <- attr(terms, "predvars.bart")
    attr(tt, "predvars") <- pred_vars
    
    # bart, when used with a . variable expansion, can end up
    # trying to build a frame with the offset variable. This
    # doesn't work with the dbarts model frame functions, so we
    # strip out anything that isn't in predvars
    term_labels <- attr(tt, "term.labels")
    extra_terms <- term_labels[term_labels %not_in% as.character(pred_vars)[-1L]]
    if (length(extra_terms) > 0L)
      tt <- strip_extra_terms(tt, extra_terms)
  }
  # Possibly re-order terms in case they apepared in a different order in
  # the original, which happens if, for example, a variable appears in the
  # fixed formula and again in the bart one (or in a . expression)
  
  vnames  <- as.character(attr(tt, "variables"))
  pvnames <- as.character(attr(tt, "predvars"))
  attr(tt, "variables")[vnames %in% pvnames] <- attr(tt, "variables")[match(pvnames, vnames)]
  
  tt
}

model.frame.stan4bartFit <- function(formula, ...)
{
  dots_list <- list(...)
  type <- match.arg(dots_list$type, c("all", "fixed", "random", "bart"))
  
  frame <- formula$frame
  
  if (type == "fixed") {
    # can return a lot of useless junk
    frame <- frame[as.character(attr(terms(formula), "predvars.fixed"))[-1L]]
  } else if (type == "random") {
    frame <- frame[as.character(attr(terms(formula), "predvars.random"))[-1L]]
  } else if (type == "bart") {
    frame <- frame[attr(terms(formula), "varnames.bart")]
  }

  frame
}

rm(`%ORifNotInReformulas%`)
if (exists("reformulasns")) rm(list = "reformulasns")
rm(`%ORifNotInLme4%`)
if (exists("lme4ns")) rm(list = "lme4ns")

make_glmerControl <- function(..., ignore_lhs = FALSE, ignore_x_scale = FALSE) {
  glmerControl(check.nlev.gtreq.5 = "ignore",
               check.nlev.gtr.1 = "stop",
               check.nobs.vs.rankZ = "ignore",
               check.nobs.vs.nlev = "ignore",
               check.nobs.vs.nRE = "ignore", 
               check.formula.LHS = if (ignore_lhs) "ignore" else "stop",
               check.scaleX = if (ignore_x_scale) "ignore" else "warning",
               ...)  
}

