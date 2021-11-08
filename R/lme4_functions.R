glFormula <- function (formula, data = NULL, subset, weights, 
    na.action, offset, contrasts = NULL,
    start, mustart, etastart, 
    control = glmerControl(), ...) 
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
    
    fixedform <- formula
    RHSForm(fixedform) <- nobart(nobars(RHSForm(fixedform)))
    RHSForm(mf$formula) <- RHSForm(fixedform)
    fixedfr <- eval(mf, parent.frame())
    fixedterms <- attr(fixedfr, "terms")
    attr(attr(fr, "terms"), "predvars.fixed") <- attr(fixedterms, "predvars")
    attr(attr(fr, "terms"), "varnames.fixed") <- names(fixedfr)
    attr(fr, "na.action.fixed") <- attr(fixedfr, "na.action")
    
    
    bartform <- formula
    RHSForm(bartform) <- allbart(nobars(RHSForm(bartform)))
    RHSForm(mf$formula) <- RHSForm(bartform)
    mf$drop.unused.levels <- FALSE
    bartfr <- eval(mf, parent.frame())
    bartterms <- attr(bartfr, "terms")
    attr(attr(fr, "terms"), "predvars.bart") <- attr(bartterms, "predvars")
    attr(attr(fr, "terms"), "varnames.bart") <- names(bartfr)
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
    attr(attr(fr, "terms"), "varnames.random") <- names(ranfr)
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
      attr(attr(fr, "terms"), "varnames.bart") <- names(bartfr)
    }
    
    bartData <- dbarts::dbartsData(bartform, bartfr)
    if (ncol(bartData@x) == 0L)
      stop("no bart component detected; consider using rstanarm instead")
    
    X <- model.matrix(fixedform, fixedfr, contrasts)
    if (is.null(rankX.chk <- control[["check.rankX"]])) 
      rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
    if (is.null(scaleX.chk <- control[["check.scaleX"]])) 
      scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind = scaleX.chk)
    intercept_col <- match("(Intercept)", colnames(X))
    if (!is.na(intercept_col))
      X <- X[, -intercept_col, drop = FALSE]
    
    terms <- attr(fr, "terms")
    
    varnames.random <- attr(ranterms,   "term.labels")[attr(ranterms,   "order") == 1]
    varnames.fixed  <- attr(fixedterms, "term.labels")[attr(fixedterms, "order") == 1]
    varnames.bart   <- attr(bartterms,  "term.labels")[attr(bartterms,  "order") == 1]
    
    varnames.mixed <- varnames.bart %in% varnames.random | varnames.bart %in% varnames.fixed
    if (any(varnames.mixed))
      warning("variable(s) '", paste0(varnames.bart[varnames.mixed], collapse = "', '"),
              "' that appear in both parametric and nonparametric are not identifiable; ",
              "model will fit but some results may be uninterpretable")
    
    result <- list(fr = fr, X = X, bartData = bartData, reTrms = reTrms, formula = formula, 
                   terms = terms)
    if (length(reTrms) > 0L)
      result$wmsgs <- c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank)
    
    result
}

tryResult <- tryCatch(lme4ns <- asNamespace("lme4"), error = function(e) e)
if (!inherits(tryResult, "error")) {
  `%ORifNotInLme4%` <- function(a, b) {
    mc <- match.call()
    if (is.symbol(mc[["a"]])) a <- deparse(mc[["a"]])
    if (!is.null(lme4ns[[a]])) lme4ns[[a]] else b
  }
} else {
  `%ORifNotInLme4%` <- function(a, b) b
}

glmerControl <- glmerControl %ORifNotInLme4% function (optimizer = c("bobyqa", "Nelder_Mead"), restart_edge = FALSE, 
    boundary.tol = 1e-05, calc.derivs = TRUE, use.last.params = FALSE, 
    sparseX = FALSE, standardize.X = FALSE, check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "stop", check.nlev.gtreq.5 = "ignore", 
    check.nlev.gtr.1 = "stop", check.nobs.vs.nRE = "stop", check.rankX = c("message+drop.cols", 
        "silent.drop.cols", "warn+drop.cols", "stop.deficient", 
        "ignore"), check.scaleX = c("warning", "stop", "silent.rescale", 
        "message+rescale", "warn+rescale", "ignore"), check.formula.LHS = "stop", 
    check.conv.grad = .makeCC("warning", tol = 0.002, relTol = NULL), 
    check.conv.singular = .makeCC(action = "message", tol = formals(isSingular)$tol), 
    check.conv.hess = .makeCC(action = "warning", tol = 1e-06), 
    optCtrl = list(), mod.type = "glmer", tolPwrss = 1e-07, compDev = TRUE, 
    nAGQ0initStep = TRUE, check.response.not.const = "stop") 
{
    stopifnot(is.list(optCtrl))
    if (mod.type == "glmer" && length(optimizer) == 1) {
        optimizer <- replicate(2, optimizer)
    }
    c.opts <- paste0(mod.type, "Control")
    merOpts <- getOption(c.opts)
    if (!is.null(merOpts)) {
        nn <- names(merOpts)
        nn.ok <- .get.checkingOpts(names(merOpts))
        if (length(nn.ignored <- setdiff(nn, nn.ok)) > 0) {
            warning("some options in ", shQuote(sprintf("getOption('%s')", 
                c.opts)), " ignored : ", paste(nn.ignored, collapse = ", "))
        }
        for (arg in nn.ok) {
            if (do.call(missing, list(arg))) 
                assign(arg, merOpts[[arg]])
        }
    }
    check.rankX <- match.arg(check.rankX)
    check.scaleX <- match.arg(check.scaleX)
    me <- sys.function()
    chk.cconv(check.conv.grad, me)
    chk.cconv(check.conv.singular, me)
    chk.cconv(check.conv.hess, me)
    if (mod.type == "glmer" && use.last.params && calc.derivs) {
        warning("using ", shQuote("use.last.params"), "=TRUE and ", 
            shQuote("calc.derivs"), "=TRUE with ", shQuote("glmer"), 
            " will not give backward-compatible results")
    }
    ret <- namedList(optimizer, restart_edge, boundary.tol, calc.derivs, 
        use.last.params, checkControl = namedList(check.nobs.vs.rankZ, 
            check.nobs.vs.nlev, check.nlev.gtreq.5, check.nlev.gtr.1, 
            check.nobs.vs.nRE, check.rankX, check.scaleX, check.formula.LHS), 
        checkConv = namedList(check.conv.grad, check.conv.singular, 
            check.conv.hess), optCtrl = optCtrl)
    if (mod.type == "glmer") {
        ret <- c(ret, namedList(tolPwrss, compDev, nAGQ0initStep))
        ret$checkControl <- c(ret$checkControl, namedList(check.response.not.const))
    }
    class(ret) <- c(c.opts, "merControl")
    ret
}

tolPwrss <- tolPwrss %ORifNotInLme4% NULL
compDev <- compDev %ORifNotInLme4% NULL
nAGQ0initStep <- nAGQ0initStep %ORifNotInLme4% NULL
check.response.not.const <- check.response.not.const %ORifNotInLme4% NULL

lmerControl <- lmerControl %ORifNotInLme4%
  function (optimizer = "nloptwrap", restart_edge = TRUE, boundary.tol = 1e-05, 
    calc.derivs = TRUE, use.last.params = FALSE, sparseX = FALSE, 
    standardize.X = FALSE, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "stop", 
    check.nlev.gtreq.5 = "ignore", check.nlev.gtr.1 = "stop", 
    check.nobs.vs.nRE = "stop", check.rankX = c("message+drop.cols", 
        "silent.drop.cols", "warn+drop.cols", "stop.deficient", 
        "ignore"), check.scaleX = c("warning", "stop", "silent.rescale", 
        "message+rescale", "warn+rescale", "ignore"), check.formula.LHS = "stop", 
    check.conv.grad = .makeCC("warning", tol = 0.002, relTol = NULL), 
    check.conv.singular = .makeCC(action = "message", tol = formals(isSingular)$tol), 
    check.conv.hess = .makeCC(action = "warning", tol = 1e-06), 
    optCtrl = list(), mod.type = "lmer") 
{
    stopifnot(is.list(optCtrl))
    if (mod.type == "glmer" && length(optimizer) == 1) {
        optimizer <- replicate(2, optimizer)
    }
    c.opts <- paste0(mod.type, "Control")
    merOpts <- getOption(c.opts)
    if (!is.null(merOpts)) {
        nn <- names(merOpts)
        nn.ok <- .get.checkingOpts(names(merOpts))
        if (length(nn.ignored <- setdiff(nn, nn.ok)) > 0) {
            warning("some options in ", shQuote(sprintf("getOption('%s')", 
                c.opts)), " ignored : ", paste(nn.ignored, collapse = ", "))
        }
        for (arg in nn.ok) {
            if (do.call(missing, list(arg))) 
                assign(arg, merOpts[[arg]])
        }
    }
    check.rankX <- match.arg(check.rankX)
    check.scaleX <- match.arg(check.scaleX)
    me <- sys.function()
    chk.cconv(check.conv.grad, me)
    chk.cconv(check.conv.singular, me)
    chk.cconv(check.conv.hess, me)
    if (mod.type == "glmer" && use.last.params && calc.derivs) {
        warning("using ", shQuote("use.last.params"), "=TRUE and ", 
            shQuote("calc.derivs"), "=TRUE with ", shQuote("glmer"), 
            " will not give backward-compatible results")
    }
    ret <- namedList(optimizer, restart_edge, boundary.tol, calc.derivs, 
        use.last.params, checkControl = namedList(check.nobs.vs.rankZ, 
            check.nobs.vs.nlev, check.nlev.gtreq.5, check.nlev.gtr.1, 
            check.nobs.vs.nRE, check.rankX, check.scaleX, check.formula.LHS), 
        checkConv = namedList(check.conv.grad, check.conv.singular, 
            check.conv.hess), optCtrl = optCtrl)
    if (mod.type == "glmer") {
        ret <- c(ret, namedList(tolPwrss, compDev, nAGQ0initStep))
        ret$checkControl <- c(ret$checkControl, namedList(check.response.not.const))
    }
    class(ret) <- c(c.opts, "merControl")
    ret
}

checkArgs <- checkArgs %ORifNotInLme4% function (type, ...) 
{
    l... <- list(...)
    if (isTRUE(l...[["sparseX"]])) 
        warning("sparseX = TRUE has no effect at present", call. = FALSE)
    if (length(l... <- list(...))) {
        if (!is.null(l...[["family"]])) {
            warning("calling lmer with family() is deprecated: please use glmer() instead", 
                call. = FALSE)
            type <- "glmer"
        }
        if (!is.null(l...[["method"]])) {
            msg <- paste("Argument", sQuote("method"), "is deprecated.")
            if (type == "lmer") 
                msg <- paste(msg, "Use the REML argument to specify ML or REML estimation.")
            else if (type == "glmer") 
                msg <- paste(msg, "Use the nAGQ argument to specify Laplace (nAGQ=1) or adaptive", 
                  "Gauss-Hermite quadrature (nAGQ>1).  PQL is no longer available.")
            warning(msg, call. = FALSE)
            l... <- l...[names(l...) != "method"]
        }
        if (length(l...)) {
            warning("extra argument(s) ", paste(sQuote(names(l...)), 
                collapse = ", "), " disregarded", call. = FALSE)
        }
    }
}

checkCtrlLevels <- checkCtrlLevels %ORifNotInLme4% function (cstr, val, smallOK = FALSE) 
{
    bvals <- c("message", "warning", "stop", "ignore")
    if (smallOK) 
        bvals <- outer(bvals, c("", "Small"), paste0)
    if (!is.null(val) && !val %in% bvals) 
        stop("invalid control level ", sQuote(val), " in ", cstr, 
            ": valid options are {", paste(sapply(bvals, sQuote), 
                collapse = ","), "}")
    invisible(NULL)
}

checkFormulaData <- checkFormulaData %ORifNotInLme4% function (formula, data, checkLHS = TRUE, checkData = TRUE, debug = FALSE) 
{
    nonexist.data <- missDataFun(data)
    wd <- tryCatch(eval(data), error = identity)
    if (wrong.data <- inherits(wd, "simpleError")) {
        wrong.data.msg <- wd$message
    }
    bad.data <- nonexist.data || wrong.data
    if (bad.data || debug) {
        varex <- function(v, env) exists(v, envir = env, inherits = FALSE)
        allvars <- all.vars(as.formula(formula))
        allvarex <- function(env, vvec = allvars) all(vapply(vvec, 
            varex, NA, env))
    }
    if (bad.data) {
        if (allvarex(environment(formula))) {
            stop("bad 'data', but variables found in environment of formula: ", 
                "try specifying 'formula' as a formula rather ", 
                "than a string in the original model", call. = FALSE)
        }
        else {
            if (nonexist.data) {
                stop("'data' not found, and some variables missing from formula environment", 
                  call. = FALSE)
            }
            else {
                stop("bad 'data': ", wrong.data.msg)
            }
        }
    }
    else {
        denv <- if (is.null(data)) {
            if (!is.null(ee <- environment(formula))) {
                ee
            }
            else {
                parent.frame(2L)
            }
        }
        else list2env(data)
    }
    if (debug) {
        cat("Debugging parent frames in checkFormulaData:\n")
        glEnv <- 1L
        while (!identical(parent.frame(glEnv), .GlobalEnv)) {
            glEnv <- glEnv + 1L
        }
        for (i in 1:glEnv) {
            OK <- allvarex(parent.frame(i))
            cat("vars exist in parent frame ", i)
            if (i == glEnv) 
                cat(" (global)")
            cat(" ", OK, "\n")
        }
        cat("vars exist in env of formula ", allvarex(denv), 
            "\n")
    }
    stopifnot(!checkLHS || length(as.formula(formula, env = denv)) == 
        3)
    return(denv)
}

mkReTrms <- mkReTrms %ORifNotInLme4% function (bars, fr, drop.unused.levels = TRUE, reorder.terms = TRUE, 
    reorder.vars = FALSE) 
{
    if (!length(bars)) 
        stop("No random effects terms specified in formula", 
            call. = FALSE)
    stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr, 
        "data.frame"))
    names(bars) <- barnames(bars)
    term.names <- vapply(bars, safeDeparse, "")
    blist <- lapply(bars, mkBlist, fr, drop.unused.levels, reorder.vars = reorder.vars)
    nl <- vapply(blist, `[[`, 0L, "nl")
    if (reorder.terms) {
        if (any(diff(nl) > 0)) {
            ord <- rev(order(nl))
            blist <- blist[ord]
            nl <- nl[ord]
            term.names <- term.names[ord]
        }
    }
    Ztlist <- lapply(blist, `[[`, "sm")
    Zt <- do.call(rbind, Ztlist)
    names(Ztlist) <- term.names
    q <- nrow(Zt)
    cnms <- lapply(blist, `[[`, "cnms")
    nc <- lengths(cnms)
    nth <- as.integer((nc * (nc + 1))/2)
    nb <- nc * nl
    if (sum(nb) != q) {
        stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)", 
            sum(nb), q))
    }
    boff <- cumsum(c(0L, nb))
    thoff <- cumsum(c(0L, nth))
    Lambdat <- Matrix::t(do.call(Matrix::sparseMatrix, do.call(rbind, lapply(seq_along(blist), 
        function(i) {
            mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
            dd <- diag(nc[i])
            ltri <- lower.tri(dd, diag = TRUE)
            ii <- row(dd)[ltri]
            jj <- col(dd)[ltri]
            data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[, 
                jj]) + boff[i], x = as.double(rep.int(seq_along(ii), 
                rep.int(nl[i], length(ii))) + thoff[i]))
        }))))
    thet <- numeric(sum(nth))
    ll <- list(Zt = Matrix::drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x), 
               Gp = unname(c(0L, cumsum(nb))))
    ll$lower <- -Inf * (thet + 1)
    ll$lower[unique(Matrix::diag(Lambdat))] <- 0
    ll$theta[] <- is.finite(ll$lower)
    Lambdat@x[] <- ll$theta[ll$Lind]
    ll$Lambdat <- Lambdat
    fl <- lapply(blist, `[[`, "ff")
    fnms <- names(fl)
    if (length(fnms) > length(ufn <- unique(fnms))) {
        fl <- fl[match(ufn, fnms)]
        asgn <- match(fnms, ufn)
    }
    else asgn <- seq_along(fl)
    names(fl) <- ufn
    attr(fl, "assign") <- asgn
    ll$flist <- fl
    ll$cnms <- cnms
    ll$Ztlist <- Ztlist
    ll
}

findbars <- findbars %ORifNotInLme4% function (term) 
{
    fb <- function(term) {
        if (is.name(term) || !is.language(term)) 
            return(NULL)
        if (term[[1]] == as.name("(")) 
            return(fb(term[[2]]))
        stopifnot(is.call(term))
        if (term[[1]] == as.name("|")) 
            return(term)
        if (length(term) == 2) 
            return(fb(term[[2]]))
        c(fb(term[[2]]), fb(term[[3]]))
    }
    expandSlash <- function(bb) {
        makeInteraction <- function(x) {
            if (length(x) < 2) 
                return(x)
            trm1 <- makeInteraction(x[[1]])
            trm11 <- if (is.list(trm1)) 
                trm1[[1]]
            else trm1
            list(substitute(foo:bar, list(foo = x[[2]], bar = trm11)), 
                trm1)
        }
        slashTerms <- function(x) {
            if (!("/" %in% all.names(x))) 
                return(x)
            if (x[[1]] != as.name("/")) 
                stop("unparseable formula for grouping factor", 
                  call. = FALSE)
            list(slashTerms(x[[2]]), slashTerms(x[[3]]))
        }
        if (!is.list(bb)) 
            expandSlash(list(bb))
        else unlist(lapply(bb, function(x) {
            if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
                lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                  bar, list(foo = x[[2]], bar = trm)))
            else x
        }))
    }
    modterm <- expandDoubleVerts(if (inherits(term, "formula")) 
        term[[length(term)]]
    else term)
    expandSlash(fb(modterm))
}

RHSForm <- RHSForm %ORifNotInLme4% function (form, as.form = FALSE) 
{
    rhsf <- form[[length(form)]]
    if (as.form) 
        reformulate(deparse(rhsf))
    else rhsf
}

factorize <- factorize %ORifNotInLme4% function (x, frloc, char.only = FALSE) 
{
    for (i in all.vars(RHSForm(x))) {
        if (!is.null(curf <- frloc[[i]])) 
            frloc[[i]] <- makeFac(curf, char.only)
    }
    return(frloc)
}

checkNlevels <- checkNLevels %ORifNotInLme4% function (flist, n, ctrl, allow.n = FALSE) 
{
    stopifnot(is.list(ctrl), is.numeric(n))
    nlevelVec <- vapply(flist, function(x) nlevels(factor(x, 
        exclude = NA)), 1)
    cstr <- "check.nlev.gtr.1"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 2)) {
        wstr <- "grouping factors must have > 1 sampled level"
        switch(cc, warning = warning(wstr, call. = FALSE), stop = stop(wstr, 
            call. = FALSE), stop(gettextf("unknown check level for '%s'", 
            cstr), domain = NA))
    }
    else wstr <- character()
    cstr <- "check.nobs.vs.nlev"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(if (allow.n) nlevelVec > n else nlevelVec >= 
        n)) {
        wst2 <- gettextf("number of levels of each grouping factor must be %s number of observations", 
            if (allow.n) 
                "<="
            else "<")
        switch(cc, warning = warning(wst2, call. = FALSE), stop = stop(wst2, 
            call. = FALSE))
    }
    else wst2 <- character()
    cstr <- "check.nlev.gtreq.5"
    checkCtrlLevels(cstr, cc <- ctrl[[cstr]])
    if (doCheck(cc) && any(nlevelVec < 5)) {
        wst3 <- "grouping factors with < 5 sampled levels may give unreliable estimates"
        switch(cc, warning = warning(wst3, call. = FALSE), stop = stop(wst3, 
            call. = FALSE), stop(gettextf("unknown check level for '%s'", 
            cstr), domain = NA))
    }
    else wst3 <- character()
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
                "number of random effects", sprintf(" for term (%s)", 
                  term.names[i]))
            if (ww$unident) {
                switch(cc, warning = warning(ww$wstr, call. = FALSE), 
                  stop = stop(ww$wstr, call. = FALSE), stop(gettextf("unknown check level for '%s'", 
                    cstr), domain = NA))
                ww$wstr
            }
            else character()
        })))
    }
    else character()
}

checkZrank <- checkZrank %ORifNotInLme4% function (Zt, n, ctrl, nonSmall = 1e+06, allow.n = FALSE) 
{
    stopifnot(is.list(ctrl), is.numeric(n), is.numeric(nonSmall))
    cstr <- "check.nobs.vs.rankZ"
    if (doCheck(cc <- ctrl[[cstr]])) {
        checkCtrlLevels(cstr, cc, smallOK = TRUE)
        d <- dim(Zt)
        doTr <- d[1L] < d[2L]
        if (!(grepl("Small", cc) && prod(d) > nonSmall)) {
            rankZ <- Matrix::rankMatrix(if (doTr) 
                t(Zt)
            else Zt, method = "qr", sval = numeric(min(d)))
            ww <- wmsg(n, rankZ, allow.n, "number of observations", 
                "rank(Z)")
            if (ww$unident) {
                switch(cc, warningSmall = , warning = warning(ww$wstr, 
                  call. = FALSE), stopSmall = , stop = stop(ww$wstr, 
                  call. = FALSE), stop(gettextf("unknown check level for '%s'", 
                  cstr), domain = NA))
                ww$wstr
            }
            else character()
        }
        else character()
    }
    else character()
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

subbars <- subbars %ORifNotInLme4% function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
    if (is.call(term) && term[[1]] == as.name("||")) 
        term[[1]] <- as.name("+")
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

chkRank.drop.cols <- chkRank.drop.cols %ORifNotInLme4% function (X, kind, tol = 1e-07, method = "qr.R") 
{
    stopifnot(is.matrix(X))
    kinds <- eval(formals(lmerControl)[["check.rankX"]])
    if (!kind %in% kinds) 
        stop(gettextf("undefined option for 'kind': %s", kind))
    if (kind == "ignore") 
        return(X)
    p <- ncol(X)
    if (kind == "stop.deficient") {
        if ((rX <- Matrix::rankMatrix(X, tol = tol, method = method)) < 
            p) 
            stop(gettextf(sub("\n +", "\n", "the fixed-effects model matrix is column rank deficient (rank(X) = %d < %d = p);\n                   the fixed effects will be jointly unidentifiable"), 
                rX, p), call. = FALSE)
    }
    else {
        qr.X <- qr(X, tol = tol, LAPACK = FALSE)
        rnkX <- qr.X$rank
        if (rnkX == p) 
            return(X)
        msg <- sprintf(ngettext(p - rnkX, "fixed-effect model matrix is rank deficient so dropping %d column / coefficient", 
            "fixed-effect model matrix is rank deficient so dropping %d columns / coefficients"), 
            p - rnkX)
        if (kind != "silent.drop.cols") 
            (if (kind == "warn+drop.cols") 
                warning
            else message)(msg, domain = NA)
        contr <- attr(X, "contrasts")
        asgn <- attr(X, "assign")
        keep <- qr.X$pivot[seq_len(rnkX)]
        dropped.names <- colnames(X[, -keep, drop = FALSE])
        X <- X[, keep, drop = FALSE]
        if (Matrix::rankMatrix(X, tol = tol, method = method) < ncol(X)) 
            stop(gettextf("Dropping columns failed to produce full column rank design matrix"), 
                call. = FALSE)
        if (!is.null(contr)) 
            attr(X, "contrasts") <- contr
        if (!is.null(asgn)) 
            attr(X, "assign") <- asgn[keep]
        attr(X, "msgRankdrop") <- msg
        attr(X, "col.dropped") <- setNames(qr.X$pivot[(rnkX + 
            1L):p], dropped.names)
    }
    X
}

checkScaleX <- checkScaleX %ORifNotInLme4% function (X, kind = "warning", tol = 1000) 
{
    kinds <- eval(formals(lmerControl)[["check.scaleX"]])
    if (!kind %in% kinds) 
        stop(gettextf("unknown check-scale option: %s", kind))
    if (is.null(kind) || kind == "ignore") 
        return(X)
    cont.cols <- apply(X, 2, function(z) !all(z %in% c(0, 1)))
    col.sd <- apply(X[, cont.cols, drop = FALSE], 2L, sd)
    sdcomp <- outer(col.sd, col.sd, "/")
    logcomp <- abs(log(sdcomp[lower.tri(sdcomp)]))
    logsd <- abs(log(col.sd))
    if (any(c(logcomp, logsd) > log(tol))) {
        wmsg <- "Some predictor variables are on very different scales:"
        if (kind %in% c("warning", "stop")) {
            wmsg <- paste(wmsg, "consider rescaling")
            switch(kind, warning = warning(wmsg, call. = FALSE), 
                stop = stop(wmsg, call. = FALSE))
        }
        else {
            wmsg <- paste(wmsg, "auto-rescaled (results NOT adjusted)")
            X[, cont.cols] <- sweep(X[, cont.cols, drop = FALSE], 
                2, col.sd, "/")
            attr(X, "scaled:scale") <- setNames(col.sd, colnames(X)[cont.cols])
            if (kind == "warn+rescale") 
                warning(wmsg, call. = FALSE)
        }
    }
    else wmsg <- character()
    structure(X, msgScaleX = wmsg)
}

missDataFun <- missDataFun %ORifNotInLme4% function(d) {
    ff <- sys.frames()
    ex <- substitute(d)
    ii <- rev(seq_along(ff))
    foundAnon <- FALSE
    for(i in ii) {
        ex <- eval(substitute(substitute(x, env=sys.frames()[[n]]),
                              env = list(x = ex, n=i)))
        if (is.symbol(ex) && grepl("^\\.\\.[0-9]+$",safeDeparse(ex))) {
            foundAnon <- TRUE
            break
        }
    }
    return(!foundAnon && is.symbol(ex) && !exists(deparse(ex)))
}

safeDeparse <- safeDeparse %ORifNotInLme4% function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)


makeFac <- makeFac %ORifNotInLme4% function (x, char.only = FALSE) 
{
    if (!is.factor(x) && (!char.only || is.character(x))) 
        factor(x)
    else x
}

expandDoubleVerts <- expandDoubleVerts %ORifNotInLme4% function (term) 
{
    expandDoubleVert <- function(term) {
        frml <- formula(substitute(~x, list(x = term[[2]])))
        newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
        if (attr(terms(frml), "intercept") != 0) 
            newtrms <- c("1", newtrms)
        as.formula(paste("~(", paste(vapply(newtrms, function(trm) paste0(trm, 
            "|", deparse(term[[3]])), ""), collapse = ")+("), 
            ")"))[[2]]
    }
    if (!is.name(term) && is.language(term)) {
        if (term[[1]] == as.name("(")) {
            term[[2]] <- expandDoubleVerts(term[[2]])
        }
        stopifnot(is.call(term))
        if (term[[1]] == as.name("||")) 
            return(expandDoubleVert(term))
        term[[2]] <- expandDoubleVerts(term[[2]])
        if (length(term) != 2) {
            if (length(term) == 3) 
                term[[3]] <- expandDoubleVerts(term[[3]])
        }
    }
    term
}

barnames <- barnames %ORifNotInLme4% function (bars) 
  vapply(bars, function(x) safeDeparse(x[[3]]), "")

mkBlist <- mkBlist %ORifNotInLme4% function (x, frloc, drop.unused.levels = TRUE, reorder.vars = FALSE) 
{
    frloc <- factorize(x, frloc)
    if (is.null(ff <- tryCatch(eval(substitute(makeFac(fac), 
        list(fac = x[[3]])), frloc), error = function(e) NULL))) 
        stop("couldn't evaluate grouping factor ", deparse(x[[3]]), 
            " within model frame:", " try adding grouping factor to data ", 
            "frame explicitly if possible", call. = FALSE)
    if (all(is.na(ff))) 
        stop("Invalid grouping factor specification, ", deparse(x[[3]]), 
            call. = FALSE)
    if (drop.unused.levels) 
        ff <- factor(ff, exclude = NA)
    nl <- length(levels(ff))
    mm <- model.matrix(eval(substitute(~foo, list(foo = x[[2]]))), 
        frloc)
    if (reorder.vars) {
        mm <- mm[colSort(colnames(mm)), ]
    }
    sm <- Matrix::fac2sparse(ff, to = "d", drop.unused.levels = drop.unused.levels)
    sm <- Matrix::KhatriRao(sm, t(mm))
    dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))
    list(ff = ff, sm = sm, nl = nl, cnms = colnames(mm))
}

doCheck <- doCheck %ORifNotInLme4% function (x) 
{
    is.character(x) && !any(x == "ignore")
}

nobars_ <- nobars_ %ORifNotInLme4% function (term) 
{
    if (!anyBars(term)) 
        return(term)
    if (isBar(term)) 
        return(NULL)
    if (isAnyArgBar(term)) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars_(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars_(term[[2]])
    nb3 <- nobars_(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}


nobars <- nobars %ORifNotInLme4% function (term) 
{
    nb <- nobars_(term)
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

`RHSForm<-` <- `RHSForm<-` %ORifNotInLme4% function (formula, value) 
{
    formula[[length(formula)]] <- value
    formula
}

reOnly <- function (f, response = FALSE) 
{
  found_bars <- findbars(f)
  if (is.null(found_bars)) return(~0)
  
  reformulate(paste0("(", vapply(findbars(f), safeDeparse, 
      ""), ")"), response = if (response && length(f) == 3L) 
      f[[2]])
}

mkVarCorr <- mkVarCorr %ORifNotInLme4% function (sc, cnms, nc, theta, nms) 
{
    ncseq <- seq_along(nc)
    thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
    if (!all(nms == names(cnms))) 
        warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n", 
            "Please report!", immediate. = TRUE)
    ans <- lapply(ncseq, function(i) {
        Li <- diag(nrow = nc[i])
        Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
        rownames(Li) <- cnms[[i]]
        val <- tcrossprod(sc * Li)
        stddev <- sqrt(diag(val))
        corr <- t(val/stddev)/stddev
        diag(corr) <- 1
        structure(val, stddev = stddev, correlation = corr)
    })
    if (is.character(nms)) {
        if (anyDuplicated(nms)) 
            nms <- make.names(nms, unique = TRUE)
        names(ans) <- nms
    }
    structure(ans, sc = sc)
}

wmsg <- wgmsg %ORifNotInLme4% function (n, cmp.val, allow.n, msg1 = "", msg2 = "", msg3 = "") 
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

isSingular <- isSingular %ORifNotInLme4% function (x, tol = 1e-04) 
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

isAnyArgBar <- isAnyArgBar %ORifNotInLme4% function (term) 
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

# can't use lme4 version, as we need to use mstan4bartFit terms and model.frame
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

formula.mstan4bartFit <- function(x, type = c("all", "fixed", "random", "bart"))
{
  if (is.null(formula <- x$formula)) {
    if (!grepl("mstan4bart$", deparse(x$call[[1]]))) 
      stop("can't find formula stored in model frame or call")
    form <- as.formula(formula(x$call))
  }
  
  type <- match.arg(type)
  
  if (type == "fixed") {
    RHSForm(formula) <- nobart(nobars(RHSForm(formula)))
  } else if (type == "random") {
    formula <- reOnly(formula, response = TRUE)
  } else if (type == "bart") {
    RHSForm(formula) <- allbart(nobars(RHSForm(formula)))
  }
  
  formula
}

terms.mstan4bartFit <- function(x, type = c("all", "fixed", "random", "bart")) 
{
  type <- match.arg(type)
  
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

model.frame.mstan4bartFit <- function (formula, type = c("all", "fixed", "random", "bart"))
{
  type <- match.arg(type)
  
  frame <- formula$frame
  
  if (type == "fixed") {
    frame <- frame[attr(terms(formula), "varnames.fixed")]
  } else if (type == "random") {
    frame <- frame[attr(terms(formula), "varnames.random")]
  } else if (type == "bart") {
    frame <- frame[attr(terms(formula), "varnames.bart")]
  }
  frame
}

rm(`%ORifNotInLme4%`)
if (exists("lme4ns")) rm("lme4ns")

make_glmerControl <- function(..., ignore_lhs = FALSE, ignore_x_scale = FALSE) {
  lme4::glmerControl(check.nlev.gtreq.5 = "ignore",
               check.nlev.gtr.1 = "stop",
               check.nobs.vs.rankZ = "ignore",
               check.nobs.vs.nlev = "ignore",
               check.nobs.vs.nRE = "ignore", 
               check.formula.LHS = if (ignore_lhs) "ignore" else "stop",
               check.scaleX = if (ignore_x_scale) "ignore" else "warning",
               ...)  
}

