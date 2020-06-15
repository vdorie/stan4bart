glFormula <- function (formula, data = NULL, subset, weights, 
    na.action, offset, contrasts = NULL, treatment = NULL, start, mustart, etastart, 
    control = glmerControl(), ...) 
{
    control <- control$checkControl
    mf <- mc <- match.call()
   
    if (!is.null(mc$treatment)) {
        if (is.symbol(mc$treatment)) treatment <- as.character(mc$treatment)
      mf$treatment <- NULL
      if (!is.character(treatment))
        stop("treament must be a character or symbol")
    }
    
    ignoreArgs <- c("start", "verbose", "devFunOnly", "optimizer", 
        "control", "nAGQ")
    l... <- list(...)
    l... <- l...[!names(l...) %in% ignoreArgs]
    do.call(checkArgs, c(list("glmer"), l...))
    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr, control[[cstr]])
    denv <- checkFormulaData(formula, data, checkLHS = control$check.formula.LHS == 
        "stop")
    mc$formula <- formula <- as.formula(formula, env = denv)
    m <- match(c("data", "subset", "weights", "na.action", "offset", 
        "mustart", "etastart"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    fr.form <- subbars(formula)
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
    n <- nrow(fr)
    
    reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)
    
    if (!is.null(treatment)) {
      if (!(treatment %in% colnames(fr)))
        stop("treament must be the name of a column in data")
      uq <- sort(unique(fr[[treatment]]))
      if (length(uq) != 2L || !all(uq == c(0, 1)))
        stop("treatment must in { 0, 1 }")
      fr.cf <- fr
      fr.cf[[treatment]] <- 1 - fr.cf[[treatment]]
      reTrms.cf <- mkReTrms(findbars(RHSForm(formula)), fr.cf)
      if (!all(sapply(seq_along(reTrms$cnms), function(i) all(reTrms$cnms[[i]] == reTrms.cf$cnms[[i]]))))
        stop("counterfactual random effect design matrix does not match observed: contact package author")
    }
    wmsgNlev <- checkNlevels(reTrms$flist, n = n, control, allow.n = TRUE)
    wmsgZdims <- checkZdims(reTrms$Ztlist, n = n, control, allow.n = TRUE)
    wmsgZrank <- checkZrank(reTrms$Zt, n = n, control, nonSmall = 1e+06, 
        allow.n = TRUE)
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))
    mf$formula <- fixedform
    fixedfr <- eval(mf, parent.frame())
    attr(attr(fr, "terms"), "predvars.fixed") <- attr(attr(fixedfr, 
        "terms"), "predvars")
    #ranform <- formula
    #RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))
    #mf$formula <- ranform
    #ranfr <- eval(mf, parent.frame())
    #attr(attr(fr, "terms"), "predvars.random") <- attr(terms(ranfr), 
    #    "predvars")
    fr.bart <- fixedfr
    if ("(offset)" %in% colnames(fr.bart)) {
      keepcols <- colnames(fr.bart) != "(offset)"
      fr.bart <- fr.bart[,keepcols,drop = FALSE]
      attr(attr(fr.bart, "terms"), "dataClasses") <- attr(attr(fr.bart, "terms"), "dataClasses")[keepcols]
    }
    bartData <- dbarts::dbartsData(fixedform, fr.bart)
    if (!is.null(treatment)) {
      bartData@x.test <- bartData@x
      bartData@x.test[, treatment] <- 1 - bartData@x.test[, treatment, drop = FALSE]
      bartData@testUsesRegularOffset <- FALSE
    }
    #X <- model.matrix(fixedform, fr, contrasts)
    #if (is.null(rankX.chk <- control[["check.rankX"]])) 
    #   rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    #X <- chkRank.drop.cols(X, kind = rankX.chk, tol = 1e-07)
    #if (is.null(scaleX.chk <- control[["check.scaleX"]])) 
    #    scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    #X <- checkScaleX(X, kind = scaleX.chk)
    
    y <- model.response(fr)
    u.y <- unique(y)
    family <- if (length(u.y) == 2L && all(sort(u.y) == c(0, 1))) binomial(link = "probit") else gaussian()
    result <- list(fr = fr, bartData = bartData, reTrms = reTrms, family = family, formula = formula, 
                   wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
    if (!is.null(treatment))
      result$reTrms.cf <- reTrms.cf
    result
}

tryResult <- tryCatch(lme4ns <- asNamespace("lme4"), error = function(e) e)
if (!is(tryResult, "error")) {
  `%ORifNotInLme4%` <- function(a, b) {
    mc <- match.call()
    if (is.symbol(mc[["a"]])) a <- deparse(mc[["a"]])
    if (!is.null(lme4ns[[a]])) lme4ns[[a]] else b
  }
} else {
  `%ORifNotInLme4%` <- function(a, b) b
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
    modterm <- expandDoubleVerts(if (is(term, "formula")) 
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
    if (is(term, "formula") && length(term) == 3 && is.symbol(nb)) {
        nb <- reformulate("1", response = deparse(nb))
    }
    if (is.null(nb)) {
        nb <- if (is(term, "formula")) 
            ~1
        else 1
    }
    nb
}

`RHSForm<-` <- `RHSForm<-` %ORifNotInLme4% function (formula, value) 
{
    formula[[length(formula)]] <- value
    formula
}

reOnly <- reOnly %ORifNotInLme4% function (f, response = FALSE) 
{
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

rm(`%ORifNotInLme4%`)
if (exists("lme4ns")) rm("lme4ns")

