"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

# Every symbol appearing in function position anywhere in an expression.
called_names <- function(expr) {
  if (!is.call(expr)) return(character())
  c(if (is.symbol(expr[[1L]])) as.character(expr[[1L]]),
    unlist(lapply(as.list(expr), called_names)))
}

quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1L]] <- as.symbol(":::")
  result[[2L]] <- as.symbol("stan4bart")
  
  result[[3L]] <- if (character.only) name else match.call()[[2]]
  result
}

strip_extra_terms_from_language <- function(expr, extra_terms) {
  if (length(expr) == 3L) {
    # binary operator
    
    # Very annoying problem where
    #   any(as.character(a ~ `(offset)`) == "`(offset)`") => TRUE
    #   as.character((a ~ `(offset)`)[3L]) == "`(offset)`" => FALSE
    #
    # Consequently, we have to compare the whole expression even though
    # we aren't interested in the operator
    terms_are_extra <- as.character(expr) %in% extra_terms
    if (terms_are_extra[2L]) return(expr[[3L]])
    if (terms_are_extra[3L]) return(expr[[2L]])
    
    expr[[2L]] <- strip_extra_terms_from_language(expr[[2L]], extra_terms)
    expr[[3L]] <- strip_extra_terms_from_language(expr[[3L]], extra_terms)
  } else if (length(expr) == 2L) {
    # unary operator
    terms_are_extra <- as.character(expr) %in% extra_terms
    if (terms_are_extra[2L]) return(quote(0))
    
    expr[[2L]] <- strip_extra_terms_from_language(expr[[2L]], extra_terms)
  }
  
  expr
}

strip_extra_terms <- function(terms, extra_terms) {
  attr(terms, "term.labels") <- setdiff(attr(terms, "term.labels"), extra_terms)
  factors <- attr(terms, "factors")
  attr(terms, "factors") <-
    factors[setdiff(rownames(factors), extra_terms), setdiff(colnames(factors), extra_terms)]
  
  variables <- attr(terms, "variables")
  attr(terms, "variables") <-
    variables[as.character(variables) %not_in% extra_terms]
  
  terms <- strip_extra_terms_from_language(terms, extra_terms)
  
  terms
}

delete.weights <- function(termobj, weights)
{
  a <- attributes(termobj)
  termobj <- strip_extra_terms_from_language(termobj, weights)
  
  w <- which(sapply(a$variables, "==", weights)) - 1L
  if (length(w) == 1L && w > 0L) {
    a$variables <- a$variables[-(1 + w)]
    a$predvars  <- a$predvars[-(1 + w)]
    if (length(a$factors) > 0L)
      a$factors <- a$factors[-w,,drop = FALSE]
    if (length(a$offset) > 0L) 
      a$offset <- ifelse(a$offset > w, a$offset - 1, a$offset)
    if (length(a$specials) > 0L) {
      for (i in seq_along(a$specials)) {
        b <- a$specials[[i]]
        a$specials[[i]] <- ifelse(b > w, b - 1, b)
      }
    }  
    attributes(termobj) <- a
  }
  
  termobj
}

# Whether lme4 does the mixed-model fits: it is installed, and the option
# 'stan4bart.lme4' is not FALSE. The option is how a test reaches the paths
# that run without lme4 on a machine that has it.
lme4_is_used <- function()
{
  !isFALSE(getOption("stan4bart.lme4")) && requireNamespace("lme4", quietly = TRUE)
}

# An argument's expression as written, and the environment it was written in.
# match.call() records an argument forwarded through a wrapper's dots as ..N,
# and evaluating that forces the wrapper caller's promise outside whatever the
# caller wants bound. Each forwarding frame's own call names the Nth dots
# element, so the walk back reaches the original expression; it stops where
# that is not possible and leaves the reference to evaluate as an ordinary one.
# Must run while the forwarding frames are still on the stack.
recover_forwarded_argument <- function(expr, env) {
  while (is_dots_reference(expr)) {
    recovered <- tryCatch(
      {
        # a closure can reference dots its enclosing function owns
        while (!exists("...", envir = env, inherits = FALSE))
          env <- parent.env(env)
        # the frame's first place on the stack is its own call; later ones are
        # evaluations in it (eval(call, env)), whose caller is not the one that
        # wrote the dots
        frame <- Position(function(f) identical(f, env), sys.frames())
        parent <- sys.parents()[frame]
        # NextMethod(name = value) replaces the method's dots, but the frame
        # still records the generic's call
        if (frame > 1L && identical(sys.function(frame - 1L), NextMethod))
          stop("dots replaced by NextMethod")
        # a caller that is no frame on the stack, do.call(envir = ), numbers as
        # the frame itself and cannot be named
        if (parent >= frame) stop("unknown caller")
        caller <- sys.frame(parent)
        dots <- match.call(sys.function(frame), sys.call(frame), expand.dots = FALSE,
                           envir = referenced_dots(caller))$...
        list(dots[[as.integer(substring(expr, 3L))]], caller)
      },
      error = function(e) NULL)
    if (is.null(recovered)) break
    expr <- recovered[[1L]]
    env <- recovered[[2L]]
  }
  list(expr = expr, env = env)
}

# An environment whose dots are those of `env` as references, ..1 to ..N under
# their own names, so that match.call() records an argument passed on through
# them as a reference even where it is a constant; `env` itself where it sees
# no dots.
referenced_dots <- function(env) {
  if (!exists("...", envir = env)) return(env)
  dots <- as.list(match.call(function(...) NULL, quote(f(...)), envir = env))
  references <- lapply(seq_along(dots[-1L]), function(i) as.name(paste0("..", i)))
  names(references) <- names(dots[-1L])
  do.call(function(...) environment(), references)
}

is_dots_reference <- function(expr)
  is.symbol(expr) && grepl("^\\.\\.[1-9][0-9]*$", expr)
