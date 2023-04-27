"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

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
