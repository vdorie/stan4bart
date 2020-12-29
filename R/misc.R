"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

quoteInNamespace <- function(name, character.only = FALSE) {
  result <- quote(a + b)
  result[[1L]] <- as.symbol(":::")
  result[[2L]] <- as.symbol("stan4bart")
  
  result[[3L]] <- if (character.only) name else match.call()[[2]]
  result
}

