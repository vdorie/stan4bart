nlist <- function (...)  {
  m <- match.call()
  out <- list(...)
  no_names <- is.null(names(out))
  has_name <- if (no_names) 
    FALSE
  else nzchar(names(out))
  if (all(has_name)) 
    return(out)
  nms <- as.character(m)[-1L]
  if (no_names) {
    names(out) <- nms
  }
  else {
    names(out)[!has_name] <- nms[!has_name]
  }
  return(out)
}

"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

format.time <- function(secs) {
  hours <- floor(secs / 3600)
  secs <- secs - hours * 3600
  mins  <- floor(secs / 60)
  secs <- round(secs - mins * 60)
  
  if (hours > 0) paste0(hours, ":", mins, ":", secs) else paste0(mins, ":", secs)
}

