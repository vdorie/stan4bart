.onUnload <- function(libpath)
{
  ## gc is necessary to collect external pointers who have not yet been collected
  ## that have finalizers pointing to the soon-to-unloaded dll
  gc(FALSE)
  if (is.loaded("stan4bart_finalize", PACKAGE = "stan4bart")) {
    .Call(C_stan4bart_finalize)
    library.dynam.unload("stan4bart", libpath)
  }
}
