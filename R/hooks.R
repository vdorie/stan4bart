.onUnload <- function(libpath)
{
  ## gc is necessary to collect external pointers who have not yet been collected
  ## that have finalizers pointing to the soon-to-unloaded dll
  gc(FALSE)
  if (is.loaded("stan4bart_finalize", PACKAGE = "stan4bart")) {
    .Call("stan4bart_finalize", PACKAGE = "stan4bart")
    library.dynam.unload("stan4bart", libpath)
  }
}

.onLoad <- function(libname, pkgname) {
  if (.Platform$OS.type %in% "windows" && !("RcppParallel" %in% sapply(.dynLibs(), "[[", "name")))
    library.dynam("stan4bart", pkgname, libname,
                  DLLpath = file.path(system.file("lib", package = "RcppParallel"), .Platform$r_arch))
  else 
    library.dynam("stan4bart", pkgname, libname)
}


