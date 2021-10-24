#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(stan4bart, 0.0.1, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/init.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

TBB_LIB=`"${R_HOME}/bin/Rscript" -e "
cat(system.file('lib', .Platform[['r_arch']], package = 'RcppParallel', mustWork = TRUE))" | grep "RcppParallel"`

PKG_CXXFLAGS=`"${R_HOME}/bin/Rscript" -e "RcppParallel::CxxFlags()" -e "
TBB <- system.file('include', package = 'RcppParallel', mustWork = TRUE)
cat(paste0('-I', shQuote(TBB), ' -D_REENTRANT -DSTAN_THREADS'), ' ')"`

RCPP_PKG_LIBS=`"${R_HOME}/bin/Rscript" -e "RcppParallel::RcppParallelLibs()"`
TBB_PKG_LIBS="-L${TBB_LIB} -Wl,-rpath,'${TBB_LIB}' -ltbb -ltbbmalloc"
PKG_LIBS="$RCPP_PKG_LIBS $TBB_PKG_LIBS"

AC_SUBST(TBB_LIB)
AC_SUBST(PKG_CXXFLAGS)
AC_SUBST(PKG_LIBS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

