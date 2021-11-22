#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(stan4bart, 0.0.2, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/init.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])
AC_CONFIG_HEADERS([src/config.h])

m4_include([tools/m4/ax_check_compile_flag.m4])
m4_include([tools/m4/ax_gcc_x86_avx_xgetbv.m4])
m4_include([tools/m4/ax_gcc_x86_cpuid.m4])
m4_include([tools/m4/ax_ext.m4])
m4_include([tools/m4/ax_func_posix_memalign.m4])

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

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX


AC_CANONICAL_HOST

if test "${host_cpu}-${host_vendor}" == "x86_64-apple"; then
  AC_PROG_CC
  AC_LANG_PUSH([C])
  AC_MSG_CHECKING([if compiling under rosetta2])
  AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <stddef.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <errno.h>
    ], [
// returns 0 if translated, 1 if verifiably not, -1 otherwise
int main(void) {
  int ret = 0;
  size_t size = sizeof(ret);
  if (sysctlbyname("sysctl.proc_translated", &ret, &size, NULL, 0) != -1) 
	  return ret == 1 ? 0 : 1;
  if (errno == ENOENT)
		return 1;
  return -1;
}
    ])], [under_rosetta2="yes"],
         [under_rosetta2="no"])
  
  AC_MSG_RESULT([$under_rosetta2])
  
  if test x"$under_rosetta2" = x"no"; then
    AX_EXT
  fi
else
  AX_EXT
fi

AC_CHECK_HEADERS([malloc.h])

AC_CHECK_SIZEOF([size_t])
AC_CHECK_ALIGNOF([void*])

AC_FUNC_ALLOCA
AX_FUNC_POSIX_MEMALIGN
AC_CHECK_FUNCS([ffs])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

