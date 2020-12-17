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

STANHEADERS_SRC=`"${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "message()" -e "cat(system.file('include', 'src', package = 'StanHeaders', mustWork = TRUE))" -e "message()" | grep "StanHeaders"`

TBB=`"${R_HOME}/bin${R_ARCH_BIN}/Rscript" -e "message()" -e "cat(system.file('lib', package = 'RcppParallel', mustWork = TRUE))" -e "message()" | grep "RcppParallel"`

CPPFLAGS_1=`"${R_HOME}/bin/Rscript" -e "RcppParallel::CxxFlags()"`
CPPFLAGS_2=`"${R_HOME}/bin/Rscript" -e "StanHeaders:::CxxFlags()"`
EXTRA_CPPFLAGS="${CPPFLAGS_1} ${CPPFLAGS_2}"

AC_SUBST(STANHEADERS_SRC)
AC_SUBST(TBB)
AC_SUBST(EXTRA_CPPFLAGS)

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

