#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(stan4bart, 0.0.4, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/init.cpp])
AC_CONFIG_AUX_DIR([tools/build-aux])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
