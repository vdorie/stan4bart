#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(mbarstan, 0.0-1, vdorie@gmail.com)
AC_CONFIG_SRCDIR([src/mbartstan.cpp])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

STANHEADERS_SRC=`${R_HOME}/bin${R_ARCH_BIN}/Rscript -e "message()" \
  -e "cat(system.file('include', 'src', package = 'StanHeaders', mustWork = TRUE))" \
  -e "message()" | grep "StanHeaders"`

AC_SUBST(STANHEADERS_SRC)

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT

