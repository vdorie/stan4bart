#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.64)
AC_INIT(stan4bart, 0.0.4, vdorie@gmail.com)
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

AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

AX_EXT

AC_CHECK_HEADERS([malloc.h])

AC_CHECK_SIZEOF([size_t])
AC_CHECK_ALIGNOF([void*])

AC_FUNC_ALLOCA
AX_FUNC_POSIX_MEMALIGN
AC_CHECK_FUNCS([ffs])

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

