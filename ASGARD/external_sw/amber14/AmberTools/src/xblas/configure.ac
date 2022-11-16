dnl Copyright (c) 2005, Regents of the University of California
AC_INIT(xblas, 1.0, extended_blas@cs.berkeley.edu)
AC_PREREQ(2.59)
AC_CONFIG_SRCDIR(src/blas_extended.h)
AC_CONFIG_AUX_DIR(config)

AC_CANONICAL_TARGET

AC_PROG_CC

CC_BASE=`basename $CC`
case $target in
  powerpc*)
    case "$CC_BASE" in
      cc* | xlc* | c99 | c89 )
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS -qstrict"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [], [CFLAGS="$save_CFLAGS"])
      ;;
    esac
  ;;
  i?86-*-* | k?*-*-* | athlon-*-* | ia64-*-* | x86_64-*-*)
    case "$CC_BASE" in
      icc | ecc)
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS -mp"
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM], [], [CFLAGS="$save_CFLAGS"])
      ;;
    esac
  ;;
esac

AC_PROG_RANLIB

AC_ARG_ENABLE([fortran], AC_HELP_STRING(--enable-fortran,
    [build Fortran wrapper.  [[default=auto]]]))
if test "x$enable_fortran" = "x"; then
  enable_fortran="auto"
fi

if test "x$enable_fortran" != "xno"; then
  AC_PROG_FC([], [Fortran 77])
  if test "x$FC" = "x"; then
    enable_fortran="no"
  else
    enable_fortran="yes"
    # Set appropriate C preprocessor definitions based on name mangling scheme.
    AC_FC_FUNC(f_func)
    case "$f_func" in
      f_func_)
        CFLAGS="$CFLAGS -DCONFIG_FC_UNDERSCORE"
      ;;
      f_func__)
        CFLAGS="$CFLAGS -DCONFIG_FC_DBL_UNDERSCORE"
      ;;
      F_FUNC_)
        CFLAGS="$CFLAGS -DCONFIG_FC_UCASE -DCONFIG_FC_UNDERSCORE"
      ;;
      F_FUNC__)
        CFLAGS="$CFLAGS -DCONFIG_FC_UCASE -DCONFIG_FC_DBL_UNDERSCORE"
      ;;
      F_FUNC)
        CFLAGS="$CFLAGS -DCONFIG_FC_UCASE"
      ;;
      f_func) ;;
      *)
        AC_ERROR([Could not figure out Fortran name mangling scheme.])
      ;;
    esac

    # If --enable-xerbla is specified we need to link in Fortran libraries
    AC_ARG_ENABLE(xerbla, AC_HELP_STRING(--enable-xerbla, [use LAPACK's xerbla in BLAS_error.  [[default=no]]]))
    AC_ARG_VAR(EXTRA_LIBS, [Extra libraries used for linking Fortran libraries with C programs])
    if test "x$enable_xerbla" = "xyes"; then
      CFLAGS="$CFLAGS -DCONFIG_USE_XERBLA"
      if test "x$FC_DUMMY_MAIN" != "xnone"; then
        CFLAGS="$CFLAGS -DFC_DUMMY_MAIN=$FC_DUMMY_MAIN"
      fi
      EXTRA_LIBS="$EXTRA_LIBS $FCLIBS"
    fi
  fi
fi

x86_fix="no"
case "$target" in
  i?86-*-* | k?*-*-* | athlon-*-* | x86_64-*-*) 
    AC_CHECK_HEADER([fpu_control.h], [CFLAGS="$CFLAGS -DHAVE_FPU_CONTROL_H"])
  ;;
esac

AC_MSG_CHECKING([whether x86 FPU control code is needed])
x86_fix="no"
case "$target" in
  i?86-*-* | k?*-*-* | athlon-*-* | x86_64-*-*) 
    CFLAGS="$CFLAGS -Dx86"
    x86_fix="yes"
  ;;
esac
AC_MSG_RESULT($x86_fix)

# Check for m4
AC_ARG_VAR(M4, [M4 program])
AC_CHECK_PROG(M4, [m4], [m4])
if test "x$M4" = "x"; then
  AC_ERROR([M4 not found])
fi

# Non-GNU m4 require expanded buffer size
# Detect GNU m4 by its behaviour on __gnu__.
dnl Note we need to quote __gnu__ since autoconf goes through GNU m4.
AC_ARG_VAR(M4_OPTS, [Flags to pass to M4])
AC_MSG_CHECKING([whether we are using GNU m4])
tmp=`echo [__gnu__] | $M4 -`
if test "x$tmp" != "x"; then
  M4_OPTS="$M4_OPTS -B65536"
  AC_MSG_RESULT([no])
else
  AC_MSG_RESULT([yes])
fi

# Check for indent.  If none found, use : (do nothing).
AC_ARG_VAR(INDENT, [Indent program])
AC_ARG_VAR(INDENT_OPTS, [Flags to pass to indent])
AC_CHECK_PROG(INDENT, [indent], [indent], [:])
if test "x$INDENT_OPTS" == "x" && test "$INDENT" != ":"; then
  INDENT_OPTS="-ce -i2 -nfc1 -br -brs -cs -npcs -nprs -npsl"
  cat >conftest.c <<EOF
int main() { }
EOF
  $INDENT $INDENT_OPTS conftest.c || INDENT_OPTS=""
  rm -f conftest.c
fi

AC_ARG_ENABLE([plain-blas], 
  AC_HELP_STRING(--enable-plain-blas, 
    [build plain (no extended or mixed precision) BLAS.  [[default=yes]]]))
if test "x$enable_plain_blas" = "xno"; then
  M4_OPTS="$M4_OPTS -D no_plain_blas"
fi

if test "x$enable_fortran" = "xno"; then
  M4_OPTS="$M4_OPTS -D no_f2c"
fi

EXTRA_LIBS="$EXTRA_LIBS -lm"

AC_CONFIG_COMMANDS([gen_makefile], 
                   [$M4 $M4_OPTS Makefile.m4 >Makefile], [M4=$M4; M4_OPTS="$M4_OPTS"])

AC_SUBST(EXTRA_LIBS, $EXTRA_LIBS)
AC_CONFIG_FILES([make.inc])
AC_OUTPUT

