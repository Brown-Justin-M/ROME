AC_INIT([ROME], [1.0], [jmbrown2@nps.edu])
AM_INIT_AUTOMAKE
AC_LANG([Fortran])


my_save_fcflags=$FCFLAGS
FCFLAGS="$FCFLAGS -cpp"
AS_IF([test "x$FC" != "x"],[
    find_compiler=no
])

AC_PROG_FC
FCFLAGS=$my_save_fcflags

AS_IF([test "x$FFTW_HOME" != "x"], [
    AC_MSG_NOTICE([FFTW_HOME set, updating paths])
    FCFLAGS="-I$FFTW_HOME/include $FCFLAGS"
    LDFLAGS="-L$FFTW_HOME/lib $LDFLAGS"
])

AS_IF([test "x$NETCDF_HOME" != "x"], [
    AC_MSG_NOTICE([NETCDF_HOME set, updating paths])
    FCFLAGS="-I$NETCDF_HOME/include $FCFLAGS"
    LDFLAGS="-L$NETCDF_HOME/lib $LDFLAGS"
])

AC_CHECK_LIB([netcdff], [nf_close],[
    AC_SUBST([LIBNC], ["-lnetcdff"])
],[
    AC_MSG_ERROR([netcdf test failed])
])

my_save_fcflags="$FCFLAGS"
FCFLAGS="$FCFLAGS -fallow-argument-mismatch"
AC_MSG_CHECKING([whether FC supports -fallow-argument-mismatch])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="-fallow-argument-mismatch"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS="$my_save_fcflags $AM_FCFLAGS -O3"
AC_MSG_CHECKING([whether FC supports -O3])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="$AM_FCFLAGS -O3"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS="$my_save_fcflags $AM_FCFLAGS -g"
AC_MSG_CHECKING([whether FC supports -g])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [AM_FCFLAGS="$AM_FCFLAGS -g"],
    [AC_MSG_RESULT([no])]
)
FCFLAGS=$my_save_fcflags
AC_SUBST([AM_FCFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CHECK_HEADER([fftw3.f],[AC_MSG_NOTICE([found fftw3.f])],[AC_MSG_ERROR([could not find fftw3.f])])

AC_CHECK_LIB([fftw3], [dfftw_plan_many_dft],
    [AC_SUBST([LIBFFTW3], ["-lfftw3"])],
    [AC_MSG_ERROR(
        [fftw3 test failed (try specifying "LDFLAGS='-Lfftw3/path/lib'")])],
    [])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT