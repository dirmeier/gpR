dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[gpR],
	[0.1],
	[simon.dirmeier@gmx.de],
	[gpR],
	[https://github.com/dirmeier/gpR])
AC_CONFIG_SRCDIR([src/gpr.f90])

AC_MSG_NOTICE([building from R])
if test -z "${R_HOME}";
then
 	AC_MSG_NOTICE([no R_HOME found])
	FC=`"${R_HOME}/bin/R" CMD config FC`
fi

AC_PROG_FC
gpr_FLAGS=

case "${FC}" in
    *gfortran*)
      gpr_FLAGS="-std=f95 -Wall"
      ;;
    ifort)
      gpr_FLAGS="-fpp -fast"
      ;;
    *)
      AC_MSG_ERROR([Unsupported Fortran compiler])
      ;;
esac

AC_SUBST([gpr_FLAGS])
AC_OUTPUT([src/Makevars])

AC_MSG_RESULT([
	Have fun with: $PACKAGE_NAME

	FC:        $FC
])
