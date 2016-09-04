dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[lvgpp],
	[0.1],
	[lvgpp@simon-dirmeier.net],
	[lvgpp],
	[https://github.com/rafstraumur/lvgpp])
AC_CONFIG_SRCDIR([src/predict_lvgpr.f95])


AC_MSG_NOTICE([building from R])
if test -z "${R_HOME}"; 
then
 	AC_MSG_NOTICE([no R_HOME found])
	FC=`"${R_HOME}/bin/R" CMD config FC`
fi	

AC_PROG_FC
lvgpp_FLAGS="-std=f2008"

AC_SUBST([lvgpp_FLAGS])

AC_OUTPUT([src/Makevars])

AC_MSG_RESULT([
	Have fun with: $PACKAGE_NAME

	FC:        $FC
])
