AC_INIT(src/BCModel.cxx)
AM_INIT_AUTOMAKE(BAT, 0.3.2.rc3)
AM_CONFIG_HEADER(config.h)
AC_PROG_CC
AC_PROG_F77
AC_PROG_CXX
AM_PROG_LIBTOOL
AC_HEADER_STDC
AC_FUNC_UTIME_NULL

dnl
dnl Check for Darwin
dnl
AM_CONDITIONAL([IS_DARWIN], [test x$(uname) = xDarwin])

dnl
dnl Check for debug mode
dnl
AC_ARG_ENABLE(
	[debug],
	[  --enable-debug          Compile in debug mode],
	[
		case "${enableval}" in
			yes) debug=true ; echo Compiling in debug mode ;;
			no)  debug=false ;;
			*) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
		esac
	],
	[debug=false]
)
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

dnl
dnl Check for Cuba
dnl
AC_ARG_WITH(cuba_path, AS_HELP_STRING([--with-cuba-path=DIR    path to cuba.h file]))
AC_CHECK_LIB(cuba,main,[with_cuba_lib=yes],[with_cuba_lib=no])
if test $with_cuba_lib = yes; then
	if test -n "$with_cuba_path"; then
		echo "checking for cuba.h in path $with_cuba_path"
		AC_CACHE_CHECK([for the directory containing cuba.h], cuba_cv_pathdir,cuba_cv_pathdir=$with_cuba_path)
		AC_SUBST(CUBADIR,$cuba_cv_pathdir)
		AC_CHECK_HEADER($cuba_cv_pathdir/cuba.h, ,[echo Could not find cuba.h; echo Compiling without cuba])
	else
		echo "checking for cuba.h in system default directories"
		AC_CHECK_HEADERS(cuba.h, ,[echo Could not find cuba.h; echo Compiling without cuba])
	fi
else
	echo Could not find libcuba.a
	echo Compiling without cuba
fi
AM_CONDITIONAL([CUBA], [test x$with_cuba_lib = xyes])

dnl
dnl Check for ROOT
dnl
m4_include(tools/build/root.m4)
ROOT_PATH(
	5.16,
	,
	[
		AC_MSG_ERROR([ROOT cannot be found. Please ensure that it is installed
			and that either \$ROOTSYS is defined corectly or root-config is in the path.
			One can pass --with-rootsys=/path/to/rootsys to configure as well.])
		exit
	]
)

AC_OUTPUT(Makefile src/Makefile)
