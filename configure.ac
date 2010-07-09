AC_INIT(src/BCModel.cxx)
: ${CXXFLAGS=""}
: ${CFLAGS=""}
: ${FFLAGS=""}
AM_INIT_AUTOMAKE(BAT, 0.4.pre5)
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
		AC_SUBST(CUBAINCLUDE,[-I$cuba_cv_pathdir])
		AC_CHECK_HEADER($cuba_cv_pathdir/cuba.h, ,[echo Could not find cuba.h; echo Compiling without cuba])
	else
		echo "checking for cuba.h in system default directories"
		AC_CHECK_HEADERS(cuba.h, ,[echo Could not find cuba.h; echo Compiling without cuba])
		AC_SUBST(CUBAINCLUDE,[])
	fi
else
	echo Could not find libcuba.a
	echo Compiling without cuba
fi
AM_CONDITIONAL([CUBA], [test x$with_cuba_lib = xyes])

dnl
dnl Check for ROOT
dnl
AC_DEFUN(
	[REQUIRE_ROOT_VERSION],
	[require_ROOT_version=$1]
)
REQUIRE_ROOT_VERSION([5.16])
dnl RooStats
AC_ARG_WITH(
	[roostats],
	[AS_HELP_STRING([--with-roostats],[compile with RooStats support (default no)])],
	[
		use_roostats=yes
		REQUIRE_ROOT_VERSION([5.24])
		echo "compiling with RooStats support: required ROOT version $require_ROOT_version or later"
	],
	[use_roostats=no]
)
dnl ROOT
m4_include(tools/build/root.m4)
ROOT_PATH(
	$require_ROOT_version,
	,
	[
		AC_MSG_ERROR([ROOT cannot be found. Please ensure that it is installed
			and that either \$ROOTSYS is defined corectly or root-config is in the path.
			One can pass --with-rootsys=/path/to/rootsys to configure as well.])
		exit
	]
)
dnl RooFit
if test $use_roostats = yes; then
	HAS_ROOSTATS(
		,
		[
			AC_MSG_ERROR([ROOT was compiled without RooFit support.
				Recompile ROOT or provide path to different version.])
			AC_MSG_ERROR([Compiling BAT without RooStats support.])
			use_roostats=no
		]
	)
fi
AM_CONDITIONAL([ROOSTATS], [test x$use_roostats = xyes])

AC_OUTPUT(Makefile src/Makefile models/Makefile tools/CreateProject.sh,[chmod a+x tools/CreateProject.sh])


