AC_INIT([BAT], 0.9.4, [], [BAT])
: ${CXXFLAGS=""}
: ${CFLAGS=""}
: ${FFLAGS=""}
AM_INIT_AUTOMAKE()
AM_CONFIG_HEADER(config.h)
AC_PROG_CXX
AM_PROG_LIBTOOL
AC_HEADER_STDC
AC_FUNC_UTIME_NULL

dnl Get filename extension for dynamic libraries:
SHLIBEXT="${shrext_cmds}"
AC_SUBST(SHLIBEXT)

DEP_CFLAGS=""
DEP_LIBS=""

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
dnl Check for parallelization mode
dnl
AC_ARG_ENABLE(
	[parallelization],
	[AS_HELP_STRING([--enable-parallelization],[compile with OpenMP thread parallelization (default no)])],
	[
		case "${enableval}" in
			yes) parallelization=1 ;
				CXXFLAGS+=" -fopenmp" ;
				LDFLAGS+=" -fopenmp" ;
				DEP_CFLAGS="$DEP_CFLAGS -fopenmp" ;
				DEP_LIBS="$DEP_LIBS -fopenmp" ;;
			no)  parallelization=0 ;;
			*) AC_MSG_ERROR([bad value ${enableval} for --enable-parallelization]) ;;
		esac
	],
	[parallelization=0]
)

AC_DEFINE_UNQUOTED([THREAD_PARALLELIZATION], [$parallelization], [OpenMP thread parallelization])
AM_CONDITIONAL([THREAD_PARALLELIZATION], [test $parallelization = 1])

dnl
dnl Check for Cuba
dnl
cuba_min_version=3.3
AC_DEFINE_UNQUOTED([CUBA_MIN_VERSION], [$cuba_min_version], [Minimum cuba version])
AC_ARG_WITH(
	[cuba],
	[AS_HELP_STRING([--with-cuba[[=DIR]]],[compile with Cuba support (default no). Optional DIR is top of Cuba installation. If DIR is not given, the location is inferred from the partview executable.])],
	[
		use_cuba=yes
		echo "compiling with Cuba support: requires version $cuba_min_version or later"
	],
	[use_cuba=no]
)

AC_ARG_WITH(cuba_include_dir, AS_HELP_STRING([--with-cuba-include-dir=DIR    path to cuba.h file]))
AC_ARG_WITH(cuba_lib_dir,     AS_HELP_STRING([--with-cuba-lib-dir=DIR        path to libcuba.a file]))

if test "x$use_cuba" != xno; then
	if test "x$with_cuba" == xyes ; then
		AC_PATH_PROG(cuba_bin_dir, [partview])
		if test -f "$cuba_bin_dir" ; then
			# up two levels: /usr/local/bin/partview -> /usr/local
			with_cuba="$(dirname "$(dirname "$cuba_bin_dir")")"
			AC_MSG_CHECKING([for cuba in $with_cuba])
			AC_MSG_RESULT(yes)
		fi
	fi

	dnl check to see if lib and include can be determined automatically from $with_cuba
	dnl if cuba just compiled, everything is in one directory

	if test -z "$with_cuba_include_dir" ; then
		if test -d "$with_cuba/include" ; then
			with_cuba_include_dir="$with_cuba/include"
		elif test -d $with_cuba ; then
			with_cuba_include_dir="$with_cuba"
		fi
	fi

	if test -z "$with_cuba_lib_dir" ; then
		if test -d "$with_cuba/lib" ; then
			with_cuba_lib_dir="$with_cuba/lib"
		elif test -d $with_cuba ; then
			with_cuba_lib_dir="$with_cuba"
		fi
	fi

	cubaincflags=""
	if test -n "$with_cuba_include_dir"; then
		echo "checking for cuba.h in $with_cuba_include_dir"
		cubaincflags=-I$with_cuba_include_dir
		CPPFLAGS+=" $cubaincflags"
		CXXFLAGS+=" $cubaincflags"
	else
		echo "checking for cuba.h in system default directories"
	fi

	AC_CHECK_HEADERS(cuba.h, [], [AC_MSG_ERROR([Could not find cuba.h; use --with-cuba[[-include-dir]]=DIR or compile without cuba])])

	cubaldflags=""
	if test -n "$with_cuba_lib_dir"; then
	echo "checking for cuba library in path $with_cuba_lib_dir"
	cubaldflags="-L$with_cuba_lib_dir"
	LDFLAGS+=" $cubaldflags"
	else
	echo "checking for cuba library in system default directories"
	fi

	dnl Check for the right CUBA interface
	AC_COMPILE_IFELSE([
	AC_LANG_SOURCE(
	[[
	#include <cuba.h>
	int integrand(const int * ndim, const double xx[], const int * ncomp, double ff[], void * userdata){return 0;}

	int main()
	{
		int * nregions, *fail, *niterations;
		double * integral, *error, *prob;
		Suave(
		1, 1,
		&integrand, 0, 1/*nvec*/,
		0.1, 0.1,
		0, 12345,
		100, 200,
		1, 2,
		"",
		nregions, niterations, fail,
		integral, error, prob
		);

		return 0;
	}
	]]
	)],
	[cuba_version=3.3
	AC_DEFINE([CUBAVERSION_33], [1], [cuba interface version 3.3 detected])],
	[])

	if test x$cuba_version = x; then
	AC_COMPILE_IFELSE([
	AC_LANG_SOURCE(
	[[
	#include <cuba.h>
	int integrand(const int * ndim, const double xx[], const int * ncomp, double ff[], void * userdata){return 0;}

	int main()
	{
		int * nregions, *fail, *niterations;
		double * integral, *error, *prob;
		Suave(
		1, 1,
		&integrand, 0, 1,
		0.1, 0.1,
		0, 12345,
		100, 200,
		1, 2,
		"", 0, /* spinning cores */
		nregions, niterations, fail,
		integral, error, prob
		);

		return 0;
	}
	]]
	)],
	[cuba_version=4.0
	AC_DEFINE([CUBAVERSION_40], [1], [cuba interface version 4.0 detected])],
	[])
	fi

	if test x$cuba_version = x; then
		AC_MSG_ERROR([This version of cuba does not match
		the interface of the supported version $cuba_min_version or later])
	else
		echo detected cuba matching the interface of version $cuba_version
	fi

	AC_CHECK_LIB(cuba,main,[],[AC_MSG_ERROR([Could not find libcuba.a; use --with-cuba-lib-dir=DIR or compile without cuba])])
	AC_SUBST(CUBAINCLUDE,[$cubaincflags])
	AC_SUBST(CUBALIB,[$cubaldflags])

	dnl passed all tests, so build with cuba
	AM_CONDITIONAL([CUBA], [test 1 -eq 1])
else
	echo compiling without cuba
	AM_CONDITIONAL([CUBA], [test 1 -eq 0])
fi

dnl
dnl Check for ROOT
dnl
AC_DEFUN(
	[REQUIRE_ROOT_VERSION],
	[require_ROOT_version=$1]
)
REQUIRE_ROOT_VERSION([5.22])

dnl RooStats
AC_ARG_ENABLE(
	[roostats],
	[AS_HELP_STRING([--enable-roostats],[compile with RooStats support (default no)])],
	[
		use_roostats=yes
		REQUIRE_ROOT_VERSION([5.27.04])
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
dnl add some libraries to set of basic root libraries
ROOTLIBS="$ROOTLIBS -lMinuit"

dnl mathmore is a requirement of roostats
AC_MSG_CHECKING(whether ROOT is compiled with MathMore support)
hasmathmore=`$ROOTCONF --has-mathmore`
AC_MSG_RESULT($hasmathmore)
AM_CONDITIONAL([ROOTMATHMORE], [test x$hasmathmore = xyes])

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
if test $use_roostats = yes; then
	ROOTLIBS="$ROOTLIBS -lRooFitCore -lRooFit -lRooStats -lFoam"
fi

DEP_CFLAGS="$DEP_CFLAGS $ROOTAUXCFLAGS $ROOTCFLAGS"
DEP_LIBS="$DEP_LIBS -L$ROOTLIBDIR $ROOTLIBS"

AC_SUBST(DEP_CFLAGS)
AC_SUBST(DEP_LIBS)

PKG_LIBS="-lBAT -lBATmodels -lBATmtf -lBATmvc"
AC_SUBST(PKG_LIBS)

AC_CHECK_PROGS(PKGCONFIG, pkg-config, false)
AM_CONDITIONAL([COND_PKGCONFIG], [test "$PKGCONFIG" != "false"])

dnl files that may not be part of the release but exist only in a checkout of the source repository
if test -f "$srcdir/doc/introduction/versions.tex.in" ; then
    AC_OUTPUT(doc/introduction/versions.tex)
fi

dnl directories needed to store output files (e.g. from unit tests) in the proper place
AC_DEFINE_UNQUOTED([BAT_BUILDDIR], "$ac_pwd", [Absolute path to the build directory])
AC_DEFINE_UNQUOTED([BAT_SRCDIR], "`readlink -f $ac_pwd/$srcdir`", [Absolute path to the top-level source directory])
AC_DEFINE_UNQUOTED([BAT_TESTDIR], "$ac_pwd/test/", [Absolute path to the test directory])

AC_OUTPUT(
	Makefile
	bat.pc
	bat-config
	src/Makefile
	models/base/Makefile
	models/mtf/Makefile
	models/mvc/Makefile
	examples/Makefile
	examples/basic/binomial/Makefile
	examples/basic/combination1d/Makefile
	examples/basic/combination2d/Makefile
	examples/basic/errorpropagation/Makefile
	examples/basic/poisson/Makefile
	examples/basic/rootOutput/Makefile
	examples/advanced/polynomialfit/Makefile
	examples/advanced/roointerface/Makefile
	examples/advanced/mtf/systematics/Makefile
	examples/advanced/mvc/PPP/Makefile
	examples/advanced/mvc/oneObservable/Makefile
	examples/advanced/mvc/nuisance/Makefile
	examples/advanced/mvc/inverseProblem/Makefile
	examples/advanced/mvc/constraints/Makefile
	examples/advanced/referencecounting/Makefile
	examples/expert/TrialFunction/Makefile
	test/Makefile
	tools/MakefileTemplate
)
