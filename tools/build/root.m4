dnl -*- mode: autoconf -*-
dnl
dnl $Id: root.m4,v 1.1 2009-01-11 21:53:12 dkollar Exp $
dnl $Author: dkollar $
dnl $Date: 2009-01-11 21:53:12 $
dnl
dnl Autoconf macro to check for existence or ROOT on the system
dnl Synopsis:
dnl
dnl  ROOT_PATH([MINIMUM-VERSION, [ACTION-IF-FOUND, [ACTION-IF-NOT-FOUND]]])
dnl
dnl Some examples:
dnl
dnl    ROOT_PATH(3.03/05, , AC_MSG_ERROR(Your ROOT version is too old))
dnl    ROOT_PATH(, AC_DEFUN(HAVE_ROOT))
dnl
dnl The macro defines the following substitution variables
dnl
dnl    ROOTCONF           full path to root-config
dnl    ROOTEXEC           full path to root
dnl    ROOTCINT           full path to rootcint
dnl    ROOTLIBDIR         Where the ROOT libraries are
dnl    ROOTINCDIR         Where the ROOT headers are
dnl    ROOTCFLAGS         Extra compiler flags
dnl    ROOTLIBS           ROOT basic libraries
dnl    ROOTGLIBS          ROOT basic + GUI libraries
dnl    ROOTAUXLIBS        Auxilary libraries and linker flags for ROOT
dnl    ROOTAUXCFLAGS      Auxilary compiler flags
dnl    ROOTRPATH          Same as ROOTLIBDIR
dnl
dnl The macro will fail if root-config and rootcint isn't found.
dnl
dnl Christian Holm Christensen <cholm@nbi.dk>
dnl

AC_DEFUN([ROOT_PATH],
[
  AC_ARG_WITH(rootsys,
  [  --with-rootsys=DIR      top of the ROOT installation directory],
    user_rootsys=$withval,
    user_rootsys="none")
  if test ! x"$user_rootsys" = xnone; then
    rootbin="$user_rootsys/bin"
  elif test ! x"$ROOTSYS" = x ; then
    rootbin="$ROOTSYS/bin"
  else
   rootbin=$PATH
  fi
  AC_PATH_PROG(ROOTCONF, root-config , no, $rootbin)
  AC_PATH_PROG(ROOTEXEC, root , no, $rootbin)
  AC_PATH_PROG(ROOTCINT, rootcint , no, $rootbin)

  if test ! x"$ROOTCONF" = "xno" && \
     test ! x"$ROOTCINT" = "xno" ; then

    # define some variables
    ROOTLIBDIR=`$ROOTCONF --libdir`
    ROOTINCDIR=`$ROOTCONF --incdir`
    ROOTCFLAGS=`$ROOTCONF --noauxcflags --cflags`
    ROOTLIBS=`$ROOTCONF --noauxlibs --noldflags --libs`
    ROOTGLIBS=`$ROOTCONF --noauxlibs --noldflags --glibs`
    ROOTAUXCFLAGS=`$ROOTCONF --auxcflags`
    ROOTAUXLIBS=`$ROOTCONF --auxlibs`
    ROOTRPATH=$ROOTLIBDIR

    if test $1 ; then
      AC_MSG_CHECKING(whether ROOT version >= [$1])
      vers=`$ROOTCONF --version | tr './' ' ' | awk 'BEGIN { FS = " "; } { printf "%d", ($''1 * 1000 + $''2) * 1000 + $''3;}'`
      requ=`echo $1 | tr './' ' ' | awk 'BEGIN { FS = " "; } { printf "%d", ($''1 * 1000 + $''2) * 1000 + $''3;}'`
      if test $vers -lt $requ ; then
        AC_MSG_RESULT(no)
        no_version="yes"
        no_root="yes"
      else
        AC_MSG_RESULT(yes)
      fi
    fi
  else
    # otherwise, we say no_root
    no_root="yes"
  fi

  AC_SUBST(ROOTLIBDIR)
  AC_SUBST(ROOTINCDIR)
  AC_SUBST(ROOTCFLAGS)
  AC_SUBST(ROOTLIBS)
  AC_SUBST(ROOTGLIBS)
  AC_SUBST(ROOTAUXLIBS)
  AC_SUBST(ROOTAUXCFLAGS)
  AC_SUBST(ROOTRPATH)

  if test "x$no_version" = "xyes" ; then
    echo "ROOT version $vers is too old."
  fi

  if test "x$no_root" = "x" ; then
    ifelse([$2], , :, [$2])
  else
    ifelse([$3], , :, [$3])
  fi
])

dnl
dnl macro to check whether ROOT is compiled with RooFit support
dnl
AC_DEFUN([HAS_ROOSTATS],
[
  AC_PATH_PROG(ROOTCONF, root-config , no, $rootbin)

  if test ! x"$ROOTCONF" = "xno" ; then

    AC_MSG_CHECKING(whether ROOT is compiled with RooFit/RooStats support)
    hasroofit=`$ROOTCONF --has-roofit`
    if test x$hasroofit = xno ; then
      AC_MSG_RESULT(no)
      no_roofit="yes"
    else
      AC_MSG_RESULT(yes)
    fi

    AC_MSG_CHECKING(whether ROOT is compiled with MathMore support)
    hasmathmore=`$ROOTCONF --has-mathmore`
    if test x$hasmathmore = xno ; then
      AC_MSG_RESULT(no)
      no_mathmore="yes"
    else
      AC_MSG_RESULT(yes)
    fi
  else
    # otherwise, we say no_root
    no_roofit="yes"
  fi

  if test "x$no_roofit" = "x" && test "x$no_mathmore" = "x"; then
    ifelse([$1], , :, [$1])
  else
    ifelse([$2], , :, [$2])
  fi
])


#
# EOF
#
