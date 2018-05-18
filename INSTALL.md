Installation instructions {#cha-install}
=========================

This document provides a short description of how to compile and use
BAT on your computer.

Platforms
----------

BAT has been developed on Linux. The installation, unit tests, and
examples are run and known to work on Linux and Mac OS X. On Linux, we
test with gcc and on Mac OS X we use clang but both compilers should
work on either platform.

Windows is not supported.

Dependencies
-------------

It is understood that all commands shown below are to be entered into
a terminal.

### Required: Basic tools

BAT itself uses only C++03 features. Compilation and tests work
fine with gcc >= 4.3 and clang >= 3.3. But recent versions of ROOT
(see below) may require a C++11 compliant compiler.

Under Debian or Ubuntu, you can install the essential requirements with

    sudo apt-get install build-essential curl

In order to use the development version of BAT instead of an official
release, some more packages are needed

    sudo apt-get install autoconf automake git-core libtool

Building and installing works with autoconf >= 2.63 and automake >= 1.10. To run
the tests, a more recent automake version is needed, v1.15 is known to be
sufficient.

### Required: ROOT

ROOT is an object-oriented data-analysis framework. At http://root.cern.ch/, you
can obtain the source code as well as binary distributions for a number of Linux
distributions and Mac OS X versions. We advise to download the latest production
release of ROOT. BAT is compatible with ROOT >=5.34.19 and ROOT 6. We regularly
run unit tests with ROOT 5 and ROOT 6 to ensure backward compatibility.

On Linux, an alternative is to check your package manager for the availability
of ROOT packages. Usually these packages are rather old but often they are good
enough to build BAT. For example on Ubuntu systems up to 16.04, you can
conveniently install the requirements with

    sudo apt-get install libroot-graf2d-postscript-dev libroot-graf3d-g3d-dev\
                         libroot-math-foam-dev libroot-math-minuit-dev\
                         libroot-math-physics-dev libroot-math-mathmore-dev\
                         libroot-roofit-dev root-system-bin

#### Note

For the interface to RooFit/RooStats, ROOT must be compiled with
support for RooFit and MathMore, the latter relies on the GNU
scientific library (GSL).

### Optional: Cuba

Cuba is a library containing general-purpose multidimensional
integration algorithms. It can be obtained from
http://www.feynarts.de/cuba/. BAT is compatible with Cuba versions 3.3
through at least 4.2.

Cuba is not necessary to run BAT. We recommend it for model comparison
where expensive integrals are needed. Cuba provides integration
routines tuned for performance, which are useful for integration in
problems with not too many dimensions (~10). By default, Cuba will
evaluate in parallel and take all idle cores; the number of cores can
be set through an environment variable. For a single core, set

    CUBACORES=1

The recommended way to get Cuba is to configure BAT with the option

    --with-cuba=download

This will download a compatible version of Cuba to the local subdirectory
`external/cuba-VERSION`, compile it, and configure BAT to use it.

If you want to compile Cuba manually, make sure it is built with
position-independent code:

    ./configure CFLAGS='-fPIC -O3 -fomit-frame-pointer -ffast-math -Wall'
    make
    make install

Building
----------------------

### Obtaining BAT

You can download the latest release of BAT from
http://mpp.mpg.de/bat/.  Open a terminal, unpack the tarball usually
named like BAT-x.x.tar.gz (here x.x is the version number) and switch
to the directory

    tar -xzf BAT-x.x.tar.gz
    cd BAT-x.x

Alternatively, you can clone the git repository https://github.com/bat/bat (we
recommend using the master branch):

    git clone https://github.com/bat/bat
    cd bat
    ./autogen.sh

Now start the configuration with

    ./configure

This will check your system for all components needed to compile BAT
and set up the paths for installation. You can add the option
`--prefix=/path/to/install/bat` to `./configure`. The BAT library
files will then be installed to `$prefix/lib` and the include files to
`$prefix/include`. The default installation prefix is `/usr/local`,
which requires super-user privileges.

You can list all available options using

    ./configure --help

In the following, we describe the most useful options in detail.

### ROOT

The configure script checks for ROOT availability in the system and
fails if ROOT is not installed. You can specify the `ROOTSYS` directory
using `--with-rootsys=/path/to/rootsys`

BAT support for RooFit/RooStats is turned off by default. The feature
can be turned on explicitly with `--enable-roostats`. The configure
script will check whether the version of ROOT is sufficient and
whether ROOT was compiled with RooFit/RooStats support.

### openMP

Support for openMP threading to run multiple Markov chains in parallel
is available through the configure option `--enable-parallel`;
it is disabled by default. This requires a version of gcc accepting
the `-fopenmp` flag, anything >= 4.2 should suffice.  Note that if
threads are enabled, the default number of threads actually used is
implementation dependent and may also depend on the current load of
the CPU. Manual control over the number of threads is achieved
entirely by openMP means such as setting the environment variable
`OMP_NUM_THREADS` before running an executable.

The default version of clang does not implement openMP.

### Cuba

If you configured BAT with the option `--with-cuba=download`, BAT will
download, compile, and use Cuba automatically.  For manual
configuration, use the configure option `--with-cuba[=DIR]` to enable
Cuba. If you installed Cuba including the `partview` executable, the
Cuba installation path will be derived from its location. Otherwise,
the configure script will search for `libcuba.a` and `cuba.h` in the
system paths. If you manually specify the Cuba install path as `DIR`,
configure will look in `DIR/lib/` and `DIR/include/` instead.  For
more fine-grained control, use
`--with-cuba-include-dir=/path/to/cuba/header` and
`--with-cuba-lib-dir=/path/to/cuba/lib`.

### Advanced options

If you want to be able to step through BAT line by line with a
debugger, use `--enable-debug`. This slows down execution as it turns
off code optimization but it improves the compilation time. Another
way to speed up the build is to create only shared libraries if you
don't need static libraries: `--disable-static`. Finally, you can
reduce the output to the terminal with `--enable-silent-rules`.

### Compile

After a successful configuration, run

    make
    make install

to compile and install BAT. Note that depending on the setting of the
installation prefix you might need super-user privileges to be able to
install BAT and run `sudo make install` instead of plain `make
install`. In the former case, you might need to run `sudo ldconfig`
just once to help the loader pick up the new libraries immediately.

System setup
------------

After installation, BAT offers two mechanisms to make BAT available:

1. The script `bat-config` returns details of the BAT installation
   directories and compilation settings; see `bat-config`.

2. The file `bat.pc` contains the same information as above and can be
   used by the more powerful `pkg-config`; e.g.,

        pkg-config --modversion bat
        pkg-config --libs bat

If you do not install BAT to the system directories, you need to
manually add the path to `bat-config`, `bat.pc`, the libraries, and
the include files to the search paths. Depending on your shell, the
set of commands on linux for bash-compatible shells is

    BATPREFIX="/bat/install/prefix"
    export PATH="$BATPREFIX/bin:$PATH"
    export LD_LIBRARY_PATH="$BATPREFIX/lib:$LD_LIBRARY_PATH"
    export CPATH="$BATPREFIX/include:$CPATH"
    export PKG_CONFIG_PATH="$BATPREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"

and for csh-compatible shells is

    set BATPREFIX = /bat/install/prefix
    setenv PATH              "${BATPREFIX}/bin:${PATH}"
    setenv LD_LIBRARY_PATH   "${BATPREFIX}/lib:${LD_LIBRARY_PATH}"
    setenv CPATH             "${BATPREFIX}/include:${CPATH}"
    setenv PKG_CONFIG_PATH   "${BATPREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}"

If you want to make BAT permanently available, add the above commands
to your login script, for example to `.profile` or to `.bashrc`.

On Mac OS X you do not have to set up `LD_LIBRARY_PATH` because we use
the `rpath` option to make BAT compatible with the SIP feature enabled
by default on Mac OS X starting with El Capitan.

Updating `$CPATH` is required if you work with interactive ROOT macros
that use BAT (both for ROOT 5 and 6).

The minimal setup does not require setting `PKG_CONFIG_PATH` to run
BAT unless you want to integrate BAT into another probject using
`pkg-config`. BAT itself does not use `pkg-config`.

Including BAT in your project
-----------------------------

The most basic way to compile and link a file `example.cxx` with BAT is

    gcc `bat-config --cflags` `bat-config --libs` example.cxx -o

In a makefile, simply query `bat-config` to set appropriate
variables. However, there will be an error at runtime, for example in
interactive ROOT macros, if

    libBAT.so, libBATmodels.so, libBATmtf.so,
    libBAT.rootmap, libBATmodels.rootmap, libBATmtf.rootmap

are not in the directories searched by the library loader; see above how
to setup the `LD_LIBRARY_PATH` and the `CPATH`.

Interactive ROOT macros
-----------------------

Due to problems in ROOT 6.02.00, it is important to create an instance
of a BAT class before calling any free function defined in the BAT
libraries. Else `cling` will emit confusing
[error messages](https://github.com/bat/bat/issues/5). For example,
the right order would be

    int main() {
        BCLog::OpenLog("log.txt");
        BCAux::SetStyle();
        ...
    }

instead of the other way around around because `OpenLog` creates a
singleton object.

Contact
-------

Please consult the BAT web page http://mpp.mpg.de/bat/ for further
information. In case of questions or problems, please don't hesitate
to create an issue at https://github.com/bat/bat/issues/ or contact
the authors directly via email through bat@mpp.mpg.de.
