BAT - Bayesian Analysis Toolkit
===============================

This document provides a short description of how to compile and use
BAT on your computer.

Platforms
----------

BAT has been developed on Linux machines running different
distributions and different versions of the kernel and gcc. As far as
we know there is nothing distribution dependent inside of BAT. A gcc
version >4.2 should suffice to compile the C++ code.

The installation and functionality of BAT has also been tested on MAC OS X.

Windows is not supported.

Dependencies
-------------

### Required: ROOT

ROOT is an object-oriented data-analysis framework. You can obtain it
from http://root.cern.ch/. Since BAT version 0.4.2 a ROOT version 5.22
or later is needed to compile and run BAT. ROOT 6 is supported as
well.

Please check your Linux distribution for the availability of
precompiled packages on your system. Many distributions offer the ROOT
packages, albeit older versions. For example in Ubuntu 14.04, you can
conveniently install the entire ROOT system through the package
`root-system`. However, if you rely on the optional `roostats`
interface, you may still have to compile ROOT yourself.

#### Note

For the interface to RooFit/RooStats, a ROOT version
5.27/04 or later is necessary and ROOT has to be compiled with
MathMore enabled.

### Optional: CUBA

CUBA is a library containing general purpose multidimensional
integration algorithms. It can be obtained from
http://www.feynarts.de/cuba/. BAT will compile and run with Cuba
version 3.3 or later. Cuba is not necessary to run BAT, however, its
use is recommended as it provides integration routines tuned for
performance, which are useful for integration in problems with not too
many dimensions (~15). You need to compile the CUBA library as
position-independent code to use it from BAT. For the `gcc`, install
CUBA as

    ./configure CFLAGS='-fPIC -O3 -fomit-frame-pointer -ffast-math -Wall'
    make
    make install

Building
----------------------

Unpack the tarball containing the BAT source usually named like
BAT-x.x.tar.gz (here x.x is the version number) using command

    tar -xzf BAT-x.x.tar.gz

A directory called BAT-x.x will be created containing the source code.
Enter the directory and run the configuration using commands

    cd BAT-x.x
    ./configure

This will check your system for all components needed to compile BAT
and set up the paths for installation. You can add option
`--prefix=/path/to/install/bat` to `./configure` to specify the the
prefix to the BAT installation path. The BAT library files will be
installed to `$prefix/lib` and the include files to
`$prefix/include`. The default installation prefix is `/usr/local`.

The configure script checks for ROOT availability in the system and
fails if ROOT is not installed. You can specify the `ROOTSYS` directory
using `--with-rootsys=/path/to/rootsys`

You can configure BAT with the RooFit/RooStats support using
`--enable-roostats`. The configure script will check whether the version of
ROOT is sufficient and whether the ROOT was compiled with RooFit/RooStats
support.

Support for openMP threading to run multiple Markov chains in parallel
is available through the configure option `--enable-parallelization`;
it is disabled by default. This requires a version of gcc accepting
the `-fopenmp` flag, anything >= 4.2 should suffice.  Note that if
threads are enabled, the default number of threads actually used is
implementation dependent and may also depend on the current load of
the CPU. Manual control over the number of threads is achieved
entirely by openMP means such as setting the environment variable
`OMP_NUM_THREADS` before running an executable.

You can compile BAT with Cuba support using option
`--with-cuba[=DIR]`.  The configure script will then search for
`libcuba.a` and `cuba.h` in the system paths, or, if `DIR` is given, in
`DIR` (correct if cuba just compiled in its source directory), and in
`DIR/lib/` and `DIR/include/` subdirectories. For more fine-grained
control, use `--with-cuba-include-dir=/path/to/cuba/header` and
`--with-cuba-lib-dir=/path/to/cuba/lib`.

Note that the cuba library still needs to be available to the loader
at runtime. This is accomplished, for example, by adding
`/path/to/cuba/lib` to `LD_LIBRARY_PATH`. Please compile the Cuba library
with position-independent code. For example,

```bash
./configure CFLAGS='-fPIC -O3 -fomit-frame-pointer -ffast-math -Wall'
```

to avoid linker errors like

    /usr/bin/ld: /usr/local/lib/libcuba.a(Vegas.o): relocation R_X86_64_32
    against `.rodata.str1.8' can not be used when making a shared object;
    recompile with -fPIC

You can list all available options using

    ./configure --help

After a successful configuration, run

    make
    make install

to compile and install BAT. Note that depending on the setting of
installation prefix you might need root privileges to be able to
install BAT and run `sudo make install` instead of plain 'make
install'. In the former case, you might need to run `sudo ldconfig`
just once to help the linker pick up the new libraries immediately.

System setup
------------

After installation, BAT offers two mechanisms to make BAT available:

1. The script `bat-config` returns the installation prefix, the
   libraries for linking, the C flags, and the version as

        bat-config --prefix
        bat-config --libs
        bat-config --cflags
        bat-config --version

2. The file `bat.pc` contains the same information as above and can be
   used by the more powerful `pkg-config`; e.g.,

        pkg-config --modversion bat
        pkg-config --libs bat

If you do not install BAT to the system directories, you need to
manually add the path to `bat-config`, `bat.pc`, the libraries, and to
the include files to the search paths. Depending on your shell you can
do that via the commands

```bash
BATPREFIX="/bat/install/prefix"
export PATH="$BATPREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$BATPREFIX/lib:$LD_LIBRARY_PATH"
export CPATH="$BATPREFIX/include:$CPATH"
export PKG_CONFIG_PATH="$BATPREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
```

or

```bash
set BATPREFIX = /bat/install/prefix
setenv PATH              "${BATPREFIX}/bin:${PATH}"
setenv LD_LIBRARY_PATH   "${BATPREFIX}/lib:${LD_LIBRARY_PATH}"
setenv CPATH             "${BATPREFIX}/include:${CPATH}"
setenv PKG_CONFIG_PATH   "${BATPREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}"
```

for bash and csh compatible shells, respectively. On Mac OS X you
might also need to setup `DYLD_LIBRARY_PATH`. If you want to make BAT
permanently available, add the above commands to your `.bashrc` or
`.tcshrc`.

Note that `bat-config` needs to be on the `PATH` to compile the
programs that ship with BAT in the `examples/` subdirectory.

The variable `CPATH` is required if you work with ROOT macros
that use BAT (both for ROOT 5 and 6)

Including BAT in your project
-----------------------------

The most basic way to compile and link a file `example.cxx` with BAT is

```bash
gcc `bat-config --cflags` `bat-config --libs` example.cxx -o
```

In makefile projects, simply add option for use in compiled programs
would also be to add `bat-config --cflags` to CXXFLAGS and `bat-config
--libs` to `LDFLAGS` in your Makefile. However, there will be an error
at runtime, for example in interactive ROOT macros, if

    libBAT.so, libBATmodels.so,
    libBATmtf.so, libBATmvc.so, libBAT.rootmap, libBATmodels.rootmap,
    libBATmtf.rootmap, libBATmvc.rootmap

are not in the directories found be the library loader; see above how
to setup the `LD_LIBRARY_PATH` and the `CPATH`.

Interactive ROOT macros
-----------------------

Due to problems in ROOT 6.02.00, it is important to create an instance
of a BAT class before calling any free function defined in the BAT
libraries. Else `cling` will emit confusing
[error messages](https://github.com/bat/bat/issues/5). For example,
the right order would be

```cpp
int main() {
    BCLog::OpenLog("log.txt");
    BCAux::SetStyle();
    ...
}
```

instead of the other way around around because `OpenLog` creates a singleton object.

Contact
-------

Please, consult the BAT web page http://mpp.mpg.de/bat/ for further
information. You can also contact the authors directly via email:
bat@mpp.mpg.de
