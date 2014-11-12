BAT - Bayesian Analysis Toolkit
===============================

Release notes for version:    0.9.4
Release date:                 2014-11-20
Urgency:                      medium

New features
-------------

* support ROOT6 (requires compiler that implements C++ 11)
* support cuba 4.0

Improvements
------------

* `CreateProject.sh` becomes `bat-project`, loses a lot of clutter, and is installed to the `bin` directory
* every subclass of BCEngineMCMC provides access to the parameters via `GetParameters`

Interface changes (for users of previous BAT versions):
-------------------------------------------------------

* cuba 3.2 had bugs and is not supported anymore
* configure option `--with-roostats` became `--enable-roostats`

Bug fixes
---------

* fix off-by-one error in BCH1D::Draw
* all examples that use the ROOT interpreter work
* bat libraries have the proper dependency on other bat libraries and external libraries encoded

Improvements of the build system
--------------------------------

* important flags from ROOT are now passed on
* fix and update installation instructions, display how-to-setup-bat steps at the end of `make install`
* support VPATH builds: you can now create a subdirectory, run `configure` from there, and keep all build products in the subdirectory
* unit tests and the latex sources now part of the distribution
* `make distcheck` passes
* dependencies properly included in the libraries
* transform INSTALL to markdown format, fix csh instructions
* add a `bat-config` script that is installed to `$PREFIX/bin/`. It
  can be queried to return libraries, compiler flags, version etc. The
  variable `$BATINSTALLDIR` is no longer used in the examples and
  `bat-project`
* support `pkgconfig`
* no longer required ROOT's `libmathmore` to build the bat libraries
