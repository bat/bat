BAT - Bayesian Analysis Toolkit
===============================

[![Build Status](https://travis-ci.org/bat/bat.svg?branch=master)](https://travis-ci.org/bat/bat)

The Bayesian Analysis Toolkit is a software analysis package which allows
to address the main goals of a typical data analysis:

 - compare model predictions with data,
 - draw conclusions on the validity of the model as a representation
   of the data and
 - extract the possible values of parameters within the context of
   a model.

The BAT is based on Bayes' Theorem and is realized with the use of Markov
Chain Monte Carlo. This gives access to the full posterior probability
distribution and enables straightforward parameter estimation, limit
setting and uncertainty propagation.

The BAT is implemented in C++ on top of CERN's ROOT with emphasis put
on flexibility and modularity in defining models while keeping in mind
the reliability and speed requirements of the numerical operations.

The information about the main features introduced in this version
can be found in `doc/releasenotes.md`.

Authors
--------

See `doc/CREDITS` file for list of contributors to BAT.

Availability
-------------

All BAT releases are available from http://mpp.mpg.de/bat/.  The
source code is managed online at https://github.com/bat/bat.  BAT is
open-source software under the GPLv3 or later.  See `doc/COPYING` and
`doc/LICENSE` for licensing details.

Contents
---------

Directory `bat/`:

* `doc/`: documentation about BAT
* `examples/`: well commented example programs, see `examples/README.md`
* `BAT/`: BAT header files
* `src/`: source files of the BAT core library
* `models`: models for specific data-analysis problems: fast fitters
  in `base/` and the multi-template fitter in `mtf`
* `benchmarks`: performance test suite with html output
* `test`: unit tests
* `tools`: template script to generate a BAT project and build utilities
* `INSTALL.md`: instructions to install BAT on your system
* `README.md`: this document

Other files distributed with BAT are part of the configuration and
build system.

Installation
-------------

See the `INSTALL.md` file for installation instructions. The
instructions for the latest development version are also available
at https://github.com/bat/bat/blob/master/INSTALL.md.


Contact
-------------

For additional information and contacting the authors, please consult
the BAT web page at http://mpp.mpg.de/bat/. If you want to report an
error or file a request, please file an issue at
https://github.com/bat/bat/issues.
