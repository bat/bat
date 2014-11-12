BAT - Bayesian Analysis Toolkit
===============================

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
can be found in `doc/releasenotes.md`. The detailed list of changes
can be found in `doc/ChangeLog`.

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
* `examples/`: well commented example programs
* `BAT/`: BAT header files
* `src/`: BAT source files
* `models`: models and interfaces for BAT
* `tools`: tools useful for BAT
* `INSTALL.md`: instructions to install BAT on your system
* `README.md`: basic information about BAT

Other files distributed with BAT are part of installation and configuration
system.

Installation
-------------

See the `INSTALL.md` file for installation instructions. The
instructions for the latest development version are also available
[online](https://github.com/bat/bat/blob/master/INSTALL.md).


Contact
-------------

For additional information and contacting the authors, please consult
the BAT web page at http://mpp.mpg.de/bat/. If you want to report an
error or file a request, please file an issue at
https://github.com/bat/bat/issues.
