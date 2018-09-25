BAT - Bayesian Analysis Toolkit
===============================

[![Build Status](https://travis-ci.org/bat/bat.svg?branch=master)](https://travis-ci.org/bat/bat)
[![Manual](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://bat.github.io/bat-docs/master/manual/html/index.html)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1322675.svg)](https://doi.org/10.5281/zenodo.1322675)

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

Authors
--------

See `doc/CREDITS` file for list of contributors to BAT.

Availability
-------------

The latest code and official BAT releases are published online at
https://github.com/bat/bat. Older BAT releases are still available from
http://mpp.mpg.de/bat/. BAT is open-source software under the LGPLv3 or later.
See `doc/COPYING` and `doc/LICENSE` for licensing details.

Documentation
-------------

For **installation instructions**, see the
[`INSTALL.md`](https://github.com/bat/bat/blob/master/INSTALL.md) file or the
corresponding [manual
section](https://bat.github.io/bat-docs/master/manual/html/cha-install.html)

**manual**: [html](https://bat.github.io/bat-docs/master/manual/html/index.html), [pdf](https://bat.github.io/bat-docs/master/manual/BAT-manual.pdf)

**reference guide**: [html](https://bat.github.io/bat-docs/master/ref-guide/html/index.html)

Contact
-------------

For additional information and contacting the authors, please consult the BAT
web page at http://mpp.mpg.de/bat/. If you want to report an error or file a
request, please file an issue at https://github.com/bat/bat/issues.

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

Other files distributed with BAT are part of the configuration and
build system.
