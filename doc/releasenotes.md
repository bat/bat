BAT - Bayesian Analysis Toolkit
===============================

Release notes for version:    1.0
Release date:                 2015-12-XX
Urgency:                      high

Overview
--------

This release is the result of a major clean up of the code and build
system. This allowed us to implement important new features but
implied that we had to break compatibility with previous releases. We
removed many inconsistencies in the interface and intend to not do any
such changes in future versions 1.X of BAT. To ease the maintenance
burden, we also removed some parts of the code that felt were untested
and not easy to keep up to date.

New features and improvements
---------------------------

### proposal function

We now have a multivariate proposal function that learns its
covariance during the prerun, it's the new default. It provides a
massive speed-up compared to the factorized proposal for problems with
more than a few parameters. Both factorized and multivariate proposal
are of the Gaussian or Student's [default] type. For Student's, the
degree of freedom is adjustable and defaults to 1; that is a
Cauchy. The mixing of chains via R values has been synced with an
update to the original paper by Gelman and Rubin.

###

### Miscellaneous

* update CUBA default settings
* use exceptions instead of an error message when there is no way to
  proceed

### testing and building

* `INSTALL.md` gives explicit advice on how to install BAT in the
  current Ubuntu and Debian releases
* we now run tests for every change to the master branch with root 5
  and 6 on Ubuntu and Mac to ensure compatibility with all supported
  platforms
* `./configure --with-cuba=download` downloads and builds the latest
  version of CUBA 4.2 and sets everything up so BAT can use it without
  any further setup required from the user
* `make` builds the library and all examples
* unit tests expanded
* `make check` runs the unit tests and all examples, without requiring
  to install BAT
* out-of-source builds now work; e.g., in a subdirectory of the source
  files, one can build BAT such that the source directory remains
  clean from build artefacts
* `make distcheck` does the above and more and now passes
* the performance test suite is overhauled

Bug fixes
---------

* `BC*Fitter` classes are now thread safe
* the output precision is now adaptively calculated
* GetCurrentChain() now returns robust results with multiple threads
* MCMCInitialize() called at the appropriate time
* LogApproxBinomial returns -inf when needed
* outputs are polished
* return -inf instead inf in `BCConstantPrior`

Known problems
--------------

* ROOT 6.05/02 has a bug in the TF1 copy constructor slowing down the
  fast fitter and leading to segfaults with multiple threads; see
  issues #74 and #81
* with ROOT 5.34 and CUBA, `BCIntegrate.TEST` produces errors on the
  shell but the program passes all checks. This does not occur with
  ROOT 6.
* `BCRooInterface` and `BATCalculator` are not thread safe; both
  classes have been updated to compile and run but we can offer no
  guarantee that they work correctly because the original author is
  not available anymore. In contrast to all other BAT classes,
  `BATCalculator` can not be used in interactive ROOT 6 sessions.
* `BCTH1Prior::swap` doesn't work correctly with root older than 5.34/25

Interface changes
-----------------

* `./configure --enable-parallelization` -> `./configure --enable-parallel`
* many methods now accept or return `std::string` references instead
  of `const char*`
* rename many getters and setters in `BCEngineMCMC` striving for
  consistency and simplicity; for example, `MCMCSetNIterationsRun` ->
  `SetNIterationsRun`
