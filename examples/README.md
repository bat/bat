BAT examples
============

The examples are meant to illustrate the features of BAT. You are
welcome to modify them in order to become acquainted with BAT.  Each
example comes with its own README and more detailed instructions.
There are two sets of examples, `basic/*` and `advanced/*`. If you are
new to BAT, try `basic/poisson`.

You can also find tutorials on the BAT webpage http://mpp.mpg.de/bat

Running the examples
--------------------

Most examples need to be compiled while the ROOT macros show how to
use BAT from within the ROOT interpreter through CINT (ROOT 5) or
cling (ROOT 6).

The standard way to run the examples requires that BAT was
successfully installed in the system. In particular, the installation
directory is queried from `bat-config`. Instructions on how to run the
examples are either inside the source code (for ROOT macros) or in the
README file (for compiled examples).

For testing purposes, we compile and run the examples via `make check`
in the build directory; i.e., no installation is required here. Issued
from the base build directory, all examples are executed. You can also
run `make check` inside an individual example directory to run just
that one example. Note, however, that you won't get the output
interactively in the shell but it is captured in the respective
`example*.log` file.


Overview of examples
----------------------

### basic

Directoy:     binomial
Type:         Compiled program
Description:  A simple example for a binomial experiment. The
              efficiency for a certain process is estimated given a
              total and selected number of "events" or observations.

Directory:    combination1d
Type:         Compiled program
Description:  A simple example for the combination of two
              measurements where one is the result of an earlier
              measurement.

Directory:    combination2d
Type:         Compiled program
Description:  A simple example for the combination of two
              measurements with two variables where one is the
              result one an earlier measurement.

Directory:    efficiencyFitter
Type:         ROOT macro
Description:  Shows how to fit a histogram ratio where
              the numerator histogram is a subset of the denominator
              histogram. Binomial uncertainties are used.

Directory:    errorPropagation
Type:         Compiled program
Description:  The probability distribution for the ratio of two
              random variables is calculated.

Directory:    graphFitter
Type:         ROOT macro
Description:  Straight line fit using GraphFitter assuming gaussian
              uncertainties.

Directory:    histogramFitter
Type:         ROOT macro
Description:  Histogram fit using poissonian uncertainties.

Directory:    poisson
Type:         Compiled program
Description:  A simple example for a Poisson counting experiment. The
              expectation value for a Poisson distribution is
              estimated given a certain number of observations.

Directory:    rootOutput
Type:         Compiled program
Description:  A simple 2D-Gaussian posterior. The example shows how to
              write the MCMC to a ROOT file.

### advanced

Directory:    advancedGraphFitter
Type:         ROOT macro
Description:  Graph fits using 4 different functions assuming
              gaussian uncertainties. This example represents the
              example from the BAT paper showing a comparison
              between different fits.

Directory:    mtf
Type:         ROOT macros and compiled programs
Description:  Examples for how to fit several templates to a data set
              using 1-dim. histograms. The examples feature the
              BCMTF which is included in the models directory. This
              code will replace the BCTemplateFitter in future
              versions.

Directory:    polynomialFit
Type:         Compiled program
Description:  Example showing how to define own complicated models and build
              a program for running with BAT. It shows fits using
              two functions when assuming asymmetric uncertainties.

Directory:    referenceCounting
Type:         Compiled program
Description:  Compute the reference prior for a counting experiment and show
              how it gets updates with data.

Directory:    rooInterface
Type:         Compiled program and ROOT macros
Description:  Examples of usage of the RooInterface and BATCalculator.
              The former shows how BAT analysis can be performed on a
              RooFit workspace. The later shows how BAT can be run from
              within RooStats.

Directory:    trialFunction
Type:         Compiled program
Description:  A simple 2D-Gaussian posterior. The example shows how to
              change the trial function for the MCMC.
