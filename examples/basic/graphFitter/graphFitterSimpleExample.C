// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a graph using a function defined by the user. TGraphErrors
// with defined y uncertainties has to be defined. In the fit the
// uncertainties are considered to be Gaussian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x graphFitterSimpleExample.C
//
// or
//
//    root[1] .L graphFitterSimpleExample.C
//    root[2] graphFitterSimpleExample()
//
// or from the command line
//
//    $ root graphFitterSimpleExample.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x graphFitterSimpleExample.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro
// the #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#define COMPILER (!defined(__CINT__) && !defined(__CLING__))

#if defined(__MAKECINT__) || defined(__ROOTCLING__) || COMPILER

#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCGraphFitter.h>

#endif

// The macro performs a straight line fit to the data. The fit function
// is defined using the ROOT TF1 object and the data to fit are stored
// in the TGraphErrors object. The y-errors have to be defined in the
// TGraphErrors for the fit to work.

// ---------------------------------------------------------
void graphFitterSimpleExample()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // -------------------------
    // define a fit function, also used to create data
    TF1 f1("f1", "[0] + [1]*x", 0.0, 100.0);
    f1.SetParNames("a0", "a1");

    // Parameter limits must be defined for every parameter
    f1.SetParLimits(0, -20.0, 30.0);
    f1.SetParLimits(1,   0.5,  1.5);
    // -------------------------

    // -------------------------
    // Create data
    // initialize random number generator
    TRandom3 random(1234);

    // The data fitted are generated randomly generated from f1 defined above
    f1.SetParameters(11.0, 1.0);

    // with Gaussian smearing parameterized by "sigma":
    double sigma  =  5.0;

    // create graph
    TGraphErrors graph;

    // fill it with 10 points (between 0 and 100),
    // using sigma as y uncertainty
    int N = 10;
    for (int i = 0; i < N; ++i) {
        double x = 100 * (i + 0.5) / N;
        graph.SetPoint(i, x, random.Gaus(f1.Eval(x), sigma));
        graph.SetPointError(i, 0, sigma);
    }
    // -------------------------

    // create a new graph fitter
    BCGraphFitter gf(graph, f1);

    // set Metropolis as marginalization method
    gf.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    gf.SetPrecision(BCEngineMCMC::kQuick);

    // perform the fit
    gf.Fit();

    // print data and the fit
    TCanvas c1("c1");
    gf.DrawFit("", true);
    c1.Print("fit.pdf");

    // print marginalized distributions
    gf.PrintAllMarginalized("distributions.pdf");

    // print results
    gf.PrintSummary();
}
