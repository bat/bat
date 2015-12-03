// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a histogram using a function defined by the user. In the fit the
// uncertainties are considered to be poissonian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x histogramFitterExample.C
//
// or
//
//    root[1] .L histogramFitterExample.C
//    root[2] histogramFitterExample()
//
// or from the command line
//
//    $ root histogramFitterExample.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x histogramFitterExample.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro
// the #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#define COMPILER (!defined(__CINT__) && !defined(__CLING__))

#if defined(__MAKECINT__) || defined(__ROOTCLING__) || COMPILER

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCHistogramFitter.h>

#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>

#endif

// The macro performs a gaussian+constant fit to the data. The fit function
// is defined using the ROOT TF1 object and the data to fit are stored
// in the TH1D object.

// ---------------------------------------------------------
void histogramFitterExample()
{
    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // -------------------------
    // Create data
    // initialize random number generator
    TRandom3 random(1234);

    // create new histogram
    TH1D hist("data", ";x;N", 100, 0.0, 100.0);
    hist.SetStats(kFALSE);

    // fill signal, 100 events distributed by Gaussian with mean = 65, sigma = 5
    for (int i = 0; i < 100; ++i)
        hist.Fill(random.Gaus(65, 5));

    // fill background, 100 events, uniformly distributed
    for (int i = 0; i < 100; ++i)
        hist.Fill(random.Uniform() * 100);
    // -------------------------

    // -------------------------
    // Define a fit function, which is also used to generate data
    TF1 f1("f1", "[0]/sqrt(2*pi)/[2] * exp(-0.5*((x-[1])/[2])^2) + [3]", 0., 100.);
    f1.SetParNames("SignalYield", "SignalMean", "SignalSigma", "BackgroundYield");
    f1.SetParLimits(0,  0.0, 200.0);
    f1.SetParLimits(1, 55.0,  75.0);
    f1.SetParLimits(2,  0.1,  10.0);
    f1.SetParLimits(3,  0.0,   2.0);
    // -------------------------

    // create a new histogram fitter
    BCHistogramFitter hf(hist, f1);

    // set Metropolis as marginalization method
    hf.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    hf.SetPrecision(BCEngineMCMC::kQuick);

    // integrate function over bin (true) or use linear interpolation
    hf.SetFlagIntegration(false);

    // set priors
    hf.GetParameters().SetPriorConstantAll();

    // perform fit
    hf.Fit();

    // calculate p values
    hf.CalculatePValueFast(hf.GetBestFitParameters());
    cout << "p value " << hf.GetPValue() << ", corrected for degrees of freedom " << hf.GetPValueNDoF() << endl;

    // print marginalized distributions
    hf.PrintAllMarginalized("distributions.pdf");

    // print data and fit
    TCanvas c1("c1");
    hf.DrawFit("", true); // draw with a legend
    c1.Print("fit.pdf");

    return;
}
