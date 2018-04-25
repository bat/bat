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
#if defined(__MAKECINT__) || defined(__ROOTCLING__) || (!defined(__CINT__) && !defined(__CLING__))

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
    BCLog::OpenLog("log.txt", BCLog::detail);

    // -------------------------
    // Create data
    // initialize random number generator
    TRandom3 random(1234);

    // create new histogram
    const unsigned nbins = 50;
    const double xmin = 0;
    const double xmax = 100;
    TH1D hist("data", ";x;N", nbins, xmin, xmax);
    hist.SetStats(kFALSE);

    // fill signal, `nbins` events following a Gaussian with mean = 65, sigma = 5
    const double mu = 65;
    const double sigma = 5;
    for (int i = 0; i < nbins; ++i)
        hist.Fill(random.Gaus(mu, sigma));

    // fill background, expect one event per bin
    for (int i = 0; i < nbins; ++i)
        hist.Fill(xmin + random.Uniform() * (xmax - xmin));

    // -------------------------
    // Define a fit function with 4 parameters
    // TF1 f1("f1", MyFunctor(xmin, xmax, nbins), xmin, xmax, 4);
    // f1.SetParNames("SignalYield", "SignalMean", "SignalSigma", "BackgroundYieldPerBin");
    TF1 f1("f1", "[0]/sqrt(2*pi)/[2] * exp(-0.5*((x-[1])/[2])^2) + [3]*0.01", 0., 100.);
    f1.SetParNames("SignalYield", "SignalMean", "SignalSigma", "BackgroundYield");
    f1.SetParLimits(0,  0.0, 2 * nbins);
    f1.SetParLimits(1, mu - 10, mu + 10);
    f1.SetParLimits(2,  0.1,  2 * sigma);
    f1.SetParLimits(3,  0.0,   2 * nbins);
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

    // calculate p value...
    double p = hf.CalculatePValueFast(hf.GetBestFitParameters());

    // and correct for the degrees of freedom to yield an approximately
    // uniformly distributed p value for the true model
    cout << "p value " << p << ", corrected for degrees of freedom " << BCMath::CorrectPValue(p, hf.GetNFreeParameters(), nbins) << endl;

    // print marginalized distributions
    hf.PrintAllMarginalized("distributions.pdf");

    // print data and fit
    TCanvas c1("c1");
    hf.DrawFit("", true); // draw with a legend
    c1.Print("fit.pdf");
}
