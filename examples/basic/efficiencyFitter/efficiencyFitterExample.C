// This ROOT macro is part of BAT and can only be run if BAT was
// installed correctly. The macro shows an example of fitting an
// efficiency with a function defined by the user. The input are two
// histograms. In the fit the uncertainties of the ratio are
// considered to be binomial.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x efficiencyFitterExample.C
//
// or
//
//    root[1] .L efficiencyFitterExample.C
//    root[2] efficiencyFitterExample()
//
// or from the command line
//
//    $ root efficiencyFitterExample.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x efficiencyFitterExample.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro
// the #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#define COMPILER (!defined(__CINT__) && !defined(__CLING__))

#if defined(__MAKECINT__) || defined(__ROOTCLING__) || COMPILER

#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCEfficiencyFitter.h>

#endif

// ---------------------------------------------------------
void efficiencyFitterExample()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // set style
    BCAux::SetStyle();

    // -------------------------
    // define the fit function, which is also used in the generation of the data
    TF1 f1("f1", "0.5*(1 + TMath::Erf((x-[0])/[1]))", 0.0, 100.0);
    f1.SetParLimits(0, 35.0, 45.0);
    f1.SetParLimits(1, 10.0, 20.0);
    f1.SetParNames("mean", "sigma");
    // -------------------------

    // -------------------------
    // Create new histograms:
    TH1D h_denom("h_denom", ";x;N", 100, 0.0, 100.0);
    h_denom.SetStats(kFALSE);

    // cloning ensures identical binning
    TH1D h_numer(h_denom);
    h_numer.SetName("h_numer");

    // set drawing options
    h_denom.SetLineColor(kBlack);
    h_denom.SetFillColor(5);
    h_denom.SetFillStyle(1001);

    h_numer.SetLineColor(kRed);
    h_numer.SetLineWidth(2);
    // -------------------------

    // -------------------------
    // Fill histograms with data:
    // initialize random number generator with seed = 1234
    gRandom = new TRandom3(1234);

    // Fill h_denom randomly from Landau distribution
    // with mean = 30 and sigma = 8; simulate 1000 events
    for (unsigned i = 0; i < 1000 ; ++i)
        h_denom.Fill( gRandom->Landau(30., 8.) );

    // Fill h_numer, with probability of accenting entry from h_denom
    // following an error function with a smearing according to the
    // binomial distribution with parameters (mean, sigma) := (40, 15)
    f1.SetParameters(40, 15);
    for (int b = 1; b <= h_numer.GetNbinsX(); ++b) {
        double eff = f1.Integral(h_numer.GetXaxis()->GetBinLowEdge(b), h_numer.GetXaxis()->GetBinUpEdge(b)) / h_numer.GetXaxis()->GetBinWidth(b);
        h_numer.SetBinContent(b, gRandom->Binomial(int(h_denom.GetBinContent(b)), eff));
    }
    // -------------------------

    // create a new efficiency fitter
    BCEfficiencyFitter hef(h_denom, h_numer, f1);
    hef.SetRandomSeed(1346);

    // set options for evaluating the fit function
    hef.SetFlagIntegration(false);

    // set Metropolis as marginalization method
    hef.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set options for MCMC
    hef.SetPrecision(BCEngineMCMC::kQuick);

    // // perform fit
    hef.Fit();

    // print data and fit
    TCanvas c("c1", "");
    c.Divide(2, 1);

    c.cd(1);
    h_denom.Draw();
    h_numer.Draw("same");

    c.cd(2);

    hef.DrawFit("", true); // draw with a legend

    // print to file
    c.Print("fit.pdf");

    // print marginalized distributions
    hef.PrintAllMarginalized("distributions.pdf");
}
