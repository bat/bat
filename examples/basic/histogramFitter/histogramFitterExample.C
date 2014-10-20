//
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

// The data fitted is generated randomly as a signal peak (gaussian)
// plus a flat background using a function CreateHistogram(nbins, ns, nb, seed)
// The arguments are 'nbins' - number of bins of the histogram,
// 'ns' - number of signal events, 'nb' - number of background events,
// 'seed' - initial seed for the random number generator. The location
// and the width of the signal peak can be set up using the variables
// 'mean' and 'sigma' below.

TH1D* CreateHistogram(int nbins, int ns, int nb, int seed = 0);

const double mean  = 65.0;
const double sigma =  5.0;

//
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

   // create data
   TH1D* hist = CreateHistogram(20, 100, 100, 132);

   // define a fit function
   TF1* f1 = new TF1("f1", "[0] / sqrt(2.0 * 3.1416) / [2] * exp(-(x-[1])*(x-[1])/2./[2]/[2]) + [3]", 0., 100.);
   f1->SetParLimits(0,  0.0, 200.0);
   f1->SetParLimits(1, 55.0,  75.0);
   f1->SetParLimits(2,  0.1,  10.0);
   f1->SetParLimits(3,  0.0,   2.0);

   // create a new histogram fitter
   BCHistogramFitter* hf = new BCHistogramFitter(hist, f1);

   // set Metropolis as marginalization method
   hf->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

   // set precision
   hf->MCMCSetPrecision(BCEngineMCMC::kMedium);

   // integrate function over bin (true) or use linear interpolation
   hf->SetFlagIntegration(false);

   // set priors
   hf->SetPriorConstant(0);
   hf->SetPriorConstant(1);
   hf->SetPriorConstant(2);
   hf->SetPriorConstant(3);

   // perform fit
   hf->Fit();

   double pvalue, pvalueCorrected;
   std::vector<double> init (4);
   init[0] = mean; init[1] = sigma; init[2] = 150; init[3] = 1;
   hf->FindMode(init);
   hf->CalculatePValueFast(hf->GetBestFitParameters());
   pvalue = hf->GetPValue();
   pvalueCorrected = hf->GetPValueNDoF();

   cout << "Pvalue " << pvalue
        << ", corrected " << pvalueCorrected
        << endl;

   // print marginalized distributions
   hf->PrintAllMarginalized("distributions.pdf");

   // print data and fit
   TCanvas* c1 = new TCanvas("c1");
   hf->DrawFit("", true); // draw with a legend
   c1->Print("fit.pdf");

   return;
}

// ---------------------------------------------------------
TH1D* CreateHistogram(int nbins, int ns, int nb, int seed)
{
   // initialize random number generator
   gRandom = new TRandom3(seed);

   // create new histogram
   TH1D* hist = new TH1D("data", ";x;N", nbins, 0.0, 100.0);
   hist->SetStats(kFALSE);

   // fill signal
   for (int i = 0; i < ns; ++i)
      hist->Fill(gRandom->Gaus(mean, sigma));

   // fill background
   for (int i = 0; i < nb; ++i)
      hist->Fill(100.0 * gRandom->Uniform());

   // return the histogram
   return hist;
}

// ---------------------------------------------------------
