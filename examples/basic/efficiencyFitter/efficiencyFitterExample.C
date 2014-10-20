//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// an efficiency with a function defined by the user. The input are
// two histograms, one being a subset of the other. In the fit the
// uncertainties of the ratio are considered to be binomial.
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

// The data fitted are generated randomly using a function
// CreateHistograms(nbins,nevents,hist1,hist2,seed), where 'nbins' is
// the number of bins in both histograms, nevents is the number of
// events generated, hist1 and hist2 are pointers to the generated TH1D
// objects containing the generated histograms with hist1 being the
// denominator and hist2 being the numerator (subset of hist1).
// The denominator histogram hist1 is generated from Landau distribution.
// The numerator histogram hist2 is created as a subset of hist1 where
// the probability of accepting the entry in both histograms is set to
// follow an error function with a smearing according to the binomial
// distribution. The parameters of the error function can be adjusted
// using the variables 'mean' and 'sigma'

void CreateHistograms(int nbins, int nevents, TH1D * hist1, TH1D * hist2, int seed = 0);

const double mean  = 40.0;
const double sigma = 15.0;

//
// The Data are fitted using an error function defined in fitfunction()
// and passed to TF1 object.

double fitfunction(double * x, double * par);

// ---------------------------------------------------------
void efficiencyFitterExample()
{
   int nbins = 100;

   // open log file
   BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

   // set style
   BCAux::SetStyle();

   // create new histograms
   TH1D * hist1 = new TH1D("data1", ";x;N", nbins, 0.0, 100.0);
   hist1->SetStats(kFALSE);
   hist1->SetLineColor(kBlack);
   hist1->SetFillColor(5);
   hist1->SetFillStyle(1001);

   TH1D * hist2 = new TH1D("data2", ";x;N", nbins, 0.0, 100.0);
   hist2->SetStats(kFALSE);
   hist2->SetLineColor(kRed);

   // create data
   CreateHistograms(nbins, 1000, hist1, hist2, 1000);

   // define a fit function
   TF1 * f1 = new TF1("f1", fitfunction, 0.0, 100.0, 2);
   f1->SetParLimits(0, 36.0, 44.0);
   f1->SetParLimits(1, 10.0, 20.0);

   // create a new efficiency fitter
   BCEfficiencyFitter * hef = new BCEfficiencyFitter(hist1, hist2, f1);

   // set options for evaluating the fit function
   hef->SetFlagIntegration(false);

   // set Metropolis as marginalization method
   hef->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

   // set options for MCMC
   hef->MCMCSetPrecision(BCEngineMCMC::kMedium);

   // set priors
   hef->SetPriorConstant(0);
   hef->SetPriorConstant(1);

   // perform fit
   hef->Fit();

   // print data and fit
   TCanvas * c = new TCanvas("c1","",800,400);
   c->Divide(2, 1);
   c->cd(1);
   hist1->Draw();
   hist2->Draw("same");
   c->cd(2);
   hef->DrawFit("", true); // draw with a legend

   // print to file
   c->Print("fit.pdf");

   // print marginalized distributions
   hef->PrintAllMarginalized("distributions.pdf");
}

// ---------------------------------------------------------
double fitfunction(double * x, double * par)
{
   double ff;

   if (x[0] < par[0])
      ff = TMath::Erfc((par[0]-x[0])/par[1])/2.;
   else
      ff = TMath::Erf((x[0]-par[0])/par[1])/2. + 0.5;

   return ff;
}

// ---------------------------------------------------------
void CreateHistograms(int nbins, int nevents, TH1D * hist1, TH1D * hist2, int seed)
{
   // initialize random number generator
   gRandom = new TRandom3(seed);

   // fill histogram 1
   for (int i = 0; i < nevents ; ++i)
   {
      double x = gRandom->Landau(30., 8.);
      hist1->Fill(x);
   }

   // fill histogram 2
   for (int i = 1; i <= nbins; ++i)
   {
      double x = hist1->GetXaxis()->GetBinCenter(i);
      double eff = 0;
      if (x < mean)
         eff = TMath::Erfc((mean-x)/sigma)/2.;
      else
         eff = TMath::Erf((x-mean)/sigma)/2. + 0.5;
      //		int n = gRandom->Poisson(nevents);
      //		hist1->SetBinContent(i, n);
      //		hist2->SetBinContent(i, gRandom->Binomial(n, eff));
      hist2->SetBinContent(i, gRandom->Binomial(int(hist1->GetBinContent(i)), eff));
   }
}

// ---------------------------------------------------------
