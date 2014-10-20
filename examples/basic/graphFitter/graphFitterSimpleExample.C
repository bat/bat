//
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

// The data fitted are generated randomly to follow a straight line
// in function CreateGraph(n,seed), where 'n' is the number of points
// to generate and 'seed' is the random seed for the random number generator.
// The parameters of the true function can be set via variables 'slope'
// and 'offset'. The sigma of the gaussian for smearing can be defined via
// the variable 'sigma'. It is also set as the y-error of the data.

TGraphErrors * CreateGraph(int n, int seed = 0);

const double slope  =  1.0;
const double offset = 11.0;
const double sigma  =  5.0;

// The macro performs a straight line fit to the data. The fit function
// is defined using the ROOT TF1 object and the data to fit are stored
// in the TGraphErrors object. The y-errors have to be defined in the
// TGraphErrors for the fit to work.

// ---------------------------------------------------------
void graphFitterSimpleExample()
{
   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // create data
   TGraphErrors * graph = CreateGraph(10, 1000);

   // define a fit function
   TF1 * f1 = new TF1("f1", "[0] + [1]*x", 0.0, 100.0);

   // allowed range has to be defined for every parameter
   f1->SetParLimits(0, -20.0, 30.0);
   f1->SetParLimits(1,   0.5,  1.5);

   // create a new graph fitter
   BCGraphFitter * gf = new BCGraphFitter(graph, f1);

   // set Metropolis as marginalization method
   gf->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

   // set precision
   gf->MCMCSetPrecision(BCEngineMCMC::kMedium);

   // perform the fit
   gf->Fit();

   // print data and the fit
   TCanvas * c1 = new TCanvas("c1");
   gf->DrawFit("", true);
   c1->Print("fit.pdf");

   // print marginalized distributions
   gf->PrintAllMarginalized("distributions.pdf");

   // print results
   gf->PrintResults("results.txt");
}

// ---------------------------------------------------------
TGraphErrors * CreateGraph(int n, int seed)
{
   // initialize random number generator
   gRandom = new TRandom3(seed);

   // define arrays
   double * x  = new double[n];
   double * y  = new double[n];
   double * ey = new double[n];

   // define x and y-values and the uncertainties on y
   for (int i = 0; i < n; ++i)
   {
      x[i] = 100. / double(n) * double(i+0.5);
      y[i] = gRandom->Gaus(offset + slope*x[i], sigma);
      ey[i] = sigma;
   }

   // create new graph
   TGraphErrors * graph = new TGraphErrors(n, x, y, 0, ey);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);

   // return the graph
   return graph;
}

// ---------------------------------------------------------
