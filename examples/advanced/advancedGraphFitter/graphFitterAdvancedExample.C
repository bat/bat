//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly. The macro shows an example of fitting
// a graph using a function defined by the user. TGraphErrors
// with defined y uncertainties has to be defined. In the fit the
// uncertainties are considered to be Gaussian.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x graphFitterAdvancedExample.C
//
// or
//
//    root[1] .L graphFitterAdvancedExample.C
//    root[2] graphFitterAdvancedExample()
//
// or from the command line
//
//    $ root graphFitterAdvancedExample.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x graphFitterAdvancedExample.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro.
// The #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#define COMPILER (!defined(__CINT__) && !defined(__CLING__))

#if defined(__MAKECINT__) || defined(__ROOTCLING__) || COMPILER

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCGraphFitter.h>
#include <BAT/BCSummaryTool.h>

#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TString.h>

#include <fstream>
#include <vector>

#endif

// The includes below need to be always present

// By default the data in the file data/datax.txt are fitted. These
// are the same data that are used for the example in the BAT paper
// (see BAT webpage). To generate new data and fit them use
// CreateDataGraph() function.
//

TGraphErrors * CreateDataGraph(int n=100, double xmin=0.1, double xmax=19.9);

//
// The data are generated according to second order polynomial plus
// a gaussian peak. The parameters and the smearing can be set below

const double p0 =  0.;
const double p1 =  0.5;
const double p2 =  0.02;
const double a  = 15.;
const double m  =  5.;
const double s  =  0.5;

const double sigmay = 4.;

//
// The macro performs fits with four different functions (models):
//   2nd order polynomial
//   gaussian peak + constant
//   gaussian peak + straight line
//   gaussian peak + 2nd order polynomial

// ----------------------------------------------------------------------
void graphFitterAdvancedExample()
{
   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // uncomment the next command if you want the data
   // to be generated randomly
   //	TGraphErrors* gr = CreateDataGraph();

   // uncomment the next command if you want the data
   // to be read data from file
   // the file supplied with BAT distribution
   // is data/datax.txt and contains the data which were
   // used in the BAT paper
   TGraphErrors* gr = new TGraphErrors("data/datax.txt","%lg %lg %lg");

   gr->SetMarkerStyle(20);
   gr->SetMarkerSize(.5);

   // prepare fitting functions
   // 2nd order polynomial
   TF1* f1 = new TF1("f1", "[0]+[1]*x+[2]*x*x", 0., 100.);
   f1->SetParLimits(0,   0.,    5.);  // offset
   f1->SetParLimits(1,   0.,    1.2); // slope
   f1->SetParLimits(2,  -0.1,   0.1); // quad

   // constant + gaussian
   TF1* f2 = new TF1("f2", "[0]+[1]/(sqrt(2.*3.141592)*[3]) * exp(-(x-[2])*(x-[2])/(2.*[3]*[3]))", 0., 100.);
   f2->SetParLimits(0,   0.,   10.); // offset
   f2->SetParLimits(1,   0.,  200.); // A_gauss
   f2->SetParLimits(2,   2.,   18.); // mean
   f2->SetParLimits(3,   .2,    4.); // sigma

   // straight line + gaussian
   TF1* f3 = new TF1("f3", "[0]+[1]*x+[2]/(sqrt(2.*3.141592)*[4]) * exp(-(x-[3])*(x-[3])/(2.*[4]*[4]))", 0., 100.);
   f3->SetParLimits(0,   0.,   10.); // offset
   f3->SetParLimits(1,   0.,    2.); // slope
   f3->SetParLimits(2,   0.,  200.); // A_gauss
   f3->SetParLimits(3,   2.,   18.); // mean
   f3->SetParLimits(4,   .2,    4.); // sigma

   // 2nd order polynomial + gaussian
   TF1* f4 = new TF1("f4", "[0]+[1]*x+[2]*x*x+[3]/(sqrt(2.*3.141592)*[5]) * exp(-(x-[4])*(x-[4])/(2.*[5]*[5]))", 0., 100.);
   f4->SetParLimits(0,   0.,   10.); // offset
   f4->SetParLimits(1,   0.,    2.); // slope
   f4->SetParLimits(2,   0.,    .5); // quad
   f4->SetParLimits(3,   0.,  200.); // A_gauss
   f4->SetParLimits(4,   2.,   18.); // mean
   f4->SetParLimits(5,   .2,    4.); // sigma

   // setup all graph fitters
   vector<BCGraphFitter *> models;
   models.push_back(new BCGraphFitter(gr,f1));
   models.push_back(new BCGraphFitter(gr,f2));
   models.push_back(new BCGraphFitter(gr,f3));
   models.push_back(new BCGraphFitter(gr,f4));

   // perform the analysis on all models
   for (unsigned int imodel=0; imodel<models.size(); ++imodel) {

      // set precision
      models[imodel]->MCMCSetPrecision(BCEngineMCMC::kMedium);

      // run the fitting (MCMC + Minuit)
      models[imodel]->Fit();

      // print all marginalized distributions into a PostScript file
      models[imodel]->PrintAllMarginalized(TString::Format("plots-%d.pdf",imodel+1));

      // draw summary plots and tables
      BCSummaryTool summary(models[imodel]);
      summary.PrintParameterPlot(TString::Format("summary_pars-f%d.pdf",imodel+1));
      summary.PrintCorrelationPlot(TString::Format("summary_corr-f%d.png",imodel+1));
      summary.PrintCorrelationMatrix(TString::Format("summary_corr_matrix-f%d.pdf",imodel+1));
      //summary.PrintKnowledgeUpdatePlots(TString::Format("summary_update-f%d.pdf",imodel+1));
      //summary.PrintParameterLatex(TString::Format("summary_pars-f%d.tex",imodel+1));
   }

   // draw fit including the error band
   TCanvas* c = new TCanvas();
   c->Divide(2,2);
   for (unsigned int imodel=0; imodel<models.size(); ++imodel) {
      c->cd(imodel+1);
      models[imodel]->DrawFit();
   }
   c->Print("data-all-band.pdf");
   delete c;

   // draw all fits in the same plot (w/o error bands)
   c = new TCanvas();
   gr->Draw("ap");
   f1->SetLineColor(1);
   f1->SetLineWidth(2);
   f1->Draw("l same");
   f2->SetLineColor(2);
   f2->SetLineWidth(2);
   f2->Draw("l same");
   f3->SetLineColor(3);
   f3->SetLineWidth(2);
   f3->Draw("l same");
   f4->SetLineColor(4);
   f4->SetLineWidth(2);
   f4->Draw("l same");
   c->Print("data-all.pdf");

}


// ----------------------------------------------------------------------
TGraphErrors* CreateDataGraph(int n, double xmin, double xmax)
{
   const double pi =  3.141592;

   // initialize random number generator
   TRandom3* ran = new TRandom3(0);

   double* xx = new double[n];
   double* yy = new double[n];
   double* err = new double[n];

   double dx = (xmax-xmin)/(double)(n-1);

   // loop over points
   for (int i=0;i<n;i++)
   {
      // get x value
      double x = xmin + (double)i*dx;

      // get y value
      double yexp = p0 + p1*x + p2*x*x + a/(sqrt(2.*pi)*s)*exp(-.5*(x-m)*(x-m)/(s*s));
      double y = ran->Gaus(yexp*100., sigmay*100.) / 100.;

      xx[i]=x;
      yy[i]=y;
      err[i]=sigmay;
   }

   TGraphErrors* g = new TGraphErrors(n,xx,yy,0,err);

   delete [] xx;
   delete [] yy;
   delete [] err;

   return g;
}

// ----------------------------------------------------------------------
