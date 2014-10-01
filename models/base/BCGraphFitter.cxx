/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <TPad.h>
#include <TLegend.h>
#include <Math/ProbFuncMathCore.h>

#include "../../BAT/BCLog.h"
#include "../../BAT/BCDataSet.h"
#include "../../BAT/BCDataPoint.h"
#include "../../BAT/BCMath.h"

#include "BCGraphFitter.h"

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter()
 : BCFitter()
 , fGraph(0)
 , fFitFunction(0)
 , fErrorBand(0)
 , fGraphFitFunction(0)
{
   // set MCMC for marginalization
   SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter(const char * name)
 : BCFitter(name)
 , fGraph(0)
 , fFitFunction(0)
 , fErrorBand(0)
 , fGraphFitFunction(0)
{
   // set MCMC for marginalization
   SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter(TGraphErrors * graph, TF1 * func)
 : BCFitter()
 , fGraph(0)
 , fFitFunction(0)
 , fErrorBand(0)
 , fGraphFitFunction(0)
{
   SetGraph(graph);
   SetFitFunction(func);

   // set MCMC for marginalization
   SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter(const char * name, TGraphErrors * graph, TF1 * func)
 : BCFitter(name)
 , fGraph(0)
 , fFitFunction(0)
 , fErrorBand(0)
 , fGraphFitFunction(0)
{
   SetGraph(graph);
   SetFitFunction(func);

   // set MCMC for marginalization
   SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
}

// ---------------------------------------------------------

int BCGraphFitter::SetGraph(TGraphErrors * graph)
{
   if(!graph) {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors not created.");
      return 0;
   }

   int npoints = graph->GetN();
   if(!npoints)
   {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors is empty.");
      return 0;
   }
   else if(npoints==1)
   {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors has only one point. Not able to fit.");
      return 0;
   }

   fGraph = graph;

   double * x = fGraph->GetX();
   double * y = fGraph->GetY();
   double * ex = fGraph->GetEX();
   double * ey = fGraph->GetEY();

   if(!ey)
   {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors has NO errors set on Y. Not able to fit.");
      return 0;
   }

   BCDataSet * ds = new BCDataSet();

   // fill the dataset
   // find x and y boundaries for the error band calculation
   double xmin=x[0];
   double xmax=x[0];
   double ymin=y[0];
   double ymax=y[0];
   for (int i = 0; i < npoints; ++i)
   {
      // if x errors are not set, set them to zero
      double errx = ex ? ex[i] : 0.;

      // create the data point
      BCDataPoint * dp = new BCDataPoint(4);
      dp->SetValue(0, x[i]);
      dp->SetValue(1, y[i]);
      dp->SetValue(2, errx);
      dp->SetValue(3, ey[i]);
      ds->AddDataPoint(dp);

      if(x[i]-errx < xmin)
         xmin = x[i]-errx;
      else if(x[i]+errx > xmax)
         xmax = x[i]+errx;

      if(y[i] - 5.*ey[i] < ymin)
         ymin = y[i] - 5.*ey[i];
      else if(y[i] + 5.*ey[i] > ymax)
         ymax = y[i] + 5.*ey[i];
   }

   SetDataSet(ds);

   // set boundaries for the error band calculation
   SetDataBoundaries(0, xmin, xmax);
   SetDataBoundaries(1, ymin, ymax);

   SetFitFunctionIndices(0, 1);

   return GetNDataPoints();
}

// ---------------------------------------------------------

int BCGraphFitter::SetFitFunction(TF1 * func)
{
   if(!func)
   {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetFitFunction() : TF1 not created.");
      return 0;
   }

   // get the new number of parameters
   int npar = func->GetNpar();
   if(!npar)
   {
      BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetFitFunction() : TF1 has zero parameters. Not able to fit.");
      return 0;
   }

   // set the function
   fFitFunction = func;

   // update the model name to contain the function name
   if(fName=="model")
      SetName(TString::Format("GraphFitter with %s",fFitFunction->GetName()));

   // reset parameters
   ClearParameters(true);

   // add parameters
   for (int i = 0; i < npar; ++i)
   {
      double xmin;
      double xmax;
      fFitFunction->GetParLimits(i, xmin, xmax);

      AddParameter(fFitFunction->GetParName(i), xmin, xmax);
   }

   // set flat prior for all parameters by default
   SetPriorConstantAll();

   return GetNParameters();
}

// ---------------------------------------------------------

BCGraphFitter::~BCGraphFitter()
{
   // todo memory leak
}

// ---------------------------------------------------------

double BCGraphFitter::LogLikelihood(const std::vector<double> & params)
{
   // initialize probability
   double logl = 0.;

   // set the parameters of the function
   // passing the pointer to first element of the vector is
   // not completely safe as there might be an implementation where
   // the vector elements are not stored consecutively in memory.
   // however it is much faster than copying the contents, especially
   // for large number of parameters
   fFitFunction->SetParameters(&params[0]);

   // loop over all data points
   for (unsigned i = 0; i < GetNDataPoints(); i++)
   {
      std::vector<double> x = GetDataPoint(i)->GetValues();

      // her we ignore the errors on x even when they're available
      // i.e. we treat them just as the region specifiers
      double y = x[1];
      double yerr = x[3];
      double yexp = fFitFunction->Eval(x[0]);

      // calculate log of probability assuming
      // a Gaussian distribution for each point
      logl += BCMath::LogGaus(y, yexp, yerr, true);
   }

   return logl;
}

// ---------------------------------------------------------

double BCGraphFitter::FitFunction(const std::vector<double> & x, const std::vector<double> & params)
{
   // set the parameters of the function
   // passing the pointer to first element of the vector is
   // not completely safe as there might be an implementation where
   // the vector elements are not stored consecutively in memory.
   // however it is much faster than copying the contents, especially
   // for large number of parameters
   fFitFunction->SetParameters(&params[0]);

   return fFitFunction->Eval(x[0]);
}

// ---------------------------------------------------------

int BCGraphFitter::Fit(TGraphErrors * graph, TF1 * func)
{
   // set graph
   if (!SetGraph(graph)) {
      BCLog::OutError("BCEfficiencyFitter::Fit : Graph not defined.");
      return 0;
   }

   // set function
   if (!SetFitFunction(func)) {
      return 0;
      BCLog::OutError("BCEfficiencyFitter::Fit : Fit function not defined.");
   }

   return Fit();
}

// ---------------------------------------------------------

int BCGraphFitter::Fit()
{
   // set graph
   if (!fGraph) {
      BCLog::OutError("BCEfficiencyFitter::Fit : Graph not defined.");
      return 0;
   }

   // set function
   if (!fFitFunction) {
      BCLog::OutError("BCEfficiencyFitter::Fit : Fit function not defined.");
      return 0;
   }

   // check setup
   BCLog::Out(BCLog::detail,BCLog::detail,
         Form("Fitting %d data points with function of %d parameters",GetNDataPoints(),GetNParameters()));
   if(GetNDataPoints() <= GetNParameters())
   {
      BCLog::Out(BCLog::warning,BCLog::warning,
            Form("Number of parameters (%d) lower than or equal to number of points (%d)."
            , GetNParameters(), GetNDataPoints()));
      BCLog::Out(BCLog::warning,BCLog::warning,"Fit doesn't have much meaning.");
   }

   // perform marginalization
   MarginalizeAll();

   // maximize posterior probability, using the best-fit values close
   // to the global maximum from the MCMC
   BCIntegrate::BCOptimizationMethod method_temp = GetOptimizationMethod();
   SetOptimizationMethod(BCIntegrate::kOptMinuit);
   FindMode( GetBestFitParameters());
   SetOptimizationMethod(method_temp);

   // calculate p-value from the chi2 probability
   // this is only valid for a product of gaussiang which is the case for
   // the BCGraphFitter
   GetPvalueFromChi2(GetBestFitParameters(), 3);
   GetPvalueFromChi2NDoF(GetBestFitParameters(), 3);

   // print summary to screen
   PrintShortFitSummary(1);

   return 1;
}

// ---------------------------------------------------------

void BCGraphFitter::DrawFit(const char * options, bool flaglegend)
{
   if (!fGraph)
   {
      BCLog::Out(BCLog::error, BCLog::error,"BCGraphFitter::DrawFit() : TGraphErrors not defined.");
      return;
   }

   if (!fFitFunction)
   {
      BCLog::Out(BCLog::error, BCLog::error,"BCGraphFitter::DrawFit() : Fit function not defined.");
      return;
   }

   // check wheather options contain "same"
   TString opt = options;
   opt.ToLower();

   // if not same, draw the histogram first to get the axes
   if(!opt.Contains("same"))
      fGraph->Draw("ap");

   // draw the error band as central 68% probability interval
   fErrorBand = GetErrorBandGraph(0.16, 0.84);
   fErrorBand->Draw("f same");

   // draw the fit function on top
   fGraphFitFunction = GetFitFunctionGraph( GetBestFitParameters() );
   fGraphFitFunction->SetLineColor(kRed);
   fGraphFitFunction->SetLineWidth(2);
   fGraphFitFunction->Draw("l same");

   // now draw the histogram again since it was covered by the band and
   // the best fit
   fGraph->Draw("p same");

   // draw legend
   if (flaglegend)
   {
      TLegend * legend = new TLegend(0.25, 0.75, 0.55, 0.95);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->AddEntry(fGraph, "Data", "P");
      legend->AddEntry(fGraphFitFunction, "Best fit", "L");
      legend->AddEntry(fErrorBand, "Error band", "F");
      legend->Draw();
   }

   gPad->RedrawAxis();
}

// ---------------------------------------------------------

double BCGraphFitter::CDF(const std::vector<double> & parameters,  int index, bool /*lower*/) {

   //format: x y error_x error_y
   std::vector<double> values = fDataSet->GetDataPoint(index)->GetValues();

   if (values.at(2))
      BCLog::OutWarning("BCGraphFitter::CDF: Non-zero errors in x-direction are ignored!");

   // get the observed value
   double yObs = values.at(1);

   // expectation value
   double yExp = FitFunction(values, parameters);

   return ROOT::Math::normal_cdf(yObs, values.at(3), yExp);
}

// ---------------------------------------------------------
