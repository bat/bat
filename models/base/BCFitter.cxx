/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCFitter.h"

#include "BAT/BCModel.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCDataSet.h"
#include "BAT/BCLog.h"
#include "BAT/BCH1D.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

// ---------------------------------------------------------
BCFitter::BCFitter() : BCModel()
                     , fErrorBand(0)
                     , fFlagFillErrorBand(true)
                     , fFitFunctionIndexX(-1)
                     , fFitFunctionIndexY(-1)
                     , fErrorBandContinuous(true)
                     , fErrorBandX(std::vector<double>(0))
                     , fErrorBandNbinsX(100)
                     , fErrorBandNbinsY(500)
                     , fErrorBandXY(0)
{
}

// ---------------------------------------------------------
BCFitter::BCFitter(const char * name) : BCModel(name)
                                      , fErrorBand(0)
                                      , fFlagFillErrorBand(true)
                                      , fFitFunctionIndexX(-1)
                                      , fFitFunctionIndexY(-1)
                                      , fErrorBandContinuous(true)
                                      , fErrorBandX(std::vector<double>(0))
                                      , fErrorBandNbinsX(100)
                                      , fErrorBandNbinsY(500)
                                      , fErrorBandXY(0)
{
}

// ---------------------------------------------------------
BCFitter::~BCFitter()
{
}

// ---------------------------------------------------------
void BCFitter::MCMCIterationInterface()
{
  // call base interface
  BCEngineMCMC::MCMCIterationInterface();

  // fill error band
  if (fFlagFillErrorBand)
    FillErrorBand();
}

// ---------------------------------------------------------
void BCFitter::MarginalizePreprocess()
{
  // prepare function fitting
  double dx = 0.;
  double dy = 0.;

  if (fFitFunctionIndexX >= 0) {
    dx = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX)
          - fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX))
      / (double) fErrorBandNbinsX;

    dy = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY)
          - fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY))
      / (double) fErrorBandNbinsY;

    fErrorBandXY
      = new TH2D(TString::Format("errorbandxy_%d", BCLog::GetHIndex()), "",
                 fErrorBandNbinsX,
                 fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX) - .5 * dx,
                 fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX) + .5 * dx,
                 fErrorBandNbinsY,
                 fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY) - .5 * dy,
                 fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY) + .5 * dy);
    fErrorBandXY->SetStats(kFALSE);

    for (unsigned ix = 1; ix <= fErrorBandNbinsX; ++ix)
      for (unsigned iy = 1; iy <= fErrorBandNbinsX; ++iy)
        fErrorBandXY->SetBinContent(ix, iy, 0.);
  }

}

// ---------------------------------------------------------
void BCFitter::FillErrorBand()
{
  // function fitting
  if (fFitFunctionIndexX < 0)
    return;

  // loop over all possible x values ...
  if (fErrorBandContinuous) {
    double x = 0;
    for (unsigned ix = 0; ix < fErrorBandNbinsX; ix++) {
      // calculate x
      x = fErrorBandXY->GetXaxis()->GetBinCenter(ix + 1);

      // calculate y
      std::vector<double> xvec;
      xvec.push_back(x);

      // loop over all chains
      for (unsigned ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
        // calculate y
        double y = FitFunction(xvec, MCMCGetx(ichain));

        // fill histogram
        fErrorBandXY->Fill(x, y);
      }

      xvec.clear();
    }
  }
  // ... or evaluate at the data point x-values
  else {
    unsigned ndatapoints = fErrorBandX.size();
    double x = 0;

    for (unsigned ix = 0; ix < ndatapoints; ++ix) {
      // calculate x
      x = fErrorBandX.at(ix);

      // calculate y
      std::vector<double> xvec;
      xvec.push_back(x);

      // loop over all chains
      for (unsigned ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
        // calculate y
        double y = FitFunction(xvec, MCMCGetx(ichain));

        // fill histogram
        fErrorBandXY->Fill(x, y);
      }

      xvec.clear();
    }
  }
}

// ---------------------------------------------------------
std::vector<double> BCFitter::GetErrorBand(double level) const
{
   std::vector<double> errorband;

   if (!fErrorBandXY)
      return errorband;

   int nx = fErrorBandXY->GetNbinsX();
   errorband.assign(nx, 0.);

   // loop over x and y bins
   for (int ix = 1; ix <= nx; ix++) {
      TH1D * temphist = fErrorBandXY->ProjectionY("temphist", ix, ix);

      int nprobSum = 1;
      double q[1];
      double probSum[1];
      probSum[0] = level;

      temphist->GetQuantiles(nprobSum, q, probSum);

      errorband[ix - 1] = q[0];
   }

   return errorband;
}

// ---------------------------------------------------------
TGraph * BCFitter::GetErrorBandGraph(double level1, double level2) const
{
   if (!fErrorBandXY)
      return 0;

   // define new graph
   int nx = fErrorBandXY->GetNbinsX();

   TGraph * graph = new TGraph(2 * nx);
   graph->SetFillStyle(1001);
   graph->SetFillColor(kYellow);

   // get error bands
   std::vector<double> ymin = GetErrorBand(level1);
   std::vector<double> ymax = GetErrorBand(level2);

   for (int i = 0; i < nx; i++) {
      graph->SetPoint(i, fErrorBandXY->GetXaxis()->GetBinCenter(i + 1), ymin[i]);
      graph->SetPoint(nx + i, fErrorBandXY->GetXaxis()->GetBinCenter(nx - i), ymax[nx - i - 1]);
   }

   return graph;
}

// ---------------------------------------------------------
TH2D * BCFitter::GetErrorBandXY_yellow(double level, int nsmooth) const
{
   if (!fErrorBandXY)
      return 0;

   int nx = fErrorBandXY->GetNbinsX();
   int ny = fErrorBandXY->GetNbinsY();

   // copy existing histogram
   TH2D * hist_tempxy = (TH2D*) fErrorBandXY->Clone(
         TString::Format("%s_sub_%f.2", fErrorBandXY->GetName(), level));
   hist_tempxy->Reset();
   hist_tempxy->SetFillColor(kYellow);

   // loop over x bins
   for (int ix = 1; ix < nx; ix++) {
      BCH1D * hist_temp = new BCH1D();

      TH1D * hproj = fErrorBandXY->ProjectionY("temphist", ix, ix);
      if (nsmooth > 0)
         hproj->Smooth(nsmooth);

      hist_temp->SetHistogram(hproj);

      TH1D * hist_temp_yellow = hist_temp->GetSmallestIntervalHistogram(level);

      for (int iy = 1; iy <= ny; ++iy)
         hist_tempxy->SetBinContent(ix, iy, hist_temp_yellow->GetBinContent(iy));

      delete hist_temp_yellow;
      delete hist_temp;
   }

   return hist_tempxy;
}

// ---------------------------------------------------------
TGraph * BCFitter::GetFitFunctionGraph(const std::vector<double> &parameters)
{
   if (!fErrorBandXY)
      return 0;

   // define new graph
   int nx = fErrorBandXY->GetNbinsX();
   TGraph * graph = new TGraph(nx);

   // loop over x values
   for (int i = 0; i < nx; i++) {
      double x = fErrorBandXY->GetXaxis()->GetBinCenter(i + 1);

      std::vector<double> xvec;
      xvec.push_back(x);
      double y = FitFunction(xvec, parameters);
      xvec.clear();

      graph->SetPoint(i, x, y);
   }

   return graph;
}

// ---------------------------------------------------------
TGraph * BCFitter::GetFitFunctionGraph(const std::vector<double> &parameters, double xmin, double xmax, int n)
{
   // define new graph
   TGraph * graph = new TGraph(n + 1);

   double dx = (xmax - xmin) / (double) n;

   // loop over x values
   for (int i = 0; i <= n; i++) {
      double x = (double) i * dx + xmin;
      std::vector<double> xvec;
      xvec.push_back(x);
      double y = FitFunction(xvec, parameters);

      xvec.clear();

      graph->SetPoint(i, x, y);
   }

   return graph;
}

// ---------------------------------------------------------
int BCFitter::ReadErrorBandFromFile(const char * file)
{
   TFile * froot = TFile::Open(file);
   if (!froot->IsOpen()) {
      BCLog::OutError(Form("BCFitter::ReadErrorBandFromFile. Couldn't open file %s.", file));
      return 0;
   }

   int r = 0;

   TH2D * h2 = (TH2D*) froot->Get("errorbandxy");
   if (h2) {
      h2->SetDirectory(0);
      h2->SetName(TString::Format("errorbandxy_%d", BCLog::GetHIndex()));
      SetErrorBandHisto(h2);
      r = 1;
   }
   else
      BCLog::OutWarning(
            Form("BCFitter::ReadErrorBandFromFile : Couldn't read histogram \"errorbandxy\" from file %s.",file));

   froot->Close();

   return r;
}

// ---------------------------------------------------------
void BCFitter::FixDataAxis(unsigned int index, bool fixed)
{
   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCFitter::FixDataAxis : Index out of range.");
      return;
   }

   if (fDataFixedValues.size() == 0)
      fDataFixedValues.assign(fDataSet->GetDataPoint(0)->GetNValues(),
            false);

   fDataFixedValues[index] = fixed;
}

// ---------------------------------------------------------
bool BCFitter::GetFixedDataAxis(unsigned int index) const
{
   // check if index is within range
   if (index > fDataSet->GetDataPoint(0)->GetNValues()) {
      BCLog::OutError("BCFitter::GetFixedDataAxis : Index out of range.");
      return false;
   }

   return fDataFixedValues[index];
}

// ---------------------------------------------------------
void BCFitter::SetErrorBandContinuous(bool flag)
{
   fErrorBandContinuous = flag;

   if (flag)
      return;

   // clear x-values
   fErrorBandX.clear();

   // copy data x-values
   for (unsigned int i = 0; i < fDataSet->GetNDataPoints(); ++i)
      fErrorBandX.push_back(fDataSet->GetDataPoint(i)->GetValue(fFitFunctionIndexX));
}

// ---------------------------------------------------------
