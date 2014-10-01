/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>

#include "BCMTFComparisonTool.h"

// ---------------------------------------------------------
BCMTFComparisonTool::BCMTFComparisonTool(const char * name)
{
   fName = name;
}

// ---------------------------------------------------------
BCMTFComparisonTool::~BCMTFComparisonTool()
{
   for (int i = 0; i < GetNContributions(); ++i) {
      if (fHistogramContainer.at(i))
         delete fHistogramContainer.at(i);
   }
}

// ---------------------------------------------------------
void BCMTFComparisonTool::AddContribution(const char * name, TH1D hist)
{
   // add name to container
   fNameContainer.push_back(name);

   // add histogram to container
   fHistogramContainer.push_back(new TH1D(hist));

   // add central value to container
   fCentralValueContainer.push_back(hist.GetMean());

   // add uncertainty to container
   fUncertaintyContainer.push_back(hist.GetRMS());
}

// ---------------------------------------------------------
void BCMTFComparisonTool::AddContribution(const char * name, double centralvalue, double uncertainty)
{
   // add name to container
   fNameContainer.push_back(name);

   // add 0 to container
   fHistogramContainer.push_back(0);

   // add central value to container
   fCentralValueContainer.push_back(centralvalue);

   // add uncertainty to container
   fUncertaintyContainer.push_back(uncertainty);
}

// ---------------------------------------------------------
void BCMTFComparisonTool::PrintHistograms(const char * filename)
{
   // get number of histograms
   int nhistograms = (int) fHistogramContainer.size();

   // create canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // loop over all histograms
   for (int i = 0; i < nhistograms; ++i) {
      // get histogram
      TH1D * hist = fHistogramContainer.at(i);

      // set color
      hist->SetLineColor(2+i);

      // draw histogram
      if (i == 0) {
         hist->Draw("HIST");
         std::cout << " here as well." << std::endl;
      }
      else {
         hist->Draw("SAMEHIST");
         std::cout << " here as well 2." << std::endl;
      }
   }

   // print canvas
   c1->Print(filename);

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFComparisonTool::DrawOverview()
{
   // get number of contributions
   int ncontributions = GetNContributions();

   // create graph
   TGraphAsymmErrors * graph_contributions = new TGraphAsymmErrors(ncontributions);
   graph_contributions->SetMarkerStyle(20);
   graph_contributions->SetMarkerSize(1);

   // coordinate system
   double xmin = 0.0;
   double xmax = 0.0;
   double xwidth = 0.0;
   double ymin = -0.5;
   double ymax = double(ncontributions)-0.5;

   // ---- fill the graph ---- //

   // loop over all contributions
   for (int i = 0; i < ncontributions; ++i) {

      // get summary information
      double centralvalue = fCentralValueContainer.at(i);
      double uncertainty  = fUncertaintyContainer.at(i);

      // update coordinate system
      if ((centralvalue-uncertainty) < xmin || i == 0)
         xmin = centralvalue-uncertainty;
      if ((centralvalue+uncertainty) > xmax || i == 0)
         xmax = centralvalue+uncertainty;
      xwidth = xmax - xmin;

      // set point and error
      graph_contributions->SetPoint(i, centralvalue, double(ncontributions-i-1));
      graph_contributions->SetPointError(i, uncertainty, uncertainty, 0, 0);
   }

   // ---- do the plotting ---- //

   // create histogram for axes
   TH2D * hist_axes = new TH2D("", Form(";%s;", GetName().c_str()), 1, xmin - 0.25 *xwidth, xmax + 1.75 * xwidth, ncontributions, ymin, ymax);
   hist_axes->SetStats(kFALSE);
   hist_axes->GetYaxis()->SetNdivisions(0);
   hist_axes->GetYaxis()->SetTitleOffset(1.0);

   // create latex
   TLatex * latex = new TLatex();
   latex->SetTextSize(0.04);
   if (ncontributions>=10)
      latex->SetTextSize(0.02);
   latex->SetTextAlign(12);

   // draw
   hist_axes->Draw();
   graph_contributions->Draw("SAMEPZ");

   // loop over all contributions and draw labels
   for (int i = 0; i < ncontributions; ++i) {
      latex->DrawLatex(xmax + 0.25 *xwidth, double(ncontributions-i-1), fNameContainer.at(i).c_str());
   }
   hist_axes->Draw("SAMEAXIS");

   // todo: delete objects properly
}

// ---------------------------------------------------------
void BCMTFComparisonTool::PrintOverview(const char * filename)
{
   // create canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw the overview
   DrawOverview();

   // print to file
   c1->Print(filename);
}

// ---------------------------------------------------------
