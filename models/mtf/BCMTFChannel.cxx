/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <iostream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>

#include "BCMTFTemplate.h"
#include "BCMTFSystematicVariation.h"

#include "BCMTFChannel.h"

// ---------------------------------------------------------
BCMTFChannel::BCMTFChannel(const char * name)
   : fData(0)
   , fFlagChannelActive(true)
   , fHistUncertaintyBandExpectation(0)
   , fHistUncertaintyBandPoisson(0)
{
   fName = name;
}

// ---------------------------------------------------------
BCMTFChannel::~BCMTFChannel()
{
   if (fData)
      delete fData;

   for (unsigned int i = 0; i < fTemplateContainer.size(); ++i)
      delete (fTemplateContainer.at(i));

   for (unsigned int i = 0; i < fSystematicVariationContainer.size(); ++i)
      delete (fSystematicVariationContainer.at(i));

   /*
     if (fHistUncertaintyBandExpectation)
     delete fHistUncertaintyBandExpectation;

     if (fHistUncertaintyBandPoisson)
     delete fHistUncertaintyBandPoisson;
   */
}

// ---------------------------------------------------------
void BCMTFChannel::PrintTemplates(std::string filename)
{
   // check if file extension is pdf
   if ( (filename.find_last_of(".") != std::string::npos) &&
        (filename.substr(filename.find_last_of(".")+1) == "pdf") ) {
      ; // it's a PDF file

   }
   else if ( (filename.find_last_of(".") != std::string::npos) &&
             (filename.substr(filename.find_last_of(".")+1) == "ps") ) {
      ; // it's a PS file
   }
   else {
      ; // make it a PDF file
      filename += ".pdf";
   }

   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // get number of templates
   int ntemplates = int(fTemplateContainer.size());

   int first_hist = -1;
   int last_hist  = 0;

    // calculate first and last existing histogram
   for (int i = 0; i < ntemplates; ++i) {
      // get histogram
      TH1D * temphist = GetTemplate(i)->GetHistogram();

      if (first_hist < 0 && temphist)
         first_hist = i;

      if (temphist)
         last_hist = i;
   }

   // Don't print if there are no histograms
   if ( first_hist < 0 ) {
      // free memory
      delete c1;

      return;
   }

  // loop over templates
   for (int i = 0; i < ntemplates; ++i) {

      // get histogram
      TH1D * temphist = GetTemplate(i)->GetHistogram();

      TLatex* l = new TLatex();

      // draw
      if (temphist) {
         temphist->Draw();
         l->DrawTextNDC(0.2, 0.9, Form("%s - %s", fName.c_str(), GetTemplate(i)->GetProcessName().c_str()));
      }

      // print
      if (i == first_hist && (first_hist != last_hist))
         c1->Print(std::string( filename + "(").c_str());
      else if (i == last_hist && (first_hist != last_hist))
         c1->Print(std::string( filename + ")").c_str());
      else {
         if (temphist)
            c1->Print(filename.c_str());
      }

      // free memory
      delete l;
   }

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFChannel::PrintTemplate(int index, const char * filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // get number of systematics
   unsigned int nsystematics = fSystematicVariationContainer.size();

   // draw template
   GetTemplate(index)->GetHistogram()->Draw();

   // print to file
   if (nsystematics == 0)
      c1->Print(filename);
   else
      c1->Print( (std::string(filename)+std::string("(")).c_str() );

   // loop over systematics
   for (unsigned int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      c1->cd();

      // check that histogram exists
      if (!GetSystematicVariation(isystematic)->GetHistogramUp(index) ||
          !GetSystematicVariation(isystematic)->GetHistogramDown(index))
         continue;

      // get histogram
      TH1D hist = TH1D(*GetTemplate(index)->GetHistogram());
      TH1D hist_up(hist);
      TH1D hist_down(hist);

      // set style
      hist.SetFillStyle(0);
      hist_up.SetFillStyle(0);
      hist_up.SetLineStyle(2);
      hist_down.SetFillStyle(0);
      hist_down.SetLineStyle(3);

      // get efficiency
      double efficiency = GetTemplate(index)->GetEfficiency();

      // scale histogram
      hist.Scale(efficiency/ hist.Integral());

      // loop over all bins
      for (int i = 1; i <= hist.GetNbinsX(); ++i) {
         hist.SetBinContent(i, GetTemplate(index)->GetHistogram()->GetBinContent(i));
         hist_up.SetBinContent(i, hist.GetBinContent(i) * (1.0 + GetSystematicVariation(isystematic)->GetHistogramUp(index)->GetBinContent(i)));
         hist_down.SetBinContent(i, hist.GetBinContent(i) * (1.0 - GetSystematicVariation(isystematic)->GetHistogramDown(index)->GetBinContent(i)));
      }

      // draw histogram
      hist_up.Draw("HIST");
      hist.Draw("HISTSAME");
      hist_down.Draw("HISTSAME");

      if (isystematic < nsystematics-1)
         c1->Print(filename);
      else
         c1->Print( (std::string(filename)+std::string(")")).c_str() );
   }

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFChannel::PrintHistUncertaintyBandExpectation(const char* filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw histogram
   fHistUncertaintyBandExpectation->Draw("COLZ");
   c1->Draw();

   // print
   c1->Print(filename);

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFChannel::CalculateHistUncertaintyBandPoisson()
{
   // calculate histogram
   int nbinsy_exp = fHistUncertaintyBandExpectation->GetNbinsY();
   int nbinsx_poisson = fHistUncertaintyBandPoisson->GetNbinsX();
   int nbinsy_poisson = fHistUncertaintyBandPoisson->GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_w = 0;
      // loop over y-axis of expectation and calculate sum of weights
      for (int iy = 1; iy <= nbinsy_exp; ++iy) {
         double w = fHistUncertaintyBandExpectation->GetBinContent(ix, iy);
         sum_w += w;
      }
      // loop over y-axis of expectation
      for (int iy = 1; iy <= nbinsy_exp; ++iy) {
         double w = fHistUncertaintyBandExpectation->GetBinContent(ix, iy)/sum_w;
         double expectation = fHistUncertaintyBandExpectation->GetYaxis()->GetBinCenter(iy);
         // loop over y-axis of observation
         for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
            double p = TMath::Poisson(double(jbin-1), expectation);
            double bincontent = 0;
            if (iy>1)
               bincontent=fHistUncertaintyBandPoisson->GetBinContent(ix, jbin);
            fHistUncertaintyBandPoisson->SetBinContent(ix, jbin, bincontent+p*w);
         }
      }
   }

}

// ---------------------------------------------------------
TH1D* BCMTFChannel::CalculateUncertaintyBandPoisson(double minimum, double maximum, int color)
{
   TH1D* hist = new TH1D(*(fData->GetHistogram()));
   hist->SetMarkerSize(0);
   hist->SetFillColor(color);
   hist->SetFillStyle(1001);

   int nbinsx_poisson = fHistUncertaintyBandPoisson->GetNbinsX();
   int nbinsy_poisson = fHistUncertaintyBandPoisson->GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_p = 0;  // sum of all probabilities inside the interval
      int limit_min = 0;
      int limit_max = nbinsx_poisson-1;

      // loop over y-axis of observation
      for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
         double p = fHistUncertaintyBandPoisson->GetBinContent(ix, jbin);
         sum_p+=p;
         if (sum_p < minimum)
            limit_min=jbin;
         if (sum_p > maximum && (sum_p - p) < maximum )
            limit_max=jbin-1;
      }
      //    hist->SetBinContent(ix, 0.5*double(limit_min+limit_max));
      //    hist->SetBinError(ix, 0.5*double(limit_max-limit_min));
      double ylimit_min = fHistUncertaintyBandPoisson->GetYaxis()->GetBinCenter(limit_min);
      double ylimit_max = fHistUncertaintyBandPoisson->GetYaxis()->GetBinCenter(limit_max);
      hist->SetBinContent(ix, 0.5*double(ylimit_min+ylimit_max));
      hist->SetBinError(ix, 0.5*double(ylimit_max-ylimit_min));
   }

   return hist;
}

// ---------------------------------------------------------
void BCMTFChannel::PrintHistCumulativeUncertaintyBandPoisson(const char* filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // calculate error band
   this->CalculateHistUncertaintyBandPoisson();

   TH2D hist(*fHistUncertaintyBandPoisson);

   int nbinsx_poisson = hist.GetNbinsX();
   int nbinsy_poisson = hist.GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_p = 0;  // sum of all probabilities inside the interval

      // loop over y-axis of observation
      for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
         double p = hist.GetBinContent(ix, jbin);
         sum_p+=p;
         hist.SetBinContent(ix, jbin, sum_p);
      }
   }

   // draw histogram
   hist.Draw("COLZ");
   c1->Draw();

   // print
   c1->Print(filename);

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFChannel::PrintHistUncertaintyBandPoisson(const char* filename, const char * options)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // calculate error band
   this->CalculateHistUncertaintyBandPoisson();

   // draw histogram
   fHistUncertaintyBandPoisson->Draw(options);
   c1->Draw();

   // print
   c1->Print(filename);

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void BCMTFChannel::PrintUncertaintyBandPoisson(const char* filename, double minimum, double maximum, int color)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // calculate error band
   TH1D* hist=CalculateUncertaintyBandPoisson(minimum, maximum, color);

   // draw histogram
   hist->Draw("E2");
   c1->Draw();

   // print
   c1->Print(filename);

   // free memory
   delete c1;
}

// ---------------------------------------------------------
