/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <iostream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
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
void BCMTFChannel::PrintTemplates(const char * filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->Divide(2, 2);

   std::string f(filename);

   // get number of templates
   unsigned int ntemplates = fTemplateContainer.size();

   // calculate number of pages
   unsigned int npages =  ntemplates / 4;
   if (ntemplates % 4 > 0)
      npages++;

   // loop over pages
   for (unsigned int i = 0; i < npages; ++i) {
      // loop over pads
      for (unsigned int j = 0; j < 4; ++j) {
         // calculate template index
         unsigned int templateindex = 4 * i + j;

         if (templateindex < ntemplates && GetTemplate(templateindex)->GetHistogram()) {
            // cd into pad
            c1->cd(j+1);

            // get histogram
            TH1D * hist = new TH1D( *( (TH1D *) GetTemplate(templateindex)->GetHistogram()->Clone() ) );

            // get efficiency
            double efficiency = GetTemplate(templateindex)->GetEfficiency();

            // scale histogram
            hist->Scale(efficiency/hist->Integral());

            // draw histogram
            hist->Draw("HIST");
         }
         else {
            // clear pad
            c1->cd(j+1)->Clear();
         }
      }

      if (npages == 1)
         c1->Print(f.c_str());
      else if (i == 0) {
         c1->Print( (f+std::string("[")).c_str() );
         c1->Print(f.c_str());
      }
      else if (i < npages-1) {
         c1->Print(f.c_str());
      }
      else {
         c1->Print(f.c_str());
         c1->Print( (f+std::string("]")).c_str() );
      }
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
		hist->SetBinContent(ix, 0.5*double(limit_min+limit_max));
		hist->SetBinError(ix, 0.5*double(limit_max-limit_min));
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
void BCMTFChannel::PrintHistUncertaintyBandPoisson(const char* filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

	 // calculate error band
	 this->CalculateHistUncertaintyBandPoisson();

	 // draw histogram
	 fHistUncertaintyBandPoisson->Draw("COLZ");
	 c1->Draw();

	 // print
	 c1->Print(filename);

	 // free memory
	 delete c1;
}

// ---------------------------------------------------------
