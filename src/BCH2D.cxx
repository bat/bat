/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCH2D.h"
#include "BCH1D.h"

#include "BCLog.h"
#include "BCMath.h"

#include <TArrow.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMarker.h>
#include <TObject.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TLine.h>

#include <math.h>
#include <algorithm>

// ---------------------------------------------------------
BCH2D::BCH2D(TH2 * h)
	: BCHistogramBase(h,2)
	, fIntegratedHistogram(0)
	, fBandType(kSmallestInterval)
	, fLogz(false)
	, fDrawProfileX(false)
	, fProfileXType(kProfileMean)
	, fProfileXLineColor(kBlack)
	, fProfileXLineStyle(2)
	, fDrawProfileY(false)
	, fProfileYType(kProfileMean)
	, fProfileYLineColor(kBlack)
	, fProfileYLineStyle(2)
{
	SetDrawLocalMode(true);
}

// ---------------------------------------------------------
BCH2D::BCH2D(const BCH2D & other)
	: BCHistogramBase(other)
	, fIntegratedHistogram(0)
{
	CopyOptions(other);
}

// ---------------------------------------------------------
BCH2D::~BCH2D() {
	if (fIntegratedHistogram)
		delete fIntegratedHistogram;
}

// ---------------------------------------------------------
void BCH2D::CopyOptions(const BCH2D & other) {
	BCHistogramBase::CopyOptions(other);
	fBandType = other.fBandType;
	fLogz = other.fLogz;
	fDrawProfileX = other.fDrawProfileX;
	fProfileXType = other.fProfileXType;
	fProfileXLineColor = other.fProfileXLineColor;
	fProfileXLineStyle = other.fProfileXLineStyle;
	fDrawProfileY = other.fDrawProfileY;
	fProfileYType = other.fProfileYType;
	fProfileYLineColor = other.fProfileYLineColor;
	fProfileYLineStyle = other.fProfileYLineStyle;
}

// ---------------------------------------------------------
void BCH2D::Draw(std::string options, std::vector<double> intervals) {
	if (!((TH2D*)GetHistogram()))
		return;

	// convert options to lowercase
	std::transform(options.begin(),options.end(),options.begin(),::tolower);

	if (intervals.empty())
		intervals = fIntervals;

	// check intervals values
	for (unsigned i=0; i<intervals.size(); ++i) {
		if (intervals[i] < 0 or intervals[i] > 1) {
			BCLog::OutWarning("BCH2D::Draw : intervals values must be in [0,1]. Using defaults.");
			intervals.clear();
			break;
		}
		if (i>1 and fBandType==kSmallestInterval) {
			BCLog::OutWarning("BCH2D::Draw : intervals must be in increasing order for specified band type. Using defaults.");
			intervals.clear();
			break;
		}
	}

	// set defaults if empty
	if (intervals.empty())
		switch (fBandType) {

		case kNoBands:
			break;

		case kSmallestInterval:
		default:
			if (fNBands>0) intervals.push_back(0.682689492137);
			if (fNBands>1) intervals.push_back(0.954499736104);
			if (fNBands>2) intervals.push_back(0.997300203937);
			if (fNBands>3) intervals.push_back(0.999936657516);
			if (fNBands>4) intervals.push_back(0.999999426697);
			if (fNBands>5) intervals.push_back(0.999999998027);
			if (fNBands>6) intervals.push_back(0.999999999997);
			if (fNBands>7) intervals.push_back(1);
			break;

		}

	Smooth(fNSmooth);

	// if option "same" is not specified, draw axes and add "same" to options
	if (options.find("same") == std::string::npos) {
		gPad -> SetLogx(fLogx);
		gPad -> SetLogy(fLogy);
		gPad -> SetLogz(fLogz);
		GetHistogram()->Draw("axis");
		gPad -> Update();
		options += "same";
	}

	if (fBandType == kNoBands or intervals.empty()) {
		GetHistogram() -> Draw("colzsame");
		gPad -> Update();
	} else {

		unsigned nbands = intervals.size();

		// make sure enough colors have been designated
		while (nbands>fBandColors.size())
			fBandColors.push_back(fBandColors.back()-1);

		// set contour levels
		std::vector<double> levels;
		std::vector<std::string> legend_text;
		switch (fBandType) {

		case kSmallestInterval:
		default:
			// fill vector with pairs of bin densities and masses
			std::vector<std::pair<double,double> > bin_dens_vol;
			bin_dens_vol.reserve(fHistogram->GetNbinsX()*fHistogram->GetNbinsY());
			for (int i=1; i<=fHistogram->GetNbinsX(); ++i)
				for (int j=1; j<=fHistogram->GetNbinsY(); ++j)
					if (fHistogram->GetBinContent(i,j)>0)
						bin_dens_vol.push_back(std::make_pair(fHistogram->GetBinContent(i,j),GetHistogram()->Integral(i,i,j,j,"width")));
			if (bin_dens_vol.empty())
				break;
			// sort bins by densities from lowest to highest
			std::stable_sort(bin_dens_vol.begin(),bin_dens_vol.end(),Compare);
			
			// set contour levels by amount of probability mass = 1-interval (in reverse order of intervals)
			// so that contours are of increasing level as needed by ROOT
			double prob_sum = 1;
			for (unsigned i=0; i<bin_dens_vol.size() and levels.size()<intervals.size() and prob_sum>=0; ++i) {
				prob_sum -= bin_dens_vol[i].second;
				while (levels.size()<intervals.size() and prob_sum <= intervals[intervals.size()-levels.size()-1]) {
					levels.push_back(bin_dens_vol[i].first);
					legend_text.push_back(TString::Format("smallest %.1f%% interval(s)",100*(prob_sum+bin_dens_vol[i].second)).Data());
				}
			}
			levels.push_back(bin_dens_vol.back().first);
			break;
		}

		if (levels.size()>1) {

			// set contour levels
			GetHistogram() -> SetContour(levels.size()-1,&levels[0]);

			// set contour colors
			std::vector<int> colors;
			for (int i=levels.size()-2; i>=0; --i)
				colors.push_back(fBandColors[i]);

			gStyle -> SetPalette(colors.size(),&colors[0]);
			
			if (fBandFillStyle<=0)
				fHistogram -> Draw("cont same");
			else {
				fHistogram -> SetFillStyle(fBandFillStyle);
				fHistogram -> Draw("col same");
			}
			gPad -> Update();
			
			// Set legend entries
			switch (fBandType) {
			case kSmallestInterval:
			default:
				for (int i=legend_text.size()-1; i>=0; --i) {
					TLegendEntry* le = fLegend->AddEntry((TObject*)0, legend_text[i].data(), "F");
					if (fBandFillStyle>0) {
						le -> SetFillColor(colors[i]);
						le -> SetFillStyle(fBandFillStyle);
					}
					le -> SetLineColor(colors[i]);
					le -> SetTextAlign(12);
					le -> SetTextFont(62);
					le -> SetTextSize(0.03);
				}
			}

		}
	}
	
	DrawGlobalMode();
	DrawLocalMode();
	DrawMean();
	DrawProfileGraphs();
	
	DrawLegend();
	gPad -> RedrawAxis();

}

// ---------------------------------------------------------
void BCH2D::CalculateIntegratedHistogram()
{
   int nz = 100;

   double zmin = fHistogram->GetMinimum();
   double zmax = fHistogram->GetMaximum();
   double dz   = (zmax-zmin);

   double nx = fHistogram->GetNbinsX();
   double ny = fHistogram->GetNbinsY();

   // create histogram
   if (fIntegratedHistogram)
      delete fIntegratedHistogram;

   fIntegratedHistogram = new TH1D(
      TString::Format("%s_int_prob_%d",fHistogram->GetName(),BCLog::GetHIndex()), "", nz, zmin, zmax+0.05*dz);
   fIntegratedHistogram->SetXTitle("z");
   fIntegratedHistogram->SetYTitle("Integrated probability");
   fIntegratedHistogram->SetStats(kFALSE);

   // loop over histogram
   for (int ix = 1; ix <= nx; ix++) {
      for (int iy = 1; iy <= ny; iy++) {
         double p = fHistogram->GetBinContent(ix, iy);
         fIntegratedHistogram->Fill(p, p);
      }
   }

   // normalize histogram
   fIntegratedHistogram->Scale(1.0/fIntegratedHistogram->GetEntries());

}

// ---------------------------------------------------------
double BCH2D::GetLevel(double p)
{
   double quantiles[1];
   double probsum[1];
   probsum[0] = p;

   fIntegratedHistogram->GetQuantiles( 1, quantiles, probsum);

   return quantiles[0];
}

// ---------------------------------------------------------
double BCH2D::GetArea(double p)
{
   // copy histograms for bands
	TH2D hist_temp(*((TH2D*)fHistogram));

   // calculate total number of bins
   int nbins = hist_temp.GetNbinsX() * hist_temp.GetNbinsY();

   // area
   double area = 0;

   // integrated probability
   double intp = 0;

   // a counter
   int counter = 0;

   // loop over bins
   while ( (intp < p) && (counter < nbins) ) {

      // find maximum bin
      int binx;
      int biny;
      int binz;
      hist_temp.GetBinXYZ(hist_temp.GetMaximumBin(), binx, biny, binz);

      // increase probability
      double dp = hist_temp.GetBinContent(binx, biny);
      intp += dp * hist_temp.GetXaxis()->GetBinWidth(binx) * hist_temp.GetYaxis()->GetBinWidth(biny);

      // reset maximum bin
      hist_temp.SetBinContent(binx, biny, 0.);

      // increase area
      area += hist_temp.GetXaxis()->GetBinWidth(binx) * hist_temp.GetYaxis()->GetBinWidth(biny);

      // increase counter
      counter++;
   }

   // return area
   return area;
}

// ---------------------------------------------------------
std::vector<int> BCH2D::GetNIntervalsY(TH2 * h, int &nfoundmax)
{
   std::vector<int> nint;

   int nx = h->GetNbinsX();
   int ny = h->GetNbinsY();

   nfoundmax=0;

   // loop over histogram bins in x
   for (int ix=1; ix<=nx; ix++)
   {
      int nfound=0;

      // loop over histogram bins in y
      // count nonzero intervals in y
      for (int iy=1; iy<=ny; iy++)
         if(h->GetBinContent(ix,iy)>0.)
         {
            while(h->GetBinContent(ix,++iy)>0.)
               ;
            nfound++;
         }

      // store maximum number of nonzero intervals for the histogram
      if(nfound>nfoundmax)
         nfoundmax=nfound;

      nint.push_back(nfound);
   }

   return nint;
}

// ---------------------------------------------------------
TGraph* BCH2D::CalculateProfileGraph(BCH2DProfileAxis axis, BCH2DProfileType bt) {

	unsigned n_i = (axis==kProfileY) ? GetHistogram()->GetNbinsY() : GetHistogram()->GetNbinsX();
	unsigned n_j = (axis==kProfileY) ? GetHistogram()->GetNbinsX() : GetHistogram()->GetNbinsY();

	TGraph * g = new TGraph();

	// loop over bins of axis to be profiled
	for (unsigned i=1; i<=n_i; ++i) {

		switch (bt) {

		case kProfileMedian: {

			// calculate 50% of total probability mass in slice
			double median_prob = 0.5 * ((axis==kProfileY) ? GetHistogram()->Integral(1,n_j,i,i,"width") : GetHistogram()->Integral(i,i,1,n_j,"width"));
			if (median_prob <= 0)
				break;

			double prob_sum = 0;
			// loop until 50% of probability mass is found
			for (unsigned j=1; j<=n_j; ++j) {
				prob_sum += (axis==kProfileY) ? GetHistogram()->Integral(j,j,i,i,"width") : GetHistogram()->Integral(i,i,j,j,"width");
				if (prob_sum > median_prob) {
					if (axis==kProfileY)
						g -> SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinLowEdge(j),GetHistogram()->GetYaxis()->GetBinCenter(i));
					else 
						g -> SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i),GetHistogram()->GetYaxis()->GetBinLowEdge(j));
					break;
				}
			}
			break;
		}

		case kProfileMode: {
			double max_val = 0;
			unsigned j_max_val = 0;
			for (unsigned j=1; j<=n_j; ++j) {
				double val = (axis==kProfileY) ? GetHistogram()->GetBinContent(j,i) : GetHistogram()->GetBinContent(i,j);
				if (val>max_val) {
					j_max_val = j;
					max_val = val;
				}
			}
			if (j_max_val > 0) {
				if (axis==kProfileY)
					g -> SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(j_max_val),GetHistogram()->GetYaxis()->GetBinCenter(i));
				else 
					g -> SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i),GetHistogram()->GetYaxis()->GetBinCenter(j_max_val));
			}
			break;
		}

		case kProfileMean:
		default: {
			// calculate total probability mass in slice
			double mass_sum = 0;
			double sum = 0 ;
			for (unsigned j=0; j<=n_j; ++j) {
				double mass = (axis==kProfileY) ? GetHistogram()->Integral(j,j,i,i,"width") : GetHistogram()->Integral(i,i,j,j,"width");
				mass_sum += mass;
				sum += mass * ((axis==kProfileY) ? GetHistogram()->GetXaxis()->GetBinCenter(j) : GetHistogram()->GetYaxis()->GetBinCenter(j));
			}
			if (mass_sum >= 0) {
				if (axis==kProfileY)
					g -> SetPoint(g->GetN(), sum/mass_sum, GetHistogram()->GetYaxis()->GetBinCenter(i));
				else
					g -> SetPoint(g->GetN(), GetHistogram()->GetXaxis()->GetBinCenter(i), sum/mass_sum);
			}
			break;
		}
		}
		
	}
	
  // return the graph
  return g;
}

// ---------------------------------------------------------
void BCH2D::DrawProfileGraphs() {
	if (fDrawProfileX) {
		TGraph * graph_profile = CalculateProfileGraph(kProfileX,fProfileXType);
		graph_profile -> SetLineColor(fProfileXLineColor);
		graph_profile -> SetLineStyle(fProfileXLineStyle);
		graph_profile -> Draw("sameL");
		fROOTObjects.push_back(graph_profile);
		std::string xtext = "profile x";
		switch (fProfileXType) {
		case kProfileMode:
			xtext += " (mode)";
			break;
		case kProfileMedian:
			xtext += " (median)";
			break;
		case kProfileMean:
		default:
			xtext += " (mean)";
			break;
		}
		fLegend -> AddEntry(graph_profile, xtext.data(), "L");
	}
	if (fDrawProfileY) {
		TGraph * graph_profile = CalculateProfileGraph(kProfileY,fProfileYType);
		graph_profile -> SetLineColor(fProfileYLineColor);
		graph_profile -> SetLineStyle(fProfileYLineStyle);
		graph_profile -> Draw("sameL");
		fROOTObjects.push_back(graph_profile);
		std::string ytext = "profile y";
		switch (fProfileYType) {
		case kProfileMode:
			ytext += " (mode)";
			break;
		case kProfileMedian:
			ytext += " (median)";
			break;
		case kProfileMean:
		default:
			ytext += " (mean)";
			break;
		}
		fLegend -> AddEntry(graph_profile, ytext.data(), "L");
	}
}

