/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCH1D.h"

#include "BCLog.h"
#include "BCMath.h"

#include <TH1D.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TPaveLabel.h>
#include <TLatex.h>
#include <TError.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TString.h>

#include <math.h>
#include <algorithm>

// ---------------------------------------------------------
BCH1D::BCH1D(TH1 * hist)
	: BCHistogramBase(hist,1)
	, fBandType(kSmallestInterval)
	, fNQuantiles(0)
	, fQuantileLineColor(kBlack)
	, fDrawMedian(false)
	, fDrawCentral68(true)
{
}

// ---------------------------------------------------------
BCH1D::BCH1D(const BCH1D & other) 
	: BCHistogramBase(other)
{
	CopyOptions(other);
}

// ---------------------------------------------------------
BCH1D::~BCH1D() {
}

// ---------------------------------------------------------
void BCH1D::CopyOptions(const BCH1D & other) {
	BCHistogramBase::CopyOptions(other);
	fBandType = other.fBandType;
	fNQuantiles = other.fNQuantiles;
	fQuantileLineColor = other.fQuantileLineColor;
	fDrawMedian = other.fDrawMedian;
	fDrawCentral68 = other.fDrawCentral68;
}

// ---------------------------------------------------------
void BCH1D::SetColorScheme(BCHColorScheme scheme) {
	BCHistogramBase::SetColorScheme(scheme);
	SetQuantileLineColor(GetMarkerColor());
}

// ---------------------------------------------------------
double BCH1D::GetQuantile(double probability)
{
   int nquantiles = 1;
   double quantiles[1];
   double probsum[1];

   probsum[0] = probability;

   // use ROOT function to calculate quantile.
   fHistogram->GetQuantiles(nquantiles, quantiles, probsum);

   return quantiles[0];
}

// ---------------------------------------------------------
void BCH1D::Draw(std::string options, std::vector<double> intervals) {
	// convert optins to lowercase.
	std::transform(options.begin(),options.end(),options.begin(),::tolower);

	if (intervals.empty())
		intervals = fIntervals;

	// check number of intervals values if user-specified
	if (fBandType==kUserSpecified and intervals.size()==1) {
		BCLog::OutWarning("BCH1D::Draw : at least two intervals values must be specified for user-specified intervals. No bands will be drawn.");
		intervals.clear();
	}

	// check interval values
	if (fBandType != kNoBands) 
		for (unsigned i=0; i<intervals.size(); ++i) {
			if (intervals[i] < 0 or intervals[i] > 1) {
				BCLog::OutWarning("BCH1D::Draw : intervals values must be in [0,1]. Using defaults.");
				intervals.clear();
				break;
			}
			if (i>1 and (fBandType==kCentralInterval or fBandType==kSmallestInterval or fBandType==kUpperLimit or fBandType==kUserSpecified) and intervals[i]-intervals[i-1]<=0) {
				BCLog::OutWarning("BCH1D::Draw : intervals must be in increasing order for specified band type. Using defaults.");
				intervals.clear();
				break;
			}
			if (i>1 and fBandType==kLowerLimit and intervals[i]-intervals[i-1]>=0) {
				BCLog::OutWarning("BCH1D::Draw : intervals must be in decreasing order for lower-limit band type. Using defaults.");
				intervals.clear();
				break;
			}
		}

	// set defaults if empty
	if (intervals.empty())
		switch (fBandType) {

		case kNoBands:
			break;

		case kUserSpecified:
			break;
		
		case kUpperLimit:
			if (fNBands>0) intervals.push_back(0.90);
			if (fNBands>1) intervals.push_back(0.95);
			if (fNBands>2) intervals.push_back(0.99);
			for (unsigned i=3; i<fNBands; ++i) // not recommended to go (far) beyond 3
				intervals.push_back(0.9+intervals.back()/10.);
			break;
		
		case kLowerLimit:
			if (fNBands>0) intervals.push_back(0.10);
			if (fNBands>1) intervals.push_back(0.05);
			if (fNBands>2) intervals.push_back(0.01);
			for (unsigned i=3; i<fNBands; ++i) // not recommended to go (far) beyond 3
				intervals.push_back(intervals.back()/10.);
			break;

		case kCentralInterval:
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

	fHistogram -> SetLineColor(kBlack);

	// if option "same" is not specified, draw axes and add "same" to options
	if (options.find("same") == std::string::npos) {
		gPad -> SetLogx(fLogx);
		gPad -> SetLogy(fLogy);
		fHistogram->Draw(options.data());
		gPad -> Update();
		options += "same";
	}

	// draw bands
	if (fBandType!=kNoBands and !intervals.empty() and fHistogram->Integral()>0) {

		unsigned nbands = intervals.size() - ((fBandType==kUserSpecified) ? 1 : 0);

		// make sure enough colors have been designated
		while (nbands>fBandColors.size())
			fBandColors.push_back(fBandColors.back()-1);

		int i0 = fROOTObjects.size();
		for (int i = nbands-1; i >= 0; --i) {

			TH1D * hist_band = 0;
			std::string legend_text;
		
			if (fBandType == kSmallestInterval) {
				hist_band = GetSmallestIntervalHistogram(intervals[i]);
				legend_text = "smallest %.1f%% interval(s)";
			}
			else {
				double p[2] = {0,1};
				switch (fBandType) {

				case kUpperLimit:
					p[1] = intervals[i];
					legend_text = "%.0f%% upper limit";
					break;

				case kLowerLimit:
					p[0] = intervals[i];
					legend_text = "%.0f%% lower limit";
					break;
				 
				case kUserSpecified:
					p[0] = intervals[i];
					p[1] = intervals[i+1];
					legend_text = "%.1f%% interval from " + TString::Format("%.1f%%%% to %.1f%%%%",p[0]*100,p[1]*100);
					break;

				case kCentralInterval:
				default:
					p[0] = 0.5 - intervals[i]/2;
					p[1] = 0.5 + intervals[i]/2;
					legend_text = "central %.1f%% interval";
					break;
				}
				double q[2];
				if (fHistogram->GetQuantiles(2,q,p)==2)
					hist_band = GetSubHisto(q[0],q[1]);
			}

			if (!hist_band)
				continue;

			// set style of band histogram
			hist_band -> SetFillStyle(GetBandFillStyle());
			hist_band -> SetFillColor(GetBandColor(i));
			hist_band -> SetLineColor(0);
			hist_band -> SetLineWidth(0);
			hist_band -> SetLineStyle(0);
			hist_band -> SetTitle(TString::Format(legend_text.data(),hist_band->Integral("width")*100).Data());
			
			// draw band
			fROOTObjects.push_back(hist_band->DrawCopy(options.data()));

			// add to legend
			delete hist_band;
		}
		for (int i=fROOTObjects.size()-1; i>=i0; --i)
			fLegend -> AddEntry(fROOTObjects[i],"","F");
	}
	 
	fHistogram -> Draw(options.data());
	gPad -> Update();

	double ymin = gPad -> GetUymin();
	double ymax = gPad -> GetUymax();
	double ymid = 0.5*(ymin+ymax);
	if (gPad->GetLogy()) {
		ymin = pow(10,ymin);
		ymax = pow(10,ymax);
		ymid = pow(10,ymid);
	}

	// Draw Quantiles
	if (fNQuantiles > 1) {

		// calculate quantile values (q)
		double q[fNQuantiles-1];
		double p[fNQuantiles-1];
		for (unsigned i=1; i<fNQuantiles; ++i)
			p[i-1] = i*1./fNQuantiles;
		if (fHistogram->GetQuantiles(fNQuantiles-1,q,p)==(int)fNQuantiles-1) {
			TLine * quantile_line = new TLine();
			quantile_line -> SetLineStyle(2);
			quantile_line -> SetLineColor(GetQuantileLineColor());
			fROOTObjects.push_back(quantile_line);

			// draw quantile lines
			for (unsigned i=0; i<fNQuantiles-1; ++i)
				quantile_line -> DrawLine(q[i], ymin, q[i], fHistogram->GetBinContent(fHistogram->FindFixBin(q[i])));
			
			std::string quantile_text;
			switch (fNQuantiles) {
			case 2:   quantile_text = "median"; break;
			case 3:   quantile_text = "terciles"; break;
			case 4:   quantile_text = "quartiles"; break;
			case 5:   quantile_text = "quintiles"; break;
			case 6:   quantile_text = "sextiles"; break;
			case 7:   quantile_text = "septiles"; break;
			case 8:   quantile_text = "octiles"; break;
			case 10:  quantile_text = "deciles"; break;
			case 12:  quantile_text = "duodeciles"; break;
			case 20:  quantile_text = "vigintiles"; break;
			case 100: quantile_text = "percentiles"; break;
			default:  quantile_text = TString::Format("%d-quantiles",fNQuantiles); break;
			}
			fLegend->AddEntry(quantile_line, quantile_text.data(), "L");
		}
	}

	DrawGlobalMode();
	DrawLocalMode();
	DrawMean();

	// Median & central 68
	if ( fDrawMedian and fNQuantiles!=2 ) {
		TMarker * marker_median = new TMarker(GetMedian(), ymid*(fLogy ? pow(ymax/ymin,-0.1) : 0.8), 21);
		marker_median -> SetMarkerColor(GetMarkerColor());
		marker_median -> SetMarkerSize(fMarkerScale*gPad->GetWNDC());
		marker_median -> Draw();
		fROOTObjects.push_back(marker_median);

		TLegendEntry * le = 0;
		double q[2], p[2] = {0.1587,0.8413};
			 
		if ( fDrawCentral68 and fHistogram -> GetQuantiles(2,q,p) == 2 ) {
			TArrow* arrow_ci = new TArrow(q[0], ymid*(fLogy ? pow(ymax/ymin,-0.1) : 0.8),
																		q[1], ymid*(fLogy ? pow(ymax/ymin,-0.1) : 0.8),
																		0.02*gPad->GetWNDC(), "<|>");
			arrow_ci -> SetLineColor(marker_median->GetMarkerColor());
			arrow_ci -> SetFillColor(marker_median->GetMarkerColor());
			arrow_ci -> Draw();
			fROOTObjects.push_back(arrow_ci);
			le = fLegend -> AddEntry(arrow_ci,"median and central 68.3% interval","PL");
			le -> SetLineColor(arrow_ci->GetLineColor());
		} else
			le = fLegend -> AddEntry(marker_median,"median","P");
		le -> SetMarkerStyle(marker_median->GetMarkerStyle());
		le -> SetMarkerSize(marker_median->GetMarkerSize());
		le -> SetMarkerColor(marker_median->GetMarkerColor());
	}

	DrawLegend();
	
	gPad -> RedrawAxis();

	return;
}

// ---------------------------------------------------------
double BCH1D::GetSmallestInterval(double & min, double & max, double content)
{
   if(content<=0. || content>= 1.)
      return 0.;

   int nbins = fHistogram->GetNbinsX();

	 if (fHistogram->GetEffectiveEntries()<=0)
		 return 0;

   double factor = fHistogram->Integral("width");
   if(factor==0)
      return 0.;

   // normalize if not yet done
   fHistogram->Scale(1./factor);

   double xmin = fHistogram->GetXaxis()->GetXmin();
   double xmax = fHistogram->GetXaxis()->GetXmax();
   double width = xmax - xmin;

   double xdown=xmin;
   double xup=xmax;

   int ndiv = 10;
   int warn=0;

   double integral_out=0.;

   // loop through the bins
   for(int i=1;i<=nbins;i++)
   {
      if(fHistogram->Integral(i,nbins,"width") < content)
         break;

      // get width of the starting bin
      double firstbinwidth = fHistogram->GetBinWidth(i);

      // divide i-th bin in ndiv divisions
      // integrate starting at each of these divisions
      for(int j=0;j<ndiv;j++)
      {
         double dxdown = (double)(ndiv-j)/(double)ndiv * firstbinwidth;
         xdown = fHistogram->GetBinLowEdge(i) + firstbinwidth - dxdown;
         double integral = dxdown * fHistogram->GetBinContent(i);

         if(integral>content)
         {
            // content fits within one bin
            xup = xdown + content / fHistogram->GetBinContent(i);
            warn = 1;
         }
         else
            for(int k=i+1;k<=nbins;k++)
            {

               double thisbin = fHistogram->GetBinContent(k) * fHistogram->GetBinWidth(k);

               if(integral+thisbin==content)
               {
                  xup = fHistogram->GetBinLowEdge(k+1);
                  break;
               }

               if(integral+thisbin>content)
               {
                  xup = fHistogram->GetBinLowEdge(k) + (content-integral)/thisbin * fHistogram->GetBinWidth(k);
                  integral += thisbin * (xup-fHistogram->GetBinLowEdge(k))/fHistogram->GetBinWidth(k);
                  break;
               }

               integral += thisbin;
            }

         if(integral < content)
            continue;

         if(xup - xdown < width)
         {
            // store necessary information
            width = xup - xdown;
            xmin  = xdown;
            xmax  = xup;
            integral_out = integral;
         }
      }
   }

   if(warn)
   {
      BCLog::OutWarning(
            Form("BCH1D::GetSmallestInterval : The requested content of %d%% fits within one bin.",(int)(content*100)));
      BCLog::OutWarning("BCH1D::GetSmallestInterval : MAKE FINER BINNING (OR CHANGE RANGE) !!!");
   }

   // restore normalization to state before calling this method
   fHistogram->Scale(factor);

   min=xmin;
   max=xmax;

   return integral_out;
}

// ---------------------------------------------------------
double BCH1D::IntegralWidth(double min, double max)
{
	if (min>max)
		return IntegralWidth(max,min);

	 if (fHistogram->GetEffectiveEntries()<=0)
		 return 0;

   int imin = fHistogram->FindBin(min);
   int imax = fHistogram->FindBin(max);
   int nbins = fHistogram->GetNbinsX();

   // if outside of histogram range, return -1.
   if ( imin<1 || imin>nbins || imax<1 || imax>nbins )
      return -1.;

   if ( imin==imax )
      return -1.;

   // calculate first bin
   double first = ( fHistogram->GetBinLowEdge(imin+1) - min ) * fHistogram->GetBinContent(imin);

   // calculate last bin
   double last = ( max - fHistogram->GetBinLowEdge(imax) ) * fHistogram->GetBinContent(imax);

   // calculate inbetween
   double inbetween=0.;
   if(imax-imin>1)
      inbetween = fHistogram->Integral(imin+1, imax-1, "width");

   return first + inbetween + last;
}

// ---------------------------------------------------------
TH1D * BCH1D::GetSubHisto(double min, double max, std::string name, bool preserve_range) {
	if (min==max or !fHistogram)
		return 0;
	if(min>max)
		return GetSubHisto(max,min,name);

	double xmin = fHistogram->GetXaxis()->GetXmin();
	double xmax = fHistogram->GetXaxis()->GetXmax();

	if (max<xmin or min>xmax)
		return 0;
	
	if (name.empty())
		name = std::string(fHistogram->GetName()) + "_subhist";

	if ( min <= xmin and max >= xmax )
		return (TH1D*) fHistogram->Clone(name.data());
	min = std::max<double>(min,xmin);
	max = std::min<double>(max,xmax);
	
	int imin = (min>xmin) ? fHistogram->FindFixBin(min) : 1;
	int imax = (max<xmax) ? fHistogram->FindFixBin(max) : fHistogram->GetNbinsX();

	// create binning for new histogram
	double bins[fHistogram->GetNbinsX()+2];
	bins[0] = (preserve_range) ? xmin : min;
	int i0 = (preserve_range) ? 2 : imin+1;
	int i1 = (preserve_range) ? fHistogram->GetNbinsX() : imax;
	unsigned n=1;
	for (int i=i0; i<=i1; ++i) {
		bins[n++] = fHistogram->GetXaxis()->GetBinLowEdge(i);
		if (min > fHistogram->GetXaxis()->GetBinLowEdge(i) and min < fHistogram->GetXaxis()->GetBinUpEdge(i))
			bins[n++] = min;
		if (max > fHistogram->GetXaxis()->GetBinLowEdge(i) and max < fHistogram->GetXaxis()->GetBinUpEdge(i))
			bins[n++] = max;
	}
	if (preserve_range or max==fHistogram->GetXaxis()->GetBinUpEdge(i1))
		bins[n++] = fHistogram->GetXaxis()->GetBinUpEdge(i1);

	// now define the new histogram
	TH1D * h0 = new TH1D(name.data(),TString::Format("%s;%s;%s",fHistogram->GetTitle(),fHistogram->GetXaxis()->GetTitle(),fHistogram->GetYaxis()->GetTitle()),n-1, bins);
	imin = h0 -> FindFixBin(min);
	imax = h0 -> FindFixBin(max);
	for(int i=imin; i<=imax; ++i)
		h0 -> SetBinContent(i,fHistogram->GetBinContent(fHistogram->FindFixBin(h0->GetBinCenter(i))));
	return h0;
}

// ---------------------------------------------------------
TH1D * BCH1D::GetSmallestIntervalHistogram(double level) {
	if (level < 0)
		return 0;

	if (fHistogram->Integral()<=0)
		return 0;

	TH1D * smallest_int_hist = (TH1D*) fHistogram->Clone(TString::Format("%s_si",fHistogram->GetName()));
	if (level >= 1)
		return smallest_int_hist;
	smallest_int_hist -> Reset(); 
	if (level == 0)
		return smallest_int_hist;

	double sumprob = 0;

	// repeat until enough probability has been gathered
	while (sumprob < level) {
		
		// find bin with max_val, so long has bin has not already been found yet.
		int bin = 0;
		double max_val = 0;
		for (int i=1; i<=fHistogram->GetNbinsX(); ++i) {
			if (smallest_int_hist->GetBinContent(i)!=0)
				continue;
			double val = fHistogram->GetBinContent(i);
			if (val < max_val)
				continue;
			bin = i;
			max_val = val;
		}
		
		smallest_int_hist -> SetBinContent(bin,max_val);
		
		// increase probability sum
		sumprob += max_val * fHistogram->GetXaxis()->GetBinWidth(bin);
	}
	
	return smallest_int_hist;
}

// ---------------------------------------------------------
std::vector<std::vector<double> > BCH1D::GetSmallestIntervals(double content)
{
	std::vector<std::vector<double> > V;
	
	TH1D * hist = GetSmallestIntervalHistogram(content);
	if (!hist)
		return V;
	
	double max = -1;

	std::vector<double> v;
	
	for (int i = 1; i <= hist->GetNbinsX(); ++i) {

		// interval starts here
		if (v.empty() and hist->GetBinContent(i) > 0.) {
			// set xmin
			v.push_back(hist->GetXaxis()->GetBinLowEdge(i));
			// init xmax
			v.push_back(hist->GetXaxis()->GetBinUpEdge(i));
			// init local maximum
			v.push_back(fHistogram->GetBinContent(i));
			v.push_back(hist->GetBinLowEdge(i));
			// init local integral
			v.push_back(0);
		}

		// interval stops here
		if (!v.empty()) {

			// find local maximum
			if (v[2] < fHistogram->GetBinContent(i)) {
				v[2] = fHistogram->GetBinContent(i);
				v[3] = hist->GetBinCenter(i);
			}
			
			// increase area
			v[4] += fHistogram->GetBinContent(i) / fHistogram->Integral();

			if (i==hist->GetNbinsX() or hist->GetBinContent(i+1) <=0) {
				v[1] = hist->GetXaxis()->GetBinUpEdge(i);
				
				// find the absolute maximum
				if (v[2] > max)
					max = v[2];
				V.push_back(v);
				v.clear();
			}
		}
	}
	
	// rescale absolute heights to relative heights
	for (unsigned i = 0; i < V.size(); ++i)
		V[i][2] *= 1. / max;

	delete hist;

	return V;
}
