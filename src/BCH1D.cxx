/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCH1D.h"

#include "BCLog.h"
#include "BCMath.h"

#include <TH1D.h>
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

// ---------------------------------------------------------

BCH1D::BCH1D()
  : fHistogram(0)
  , fROOTObjects(std::vector<TObject*>(0))
{
   fDefaultCLLimit = 95.; // in percent

   fModeFlag = 0;
}

// ---------------------------------------------------------

BCH1D::~BCH1D()
{
   if (fHistogram) delete fHistogram;

   // clear memory
   int nobjects = (int) fROOTObjects.size();
   for (int i = 0; i < nobjects; ++i) 
     if (fROOTObjects[i])
       delete fROOTObjects[i];
}

// ---------------------------------------------------------

double BCH1D::GetMode()
{
   return fHistogram->GetBinCenter(fHistogram->GetMaximumBin());
}

// ---------------------------------------------------------

double BCH1D::GetQuantile(double probability)
{
   int nquantiles = 1;
   double quantiles[1];
   double probsum[1];

   probsum[0] = probability;

   // use ROOT function to calculat quantile.
   fHistogram->GetQuantiles(nquantiles, quantiles, probsum);

   return quantiles[0];
}

// ---------------------------------------------------------

double BCH1D::GetIntegral(double valuemin, double valuemax)
{
   double integral = 0;

   int binmin = fHistogram->FindBin(valuemin);
   int binmax = fHistogram->FindBin(valuemax);

   // use ROOT function to calculate integral.
   integral = fHistogram->Integral(binmin, binmax);

   return integral;
}

// ---------------------------------------------------------

double BCH1D::GetPValue(double probability)
{
   // use ROOT function to calculate the integral from 0 to "probability".
   return fHistogram->Integral(1, fHistogram->FindBin(probability));
}

// ---------------------------------------------------------

void BCH1D::SetDefaultCLLimit(double limit)
{
   // check if limit is between 68% and 100%. Give a warning if not ...
   if(limit>=100. || limit<68.)
      BCLog::OutWarning(
         Form("BCH1D::SetDefaultCLLimit : Value %.1f out of range. Keeping default %.1f%% CL limit.",limit,fDefaultCLLimit));

   // ... or set value if true.
   else
      fDefaultCLLimit = limit;
}

// ---------------------------------------------------------
void BCH1D::SetColorScheme(int scheme)
{
  fColors.clear();
	
  // color numbering
  // 0,1,2 : intervals
  // 3 : quantile lines
	// 4 : mean, mode, median

  if (scheme == 0) {
    fColors.push_back(18);
    fColors.push_back(16);
    fColors.push_back(14);
    fColors.push_back(kBlack);
    fColors.push_back(kBlack);
  }
  else if (scheme == 1) {
    fColors.push_back(kGreen);
    fColors.push_back(kYellow);
    fColors.push_back(kRed);
    fColors.push_back(kBlack);
    fColors.push_back(kBlack);
  }
  else if (scheme == 2) {
    fColors.push_back(kBlue+4);
    fColors.push_back(kBlue+2);
    fColors.push_back(kBlue);
    fColors.push_back(kOrange);
    fColors.push_back(kOrange);
  }
  else if (scheme == 3) {
    fColors.push_back(kRed+4);
    fColors.push_back(kRed+2);
    fColors.push_back(kRed);
    fColors.push_back(kGreen);
    fColors.push_back(kGreen);
  }
  else {
    SetColorScheme(1);
  }
}

// ---------------------------------------------------------
void BCH1D::Print(const char * filename, int options, double ovalue, int ww, int wh)
{
   char file[256];
   int i=0;
   while(i<255 && *filename!='\0')
      file[i++]=*filename++;
   file[i]='\0';

// temporary
   fHistogram->SetLineWidth(1);

   // create temporary canvas
   TCanvas * c;
   if(ww > 0 && wh > 0)
      c = new TCanvas("c","c",ww,wh);
   else
      c = new TCanvas("c","c");

   c->cd();
   Draw(options, ovalue);

   gPad->RedrawAxis();

   // print to file.
   c->Print(file);
}

// ---------------------------------------------------------

void BCH1D::myPrint(const char* filename, std::string options, std::vector<double> intervals, int ww, int wh)
{
  char file[256];
  int i=0;
  while(i<255 && *filename!='\0')
    file[i++]=*filename++;
  file[i]='\0';

  // option flags
  bool flag_logx = false;
  bool flag_logy = false;

  // check content of options string
  if (options.find("logx") < options.size()) {
    flag_logx = true;
  }

  if (options.find("logy") < options.size()) {
    flag_logy = true;
  }

  // create temporary canvas
  TCanvas * ctemp;
  if(ww > 0 && wh > 0) {
    ctemp = new TCanvas("","",ww,wh);
  }
  else
    ctemp = new TCanvas("","");

  // add ctemp to list of objects
  fROOTObjects.push_back(ctemp);

  ctemp->cd();

  // set log axis 
  if (flag_logx) {
    ctemp->SetLogx();
  }

  // set log axis
  if (flag_logy) {
    ctemp->SetLogy();
  }

  myDraw(options, intervals);

	// debugKK
	double top = gPad->GetTopMargin();
	double bottom = gPad->GetBottomMargin();
	double left = gPad->GetLeftMargin();
	double right = gPad->GetRightMargin();
	
	double dx = 1.-right - left;
	double dy = 1.-top-bottom;
	double ratio = dy/dx;
	double ynew = ratio * ctemp->GetWindowWidth();
	ctemp->SetWindowSize(ctemp->GetWindowWidth(), ynew);

  gPad->RedrawAxis();

  // print to file.
  ctemp->Print(file);
}

// ---------------------------------------------------------
void BCH1D::myPrint(const char* filename, std::string options, double interval, int ww, int wh)
{
  std::vector<double> tempvec;
  tempvec.push_back(interval);
  myPrint(filename, options, tempvec, ww, wh);
}

// ---------------------------------------------------------
void BCH1D::Draw(int options, double ovalue)
{
   double min, max;
   double mode;
   double thismode = GetMode();

   int nbins = fHistogram->GetNbinsX();

   fHistogram->Scale(1./fHistogram->Integral("width"));

	 // todoKK:
   // - fix legend size, position and style
   // - maybe write a painer class?
   // options should start with letter indicating type, e.g.
	 // L: legend
	 //    - add legend
   //    - add written summary on plot (mean = xyz, ...)
   // B: band (or shaded area)
   //    - band type: smallest, central, ...
   //    - draw shaded band between two arbitrary values
   //    - add second (and third) band with arbitrary values, e.g. 90%, 95% prob. regions
   // D: drawing
	 //    - draw histogram as normal histogram
	 //    - draw histogram as smooth curve
   //    - black&white version with hatching
   //    - draw in log-scale
   //    - draw cumulative pdf
   // S: summary values (mean, rms, ...)
   //    - add indicator for mode, mean, median individually
   //    - add indicator for rms, smallest interval, central interval individually 
   //    - add indicator for quantiles (decentiles, quartiles) e.g. as red lines going through the whole histogram

   if(fModeFlag)
      mode=fMode;
   else
      mode=thismode;

   // define temporary line.
   TLine * line;

   //control if legend is required
   bool flagLegend=false;
   char confidenceLabel[25];

   // reset histogram
   fHistogram->SetLineColor(kBlack);
   fHistogram->SetLineWidth(1);
   fHistogram->SetFillStyle(0);

   // check drawing option.
   switch(options)
   {
      // Draw a band between 16% and 84% probability.
      // If the mode is outside the band only draw a limit.
      case 0:
         if (fabs(ovalue) >= 100 || ovalue==0.)
         {//default case if no args to Draw() supplied

            min = GetQuantile(.16);
            max = GetQuantile(.84);

            //draw a legend later
            flagLegend = true;
            sprintf(confidenceLabel, "Central 68%%");

            if ( fHistogram->FindBin(thismode) == fHistogram->GetNbinsX() )
            {
               min = GetQuantile(1.-(double)fDefaultCLLimit/100.);
               max = fHistogram->GetXaxis()->GetXmax();
               ovalue = fDefaultCLLimit;
               sprintf(confidenceLabel, "%g%% region", fDefaultCLLimit);
            }
            else if ( fHistogram->FindBin(thismode)==1)
            {
               min = fHistogram->GetXaxis()->GetXmin();
               max = GetQuantile((double)fDefaultCLLimit/100.);
               ovalue = -fDefaultCLLimit;
               sprintf(confidenceLabel, "%g%% region", fDefaultCLLimit);
            }
         }

         else if(ovalue < 0)
         {
            min = fHistogram->GetXaxis()->GetXmin();
            max = GetQuantile(-ovalue/100.);
         }
         else
         {
            min = GetQuantile(1. - ovalue/100.);
            max = fHistogram->GetXaxis()->GetXmax();
         }

         // do the drawing
         DrawShadedLimits(mode, min, max, ovalue);

         // add legend for the symbols mean, mode, median, confidence band
         if(flagLegend)
            this ->DrawLegend(confidenceLabel);

         break;

      // Draw a line at "ovalue".
      case 1:

         fHistogram->Draw();
         min = fHistogram->GetBinLowEdge(1);
         max = fHistogram->GetBinLowEdge(nbins+1);
         if(min<=ovalue && ovalue<=max)
         {
            line = new TLine();
            line->SetLineColor(kRed);
            line->DrawLine(ovalue, 0., ovalue, 1.05 * fHistogram->GetMaximum());
         }

         break;

      // Draw a shaded band at the smallest interval.
      case 2:

         if(ovalue<50) // just to ensure there's some sense in the number
            ovalue = 68.; // default is 68%

         GetSmallestInterval(min, max, ovalue/100.);
         DrawShadedLimits(mode, min, max, 0.);

         break;

      // Draw a shaded band at the smallest intervals
      case 3:

         if(ovalue<50) // just to ensure there's some sense in the number
            ovalue = 68.; // default is 68%

         DrawSmallest(mode,ovalue);

         break;

      // Draw just one bin for fixed delta prior
      case 4:

          DrawDelta(ovalue);

          break;

      // Sort out bad options and warn.
      default:

         BCLog::OutError(Form("BCH1D::Draw : Invalid option %d",options));
         break;
   }
}

// ---------------------------------------------------------
void BCH1D::myDraw(std::string options, std::vector<double> intervals)
{
  // todoKK:
  // plot cumulative pdf
  // plot probability density/probability

  // option flags
  bool flag_pdf0 = true;
  bool flag_pdf1 = false;
  bool flag_legend = false;
  bool flag_logx;
  bool flag_logy;
  bool flag_mode = false;
  bool flag_median = false;
  bool flag_mean = false;
  bool flag_quartiles = false;
  bool flag_deciles = false;
  bool flag_percentiles = false;

  // band type
  int bandtype = 0;

  // number of bands
  int nbands = 0; // number of shaded bands

	// define draw options called in TH1D::Draw(...)
	std::string draw_options = ""; 

  // check content of options string
  if (options.find("logx") < options.size()) {
    flag_logx = true;
  }

  if (options.find("logy") < options.size()) {
    flag_logy = true;
  }

  if (options.find("L") < options.size()) {
    flag_legend = true;
  }

  if (options.find("D1") < options.size()) {
    draw_options.append("C");
  }

  if (options.find("BTsi") < options.size()) {
    bandtype = 1;
  }
  else if (options.find("BTul") < options.size()) {
    bandtype = 2;
  }
  else if (options.find("BTll") < options.size()) {
    bandtype = 3;
  }
  else {
    bandtype = 0;
  }

  if (options.find("B1") < options.size()) {
    nbands = 1;
    if (bandtype == 0 && intervals.size() != 2) {
      intervals.clear();
      intervals.push_back(0.1587);
      intervals.push_back(0.8413);
    }
    else if (bandtype == 1 && intervals.size() != 1) {
      intervals.clear();
      intervals.push_back(0.6827);
    }
    else if (bandtype == 2 && intervals.size() != 1) {
      intervals.clear();
      intervals.push_back(0.90);
    }
    else if (bandtype == 3 && intervals.size() != 1) {
      intervals.clear();
      intervals.push_back(0.10);
    }
  }

  if (options.find("B2") < options.size()) {
    nbands = 2;
    if (bandtype == 0 && intervals.size() != 4) {
      intervals.clear();
      intervals.push_back(0.1587);
      intervals.push_back(0.8413);
      intervals.push_back(0.0228);
      intervals.push_back(0.9772);
    }
    else if (bandtype == 1 && intervals.size() != 2) {
      intervals.clear();
      intervals.push_back(0.6827);
      intervals.push_back(0.9545);
    }
    else if (bandtype == 2 && intervals.size() != 2) {
      intervals.clear();
      intervals.push_back(0.90);
      intervals.push_back(0.95);
    }
    else if (bandtype == 3 && intervals.size() != 2) {
      intervals.clear();
      intervals.push_back(0.10);
      intervals.push_back(0.05);
    }
  }
	
  if (options.find("B3") < options.size()) {
    nbands = 3;
    if (bandtype == 0 && intervals.size() != 6) {
      intervals.clear();
      intervals.push_back(0.1587);
      intervals.push_back(0.8413);
      intervals.push_back(0.0228);
      intervals.push_back(0.9772);
      intervals.push_back(0.0013);
      intervals.push_back(0.9987);
    }
    else if (bandtype == 1 && intervals.size() != 3) {
      intervals.clear();
      intervals.push_back(0.6827);
      intervals.push_back(0.9545);
      intervals.push_back(0.9973);
    }
    else if (bandtype == 2 && intervals.size() != 3) {
      intervals.clear();
      intervals.push_back(0.90);
      intervals.push_back(0.95);
      intervals.push_back(0.99);
    }
    else if (bandtype == 3 && intervals.size() != 3) {
      intervals.clear();
      intervals.push_back(0.10);
      intervals.push_back(0.05);
      intervals.push_back(0.01);
    }
  }

  if (options.find("CS0") < options.size()) {
    SetColorScheme(0);
  }
  else if (options.find("CS1") < options.size()) {
    SetColorScheme(1);
  }
  else if (options.find("CS2") < options.size()) {
    SetColorScheme(2);
  }
  else if (options.find("CS3") < options.size()) {
    SetColorScheme(3);
  }
  else  {
    SetColorScheme(1);
  }
	
  if (options.find("pdf0") < options.size()) {
    flag_pdf0 = true;
    flag_pdf1 = false;
  }

  if (options.find("pdf1") < options.size()) {
    flag_pdf0 = false;
    flag_pdf1 = true;
  }

  if (options.find("mode") < options.size()) {
    if (fModeFlag)
			flag_mode = true;
  }
  if (options.find("median") < options.size()) {
    flag_median = true;
  }
  if (options.find("mean") < options.size()) {
    flag_mean = true;
  }
  if (options.find("quartiles") < options.size()) {
    flag_quartiles = true;
  }
  if (options.find("deciles") < options.size()) {
    flag_deciles = true;
  }
  if (options.find("percentiles") < options.size()) {
    flag_percentiles = true;
  }

  // normalize histogram to unity
  fHistogram->Scale(1./fHistogram->Integral("width"));
	
  // prepare legend
  TLegend* legend = new TLegend();
  legend->SetBorderSize(0);
	legend->SetFillColor(kWhite);
  legend->SetTextAlign(12);
  legend->SetTextFont(62);
  legend->SetTextSize(0.03);

  // add legend to list of objects
  fROOTObjects.push_back(legend);

  // draw histogram
  if (flag_pdf0)
    fHistogram->Draw(draw_options.c_str());

  // draw bands
  for (int i = 0; i < nbands; ++i) {
    int col = GetColor(nbands-i-1);

    double prob_low  = 0;
    double prob_high = 0;
    double prob_interval = 0;
    double xlow  = 0;
    double xhigh = 0;

    TH1D * hist_band = 0;

    if (bandtype == 0) {
      prob_low  = intervals[2*(nbands-i)-2];
      prob_high = intervals[2*(nbands-i)-1];
      xlow  = GetQuantile(prob_low);
      xhigh = GetQuantile(prob_high);
      prob_interval = prob_high - prob_low;

      hist_band = GetSubHisto(xlow, xhigh, TString::Format("%s_sub_%d", fHistogram->GetName(), BCLog::GetHIndex()));
    }
    else if (bandtype == 1) {
      prob_interval = GetSmallestInterval(xlow, xhigh, intervals[nbands-1-i]);
      hist_band = GetSmallestIntervalHistogram(intervals[nbands-1-i]);
      for (int ibin = 1; ibin < hist_band->GetNbinsX(); ++ibin)
	hist_band->SetBinContent(ibin, hist_band->GetBinContent(ibin)*fHistogram->GetBinContent(ibin));
    }
    else if(bandtype == 2) {
      xlow = 0.;
      xhigh = GetQuantile(intervals[nbands-1-i]);
      hist_band = GetSubHisto(xlow, xhigh, TString::Format("%s_sub_%d", fHistogram->GetName(), BCLog::GetHIndex()));
      prob_interval = intervals[nbands-1-i];
    }
    else if(bandtype == 3) {
      xlow = GetQuantile(intervals[nbands-1-i]);
      xhigh = GetQuantile(1.);
      hist_band = GetSubHisto(xlow, xhigh, TString::Format("%s_sub_%d", fHistogram->GetName(), BCLog::GetHIndex()));
      prob_interval = 1.-intervals[nbands-1-i];
    }

    // set style of band histogram
		hist_band->SetFillStyle(1001);
    hist_band->SetFillColor(col);
    hist_band->SetLineColor(col);

    // draw shaded histogram
    hist_band->Draw(std::string(draw_options+std::string("same")).c_str());

		// draw histogram again
		if (flag_pdf0)
			fHistogram->Draw(std::string(std::string("SAME")+draw_options).c_str());

    // add to legend
    std::string legend_label;
    if (bandtype == 0)
      legend_label.append(Form("central %.1f%% interval ", prob_interval*100));
    else if (bandtype == 1)
      legend_label.append(Form("smallest %.1f%% interval(s)", prob_interval*100));
    else if (bandtype == 2)
      legend_label.append(Form("%.0f%% upper limit", prob_interval*100));
    else if (bandtype == 3)
      legend_label.append(Form("%.0f%% lower limit", prob_interval*100));
    
    legend->AddEntry(hist_band, legend_label.c_str(), "F");

    // add hist_band to list of objects
    fROOTObjects.push_back(hist_band);
  }

  gPad->RedrawAxis();

  // prepare size of histogram
	double xmin     = fHistogram->GetXaxis()->GetXmin();
	double xmax     = fHistogram->GetXaxis()->GetXmax();
  double ymin     = 0;
  double ymaxhist = fHistogram->GetBinContent(fHistogram->GetMaximumBin());
  double ymax     = ymaxhist;
  double xfraction = 1.-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin();
  double yfraction = 1.-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin(); 

  // check if log axis
  if (flag_logy)
    ymin = 1e-4*ymaxhist;

  // quantiles
  TLine* line_quantiles = new TLine();
  line_quantiles->SetLineStyle(2);
  line_quantiles->SetLineColor(GetColor(3));

  if (flag_quartiles) {
    for (int i = 1; i < 4; ++i) {
      double quantile_x = GetQuantile(0.25*i);
      int quantile_xbin = fHistogram->FindBin(quantile_x);
      double quantile_y = fHistogram->GetBinContent(quantile_xbin);
      double quantile_ymin = 0;
      if (flag_logy)
				quantile_ymin = 1e-4*ymaxhist;
      line_quantiles->DrawLine(quantile_x, quantile_ymin, quantile_x, quantile_y);
    }
    TLegendEntry* le = legend->AddEntry(line_quantiles, "quartiles", "L");
		if (nbands>0) 
			le->SetFillColor(GetColor(0));
		else
			le->SetFillColor(kWhite);
		le->SetFillStyle(1001);
  }

  if (flag_deciles) {
    for (int i = 1; i < 10; ++i) {
      double quantile_x = GetQuantile(0.10*i);
      int quantile_xbin = fHistogram->FindBin(quantile_x);
      double quantile_y = fHistogram->GetBinContent(quantile_xbin);
      double quantile_ymin = 0;
      if (flag_logy)
				quantile_ymin = 1e-4*ymaxhist;
      line_quantiles->DrawLine(quantile_x, quantile_ymin, quantile_x, quantile_y);
    }
    TLegendEntry* le = legend->AddEntry(line_quantiles, "deciles", "FL");
		le->SetFillColor(GetColor(0));
		le->SetFillStyle(1001);
  }
  
  if (flag_percentiles) {
    for (int i = 1; i < 100; ++i) {
      double quantile_x = GetQuantile(0.01*i);
      int quantile_xbin = fHistogram->FindBin(quantile_x);
      double quantile_y = fHistogram->GetBinContent(quantile_xbin);
      double quantile_ymin = 0;
      if (flag_logy)
				quantile_ymin = 1e-4*ymaxhist;
      line_quantiles->DrawLine(quantile_x, quantile_ymin, quantile_x, quantile_y);
    }
    TLegendEntry* le = legend->AddEntry(line_quantiles, "percentiles", "L");		
		le->SetFillColor(GetColor(0));
		le->SetFillStyle(1001);
  }

  // add line_quantiles to list of ROOT objects
  fROOTObjects.push_back(line_quantiles);

  // mean, mode, median
	TMarker* marker_mode = new TMarker(fMode, 0.50*ymaxhist, 24);
	marker_mode->SetMarkerColor(GetColor(4));
	marker_mode->SetMarkerSize(1.5);

	TMarker* marker_mean = new TMarker(GetMean(), 0.55*ymaxhist, 20);
	marker_mean->SetMarkerColor(GetColor(4));
	marker_mean->SetMarkerSize(1.5);

	TMarker* marker_median = new TMarker(GetMedian(), 0.45*ymaxhist, 21);
	marker_median->SetMarkerColor(GetColor(4));
	marker_median->SetMarkerSize(1.5);

	// standard deviation
	TArrow* arrow_std = new TArrow(GetMean()-GetRMS(), 0.55*ymaxhist,
																 GetMean()+GetRMS(), 0.55*ymaxhist,
																 0.02, "<|>");
	arrow_std->SetLineColor(GetColor(4));
	arrow_std->SetFillColor(GetColor(4));

	// central interval
	TArrow* arrow_ci = new TArrow(GetQuantile(0.1587), 0.45*ymaxhist,
																GetQuantile(0.8413), 0.45*ymaxhist,
																0.02, "<|>");
	arrow_ci->SetLineColor(GetColor(4));
	arrow_ci->SetFillColor(GetColor(4));

	// add marker_mean and arrow_std to list of ROOT objects
	fROOTObjects.push_back(marker_mean);
	fROOTObjects.push_back(marker_median);
	fROOTObjects.push_back(arrow_std);
	fROOTObjects.push_back(arrow_ci);

	if (flag_mode) {
		marker_mode->Draw();
		TLegendEntry* le = legend->AddEntry(marker_mode, "global mode", "P");
		le->SetMarkerStyle(24);
		le->SetMarkerSize(1.5);
		le->SetMarkerColor(GetColor(4));
	}

	if (flag_mean) {
		arrow_std->Draw();
		marker_mean->Draw();
		TLegendEntry* le = legend->AddEntry(arrow_std, "mean and standard deviation", "PL");
		le->SetLineColor(GetColor(4));
		le->SetMarkerStyle(20);
		le->SetMarkerSize(1.5);
		le->SetMarkerColor(GetColor(4));
	}
	
	if (flag_median) {
		arrow_ci->Draw();
		marker_median->Draw();
		TLegendEntry* le = legend->AddEntry(arrow_ci, "median and central 68.3% interval", "PL");
		le->SetLineColor(GetColor(4));
		le->SetMarkerStyle(21);
		le->SetMarkerSize(1.5);
		le->SetMarkerColor(GetColor(4));
	}
	
  // calculate legend height in NDC coordinates
  double height = 0.05*legend->GetNRows();

  // make room for legend
  if (flag_legend)
    ymax*=(1.15+height);
  else 
    ymax*=1.1;

  fHistogram->GetYaxis()->SetRangeUser(ymin, 1.05*ymaxhist);

  // calculate dimensions in NDC variables
	/*
  double xlegend1 = gStyle->GetPadLeftMargin()+0.05*xfraction;
  double xlegend2 = gStyle->GetPadLeftMargin()+0.95*xfraction;
  double ylegend1 = gStyle->GetPadBottomMargin() + 1.10*ymaxhist/ymax*yfraction;
	double ylegend2 = gStyle->GetPadBottomMargin() + (ymax-0.05*ymaxhist)/ymax*yfraction;
	*/
	// debugKK
	if (flag_legend)
		gStyle->SetPadTopMargin(0.02);

  double xlegend1 = gStyle->GetPadLeftMargin();
  double xlegend2 = 1.0-gStyle->GetPadRightMargin();
  double ylegend1 = 1.-gStyle->GetPadTopMargin()-height;
	double ylegend2 = 1.-gStyle->GetPadTopMargin();

  // place legend on top of histogram
  legend->SetX1NDC(xlegend1);
  legend->SetX2NDC(xlegend2);
  legend->SetY1NDC(ylegend1);
  legend->SetY2NDC(ylegend2);

  // draw legend
  if (flag_legend) {
    legend->Draw();
  }

	// draw line to separate legend
	/*
	if (flag_legend) {
		TLine* line_boundary = new TLine();
		line_boundary->SetLineColor(kBlack);
		line_boundary->DrawLine(xmin, 1.05*ymaxhist,
														xmax, 1.05*ymaxhist);
		fROOTObjects.push_back(line_boundary);
	}
	*/

	// rescale top margin
	double cm_top = gPad->GetTopMargin();
	double margin_fraction = (ymax - 1.05*ymaxhist)/ymax*yfraction; 

	gPad->SetTopMargin(1.-ylegend1+0.01);

  gPad->RedrawAxis();

  return;
}

// ---------------------------------------------------------
void BCH1D::myDraw(std::string options, double interval)
{
  std::vector<double> tempvec;
  tempvec.push_back(interval);
  myDraw(options, tempvec);  
}

// ---------------------------------------------------------
void BCH1D::DrawDelta(double value) const
{
    // draw histogram with axes first
    double xmin = fHistogram->GetXaxis()->GetXmin();
    double xmax = fHistogram->GetXaxis()->GetXmax();
    double ysize = 1.3 * fHistogram->GetMaximum();

    TH2D*  hsc = new TH2D(
          TString::Format("h2scale_%s_%d",fHistogram->GetName(),BCLog::GetHIndex()),"",
          50, xmin, xmax, 10, 0., ysize);
    hsc->SetStats(0);
    hsc->SetXTitle(fHistogram->GetXaxis()->GetTitle());
    hsc->SetYTitle(fHistogram->GetYaxis()->GetTitle());
    hsc->Draw();

    // write mode location

   TLatex * tmax_text = new TLatex();
   tmax_text->SetTextSize(0.035);
   tmax_text->SetTextFont(62);
   tmax_text->SetTextAlign(22); // center

   double xprint=(xmax+xmin)/2.;
   double yprint=ysize*(1-1.4*tmax_text->GetTextSize());

   tmax_text->DrawLatex(xprint,yprint, TString::Format("%s = %g", fHistogram->GetXaxis()->GetTitle(), value));
   delete tmax_text;

   // create a temporary histogram, to hold only one non-zero bin
   TH1D * hist_temp = new TH1D(TString::Format("h1scale_%s_%d", fHistogram->GetName(), BCLog::GetHIndex()), "", 100, xmin, xmax);
   hist_temp->SetBinContent(hist_temp->FindBin(value), fHistogram->GetMaximum());
   hist_temp->Draw("same");
}

// ---------------------------------------------------------
void BCH1D::DrawLegend(const char* text)
{
   //draw on top right corner

   TLegend* legend = new TLegend(0.73, 0.72, 0.86, 0.85);
   legend->SetFillColor(kWhite);
   legend->SetBorderSize(1);

   TMarker* triangle = new TMarker(0, 0, 23);
   triangle->SetMarkerColor(kRed);
   legend->AddEntry(triangle, "Global mode", "P");

   TMarker* diamond = new TMarker(0, 0, 27);
   diamond->SetMarkerColor(kBlue);
   legend->AddEntry(diamond, "Mean", "P");

   TLine * line;
   line = new TLine();
   line->SetLineStyle(2);
   line->SetLineColor(kGreen + 2);
   legend->AddEntry(line, "Median", "l");

   TLegend* band = new TLegend(0, 0, 1, 1);
   band->SetFillColor(kYellow);
   legend->AddEntry(band, text, "F");

   legend->SetTextAlign(12);

   //fine tuned by hand so text legible even with several plots in a row
   legend->SetTextSize(0.016);

   legend->Draw();
}

// ---------------------------------------------------------

// TODO Are graphics objects ever deleted from the heap? search for new TH*
// In most cases, plotting won't work as expected if the histograms are deleted in here
// as the TCanvas then remains empty.

void BCH1D::DrawShadedLimits(double mode, double min, double max, double limit)
{
   double maximum = fHistogram->GetMaximum();

   double x0 = mode;
   double y0 = fHistogram->GetBinContent( fHistogram->FindBin(mode) );

   double x1 = fHistogram->GetMean();
   double y1 = fHistogram->GetBinContent( fHistogram->FindBin(x1) );

   double x2 = GetQuantile(.5); // median
   double y2 = fHistogram->GetBinContent( fHistogram->FindBin(x2) );

   double ysize = maximum*1.3;

   double xmin = fHistogram->GetXaxis()->GetXmin();
   double xmax = fHistogram->GetXaxis()->GetXmax();

   // draw histogram with axes first
   TH2D * hsc = new TH2D(
         TString::Format("h2scale_%s_%d",fHistogram->GetName(),BCLog::GetHIndex()),"",
         10, xmin, xmax, 10, 0., ysize);
   hsc->SetStats(0);
   hsc->SetXTitle(fHistogram->GetXaxis()->GetTitle());
   hsc->SetYTitle(fHistogram->GetYaxis()->GetTitle());
   hsc->Draw();

   // draw histogram
   fHistogram->SetLineWidth(1);
   fHistogram->Draw("same");

   // draw yellow shaded region between min and max
   TH1D * hist_shaded = GetSubHisto(min,max,TString::Format("%s_sub_%d",fHistogram->GetName(),BCLog::GetHIndex()));
   hist_shaded->SetFillStyle(1001);
   hist_shaded->SetFillColor(kYellow);

   // draw shaded histogram
   hist_shaded->Draw("same");

   gPad->RedrawAxis();

   // draw triangle for mode
   TPolyLine * tmax;

   double dx = 0.01*(xmax-xmin);
   double dy = 0.04*(ysize);
   y0+=dy/5.;
   double tmax_x[] = {x0, x0 + dx, x0 - dx, x0};
   double tmax_y[] = {y0, y0 + dy, y0 + dy, y0};
   tmax = new TPolyLine(4,tmax_x,tmax_y);
   tmax->SetLineColor(kRed);
   tmax->SetLineWidth(1);
   tmax->SetFillColor(kRed);
   tmax->Draw();
   tmax->Draw("f");

   // draw diamond for mean
   TPolyLine * tmean;

   y1+=dy/5.;
//   double tmean_x[] = {x1, x1 + dx, x1 - dx, x1};
//   double tmean_y[] = {y1, y1 + dy, y1 + dy, y1};
   double tmean_x[] = {x1, x1 + dx, x1 , x1 - dx, x1};
   double tmean_y[] = {y1, y1 + dy/2., y1 + dy, y1 + dy/2., y1};
   tmean = new TPolyLine(5,tmean_x,tmean_y);
   tmean->SetLineColor(kBlue);
//   tmean->SetFillColor(kWhite);
   tmean->SetLineWidth(1);
//   tmean->SetLineStyle(1);
   tmean->Draw();

/*
   // draw arrow for median
   TPolyLine * tmed;
   TPolyLine * tmed2;

//   y2+=dy/5.;
   y2=0.+dy/5.;
   double tmed_x[] = {x2 + dx, x2, x2 - dx};
   double tmed_y[] = {y2 + dy, y2, y2 + dy};
   double tmed2_x[] = {x2, x2};
   double tmed2_y[] = {y2, y2 + dy*2.};
   tmed = new TPolyLine(3,tmed_x,tmed_y);
   tmed2 = new TPolyLine(2,tmed2_x,tmed2_y);
   tmed->SetLineColor(kGreen+2);
   tmed->SetLineWidth(1);
   tmed2->SetLineColor(kGreen+2);
   tmed2->SetLineWidth(1);
   tmed->Draw();
   tmed2->Draw();
*/
   // draw dashed line for median
   TLine * line;
   line = new TLine();
   line->SetLineStyle(2);
   line->SetLineColor(kGreen+2);
   line->DrawLine(x2, 0., x2, y2);


   // write mode location and shaded band

   // format of the number
   double delta_max = fmax(fabs(max-x1),fabs(x1-min));

   int sd = 2 + (int)log10(fabs(x1/delta_max));

   if( (int)log10(x1) > (int)log10(delta_max) )
      sd++;

   TLatex * tmax_text = new TLatex();
   tmax_text->SetTextSize(0.035);
   tmax_text->SetTextFont(62);
   tmax_text->SetTextAlign(22); // center
//   tmax_text->SetTextAlign(13); // top-left

   double xprint=(xmax+xmin)/2.;
//   double xprint=xmin+(xmax-xmin)/20.;
   double yprint=ysize*(1-1.4*tmax_text->GetTextSize());

   if(fabs(limit)<50) // just to ensure there's some sense in the number
      tmax_text->DrawLatex(xprint,yprint,
         TString::Format( TString::Format("%%s^{med} = %%.%dg ^{+%%.2g}_{ -%%.2g}",sd),
            fHistogram->GetXaxis()->GetTitle(), x2, max-x2, x2-min));

   else if (limit<0)
      tmax_text->DrawLatex(xprint,yprint,
         TString::Format( TString::Format("%%s (%d%%%% prob.) < %%.4g",-(int)limit),
            fHistogram->GetXaxis()->GetTitle(), max));

   else if (limit>0)
      tmax_text->DrawLatex(xprint,yprint,
         TString::Format( TString::Format("%%s (%d%%%% prob.) > %%.4g",(int)limit),
            fHistogram->GetXaxis()->GetTitle(), min));

/*
   TLegend * leg = new TLegend(.61,.7,.9,.88);
   leg->SetFillColor(kWhite);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->AddEntry(line,"Median", "l");
   TH1D * hh0 = new TH1D("hh0","",10,0,10);
   hh0->SetMarkerStyle(23);
   hh0->SetMarkerColor(kBlue);
   hh0->SetMarkerSize(1);
   leg->AddEntry(hh0,"Mean","p");
   TH1D * hh1 = new TH1D("hh1","",10,0,10);
   hh1->SetMarkerStyle(23);
   hh1->SetMarkerColor(kRed);
   hh1->SetMarkerSize(1);
   leg->AddEntry(hh1,"Mode","p");
   leg->AddEntry(hist_shaded, "Central 68%", "f");
   leg ->Draw();
*/

}

// ---------------------------------------------------------

void BCH1D::DrawSmallest(double mode, double prob, bool drawmean)
{
   // prepare triangle for mode
   TPolyLine * tmax;

   double x0 = mode;
   double y0 = fHistogram->GetBinContent( fHistogram->FindBin(mode) );
   double xmin = fHistogram->GetXaxis()->GetXmin();
   double xmax = fHistogram->GetXaxis()->GetXmax();
   double ysize = 1.3 * fHistogram->GetMaximum();

   double x1 = fHistogram->GetMean();
   double y1 = fHistogram->GetBinContent( fHistogram->FindBin(x1) );

   double x2 = GetQuantile(.5); // median
   double y2 = fHistogram->GetBinContent( fHistogram->FindBin(x2) );

   double dx = 0.01*(xmax-xmin);
   double dy = 0.04*(ysize);
   double tmax_x[] = {x0, x0 + dx, x0 - dx, x0};
   double tmax_y[] = {y0, y0 + dy, y0 + dy, y0};
   tmax = new TPolyLine(4,tmax_x,tmax_y);
   tmax->SetLineColor(kRed);
   tmax->SetFillColor(kRed);

   // draw histogram with axes first
   TH2D * hsc = new TH2D(
         TString::Format("h2scale_%s_%d",fHistogram->GetName(),BCLog::GetHIndex()),"",
         10, xmin, xmax, 10, 0., ysize);
   hsc->SetStats(0);
   hsc->SetXTitle(fHistogram->GetXaxis()->GetTitle());
   hsc->SetYTitle(fHistogram->GetYaxis()->GetTitle());
   hsc->Draw();

   // histogram to be filled with band
   TH1D * hist_temp1 = (TH1D*) fHistogram->Clone();
   hist_temp1->Scale(1.0/fHistogram->Integral("width"));
   hist_temp1->SetFillColor(kYellow);
   hist_temp1->SetFillStyle(1001);

   // temporary histogram
   TH1D * hist_temp2 = (TH1D*) fHistogram->Clone();
   hist_temp2->Scale(1.0/fHistogram->Integral("width"));

   // clear content
   hist_temp1->Reset();

   // loop over original histogram and copy bin untils integral equals
   // "prob"
   prob /= 100.;

   // integral
   double sum = 0.0;

   int n=0;
   int nbins=fHistogram->GetNbinsX();

   // loop
   while (sum < prob && n < nbins)
   {
      n++;

      // find bin with maximum
      int bin = hist_temp2->GetMaximumBin();

      // copy bin to new histogram
      double val = hist_temp2->GetBinContent(bin);
      hist_temp1->SetBinContent(bin, val);

      // remove maximum from temporary histogram
      hist_temp2->SetBinContent(bin, 0.0);

      // integrate by summing
      sum += val * hist_temp2->GetBinWidth(bin);
   }

   // scale histogram
   hist_temp1->Scale(fHistogram->Integral("width"));

   // draw histograms
   fHistogram->Draw("same");
   hist_temp1->Draw("same");

   // draw triangle for mode
   tmax->Draw("f");

   if(drawmean)
   {
      // draw triangle for mean
      // draw diamond for mean
      TPolyLine * tmean;

      y1+=dy/5.;
//      double tmean_x[] = {x1, x1 + dx, x1 - dx, x1};
//      double tmean_y[] = {y1, y1 + dy, y1 + dy, y1};
      double tmean_x[] = {x1, x1 + dx, x1 , x1 - dx, x1};
      double tmean_y[] = {y1, y1 + dy/2., y1 + dy, y1 + dy/2., y1};
      tmean = new TPolyLine(5,tmean_x,tmean_y);
      tmean->SetLineColor(kBlue);
//      tmean->SetFillColor(kWhite);
      tmean->SetLineWidth(1);
//      tmean->SetLineStyle(1);
      tmean->Draw();

      // draw dashed line for median
      TLine * line;
      line = new TLine();
      line->SetLineStyle(2);
      line->SetLineColor(kGreen+2);
      line->DrawLine(x2, 0., x2, y2);
   }

   // free memory
   delete hist_temp2;
}

// ---------------------------------------------------------

double BCH1D::GetSmallestInterval(double & min, double & max, double content)
{
   if(content<=0. || content>= 1.)
      return 0.;

   int nbins = fHistogram->GetNbinsX();

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
//   if(nbins<100)
//      ndiv = 1000;
//   if(nbins>1000)
//      ndiv = 10;

   int warn=0;

   double integral_out=0.;

   // loop through the bins
   for(int i=1;i<nbins+1;i++)
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
            for(int k=i+1;k<nbins+1;k++)
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
   int imin = fHistogram->FindBin(min);
   int imax = fHistogram->FindBin(max);

   int nbins = fHistogram->GetNbinsX();

   // if outside of histogram range, return -1.
   if ( imin<1 || imin>nbins || imax<1 || imax>nbins )
      return -1.;

   if ( imin==imax )
      return -1.;

   // swap values if necessary
   if (imin>imax)
   {
      int i=imin;
      double x=min;
      imin=imax, min=max;
      imax=i, max=x;
   }

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

TH1D * BCH1D::GetSubHisto(double min, double max, const char * name)
{
   if(min>max)
   {
      double t=min;
      min=max;
      max=t;
   }

   int imin = fHistogram->FindBin(min);
   int imax = fHistogram->FindBin(max);

   int nbins = fHistogram->GetNbinsX();
   double xmin = fHistogram->GetXaxis()->GetXmin();
   double xmax = fHistogram->GetXaxis()->GetXmax();

   if( min==max || (min<=xmin && max>=xmax) )
   {
      TH1D * h0 = (TH1D*) fHistogram->Clone(name);
      return h0;
   }

   if (imin<1)
   {
      imin=1;
      min=xmin;
   }
   if (imax>nbins)
   {
      imax=nbins;
      max=xmax;
   }

   double * xb = new double[nbins+3]; // nbins+1 original bin edges + 2 new bins
   int n=0; // counter

   int domin=1;
   int domax=1;

   for (int i=1;i<nbins+2;i++)
   {
      double x0 = fHistogram->GetBinLowEdge(i);

      if (min<x0 && domin)
      {
         xb[n++]=min;
         domin=0;
      }
      else if (min==x0)
         domin=0;

      if (max<x0 && domax)
      {
         xb[n++]=max;
         domax=0;
      }
      else if (max==x0)
         domax=0;

      xb[n++]=x0;
   }

   // now define the new histogram
   TH1D * h0 = new TH1D(name,"",n-1,xb);
   for(int i=1;i<n;i++)
   {
      double x0 = h0->GetBinCenter(i);
      if(x0<min || x0>max)
         continue;

      int bin=fHistogram->FindBin(x0);
      h0->SetBinContent(i, fHistogram->GetBinContent(bin));
   }

   return h0;
}

// ---------------------------------------------------------

TH1D * BCH1D::GetSmallestIntervalHistogram(double level)
{
   // create a new histogram which will be returned and all yellow
   TH1D * hist_yellow = (TH1D*) fHistogram->Clone();
   hist_yellow->Reset();
   hist_yellow->SetFillStyle(1001);
   hist_yellow->SetFillColor(kYellow);

   // copy a temporary histogram
   TH1D * hist_temp = (TH1D*) fHistogram->Clone(TString::Format("%s_%i",fHistogram->GetName(),BCLog::GetHIndex()));
   double factor = hist_temp->Integral("");

   if(factor == 0)
      return 0;

   hist_temp->Scale(1. / factor);

   // here's the algorithm:
   // 1. find the maximum bin in the temporary histogram and copy it to
   // the yellow histogram.
   // 2. remove this bin from the temporary histogram.
   // 3. repeat this until a total of "level" probability is copied to
   // the yellow histogram.

   double sumprob = 0.0;

   while (sumprob < level)
   {
      // find maximum bin and value
      int bin = hist_temp->GetMaximumBin();
      double value = hist_temp->GetMaximum();

      // copy "1" into the corresponding bin in the yellow histogram
      hist_yellow->SetBinContent(bin, 1.0);

      // set the bin value in the temporary histogram to zero
      hist_temp->SetBinContent(bin, 0.0);

      // increase probability sum
      sumprob += value;
   }

   delete hist_temp;

   return hist_yellow;
}

// ---------------------------------------------------------

std::vector<double> BCH1D::GetSmallestIntervals(double content)
{
   std::vector<double> v;

   TH1D * hist = GetSmallestIntervalHistogram(content);

   int nbins = hist->GetNbinsX();
   int ninter = 0;
   int lastbin = -1;

   double max = -1;
   double localmax = -1;
   double localmaxpos = -1;
   double localint = 0;
   bool flag = false;

   for (int i = 1; i <= nbins; ++i)
   {
      // interval starts here
      if (!flag && hist->GetBinContent(i) > 0.)
      {
         flag = true;
         v.push_back(hist->GetBinLowEdge(i));

         // remember start position of the interval
         lastbin = i;

         // increase number of intervals
         ninter++;

         // reset local maximum
         localmax = fHistogram->GetBinContent(i);
         localmaxpos = hist->GetBinLowEdge(i);

         // reset local integral
         localint = 0;
      }

      // interval stops here
      if ((flag && !(hist->GetBinContent(i) > 0.)) || (flag && i == nbins))
      {
         flag = false;
         v.push_back(hist->GetBinLowEdge(i) + hist->GetBinWidth(i));

         // set right bin to maximum if on edge
         if (i == nbins && localmax < fHistogram->GetBinContent(i))
            localmaxpos = hist->GetBinCenter(i) + 0.5 * hist->GetBinWidth(i);

         // find the absolute maximum
         if (localmax > max)
            max = localmax;

         // save local maximum
         v.push_back(localmax);
         v.push_back(localmaxpos);

         // save local integral
         v.push_back(localint);
      }

      // find local maximum
      if (i < nbins && localmax < fHistogram->GetBinContent(i))
      {
         localmax = fHistogram->GetBinContent(i);
         localmaxpos = hist->GetBinCenter(i);
      }

      // increase area
      localint += fHistogram->GetBinContent(i) / fHistogram->Integral();
   }

   // rescale absolute heights to relative heights
   for (int i = 0; i < ninter; ++i)
      v[i*5+2] = v.at(i*5+2) / max;

   return v;
}

// ---------------------------------------------------------
