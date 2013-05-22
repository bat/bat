/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCH2D.h"

#include "BCMath.h"
#include "BCLog.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TObject.h>
#include <TArrow.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TString.h>

#include <math.h>

unsigned int BCH2D::fHCounter=0;

// ---------------------------------------------------------

BCH2D::BCH2D()
: fHistogram(0)
, fIntegratedHistogram(0)
, fROOTObjects(std::vector<TObject*>(0))
{
   fModeFlag = 0;
}

// ---------------------------------------------------------

BCH2D::BCH2D(TH2D * h)
: fIntegratedHistogram(0)
, fROOTObjects(std::vector<TObject*>(0))
{
   fModeFlag = 0;

   SetHistogram(h);
}

// ---------------------------------------------------------

BCH2D::~BCH2D()
{
   if (fHistogram)
      delete fHistogram;
   if (fIntegratedHistogram)
      delete fIntegratedHistogram;

   // clear memory
   int nobjects = (int) fROOTObjects.size();
   for (int i = 0; i < nobjects; ++i)
      if (fROOTObjects[i])
         delete fROOTObjects[i];
}

// ---------------------------------------------------------
void BCH2D::SetColorScheme(int scheme)
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
void BCH2D::SetHistogram(TH2D * hist)
{
   fHistogram = hist;
   if (fHistogram)
      if (fHistogram->Integral()>0)
         fHistogram->Scale(1.0/fHistogram->Integral("width"));
}

// ---------------------------------------------------------
void BCH2D::PrintOld(const char * filename, int options, int ww, int wh)
{
   // create temporary canvas
   TCanvas * c;
   if(ww >0 && wh > 0)
      c = new TCanvas("c","c",ww,wh);
   else
      c = new TCanvas("c","c",700,700);
   c->cd();

   // draw histogram
   DrawOld(options);

   gPad->RedrawAxis();

   // print to file
   c->Print(filename);
}

// ---------------------------------------------------------
void BCH2D::Print(const char * filename, std::string options, std::vector<double> intervals, int ww, int wh)
{
   // option flags
   bool flag_logz = false;
   bool flag_rescale = false;

   // check content of options string
   if (options.find("logz") < options.size()) {
      flag_logz = true;
   }

   if (options.find("R") < options.size()) {
      flag_rescale = true;
   }

   // create temporary canvas
   TCanvas * c;
   unsigned int cindex = getNextIndex();
   if(ww > 0 && wh > 0)
      c = new TCanvas(TString::Format("c_bch2d_%d",cindex), TString::Format("c_bch2d_%d",cindex), ww, wh);
   else
      c = new TCanvas(TString::Format("c_bch2d_%d",cindex));

   // add c to list of objects
   fROOTObjects.push_back(c);

   // set log axis
   if (flag_logz) {
      c->SetLogz();
   }

   // draw histogram
   Draw(options, intervals);

   if (flag_rescale) {
      double top = gPad->GetTopMargin();
      double bottom = gPad->GetBottomMargin();
      double left = gPad->GetLeftMargin();
      double right = gPad->GetRightMargin();

      double dx = 1.-right - left;
      double dy = 1.-top-bottom;
      double ratio = dy/dx;
      double ynew = c->GetWindowWidth()/ratio;
      c->SetCanvasSize(c->GetWindowWidth(), (int) ynew);
      gPad->RedrawAxis();

      c->Modified();
      c->Update();
   }

   // print to file
   c->Print(filename);
}

// ---------------------------------------------------------
void BCH2D::Print(const char* filename, std::string options, double interval, int ww, int wh)
{
   std::vector<double> tempvec;
   tempvec.push_back(interval);
   Print(filename, options, tempvec, ww, wh);
}

// ---------------------------------------------------------
void BCH2D::Draw(std::string options, std::vector<double> intervals)
{
   // option flags
   bool flag_legend = true;
   bool flag_mode_global = false;
   bool flag_mode_local = false;
   bool flag_mean = false;
   bool flag_smooth1 = false;
   bool flag_smooth3 = false;
   bool flag_smooth5 = false;
   bool flag_smooth10 = false;

   // band type
   int bandtype = 0;

   // number of bands
   int nbands = 1; // number of shaded bands
   intervals.push_back(0.6827);

   // define draw options called in TH1D::Draw(...)
   std::string draw_options = "COLZ";

   // check content of options string
   if (options.find("smooth1") < options.size()) {
      flag_smooth1 = true;
   }

   if (options.find("smooth3") < options.size()) {
      flag_smooth3 = true;
   }

   if (options.find("smooth5") < options.size()) {
      flag_smooth5 = true;
   }

   if (options.find("smooth10") < options.size()) {
      flag_smooth10 = true;
   }

   if (options.find("nL") < options.size()) {
      flag_legend = false;
   }

   if (options.find("BTf") < options.size()) {
      bandtype = 0;
   }
   else if (options.find("BTc") < options.size()) {
      bandtype = 1;
   }
   else {
      bandtype = 0;
   }

   if (options.find("gmode") < options.size()) {
      if (fModeFlag)
         flag_mode_global = true;
   }

   if (options.find("lmode") < options.size()) {
      flag_mode_local = true;
   }

   if (options.find("mean") < options.size()) {
      flag_mean = true;
   }

   if (options.find("B1") < options.size()) {
      nbands = 1;
      if (intervals.size() != 1) {
         intervals.clear();
         intervals.push_back(0.6827);
      }
   }

   if (options.find("B2") < options.size()) {
      nbands = 2;
      if (intervals.size() != 2) {
         intervals.clear();
         intervals.push_back(0.6827);
         intervals.push_back(0.9545);
      }
   }

   if (options.find("B3") < options.size()) {
      nbands = 3;
      if (intervals.size() != 3) {
         intervals.clear();
         intervals.push_back(0.6827);
         intervals.push_back(0.9545);
         intervals.push_back(0.9973);
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

   // prepare size of histogram
   double xmin     = fHistogram->GetXaxis()->GetXmin();
   double xmax     = fHistogram->GetXaxis()->GetXmax();
   double ymin     = fHistogram->GetYaxis()->GetXmin();
   double ymaxhist = fHistogram->GetYaxis()->GetXmax();
   double ymax     = ymaxhist;

   // prepare legend
   TLegend* legend = new TLegend();
   legend->SetBorderSize(0);
   legend->SetFillColor(kWhite);
   legend->SetTextAlign(12);
   legend->SetTextFont(62);
   legend->SetTextSize(0.03);

   // add legend to list of objects
   fROOTObjects.push_back(legend);

   // copy histograms for bands
   TH2D* hist_band = new TH2D(*fHistogram);

   // add hist_band to list of ROOT objects
   fROOTObjects.push_back(hist_band);

   // calculate integrated histogram
   CalculateIntegratedHistogram();

   for (int ix = 1; ix <= fHistogram->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= fHistogram->GetNbinsY(); ++iy) {
         double p = fHistogram->GetBinContent(ix, iy);
         hist_band->SetBinContent(ix, iy, p);
      }
   }

   // define levels and colors
   std::vector<double> levels(nbands+2);
   levels[0] = 0.;

   std::vector<int> colors(nbands+1);
   colors[0] = kWhite;

   for (int i = 1; i <= nbands; ++i) {
      levels[i] = GetLevel((1.-intervals[nbands-i]));
      colors[i] = GetColor(nbands-i);
   }
   levels[nbands+1] = fIntegratedHistogram->GetXaxis()->GetXmax();

   // set contour
   hist_band->SetContour(nbands+2, &levels[0]);

   // set colors
   gStyle->SetPalette(nbands+1, &colors[0]);

   // add hist_band to legend
   for (int i = 0; i < nbands; ++i) {
      if (bandtype == 0) {
         TLegendEntry* le = legend->AddEntry((TObject*)0, Form("smallest %.1f%% interval(s)", intervals[nbands-1-i]*100), "F");
         le->SetFillColor(GetColor(nbands-1-i));
         le->SetFillStyle(1001);
         le->SetLineColor(GetColor(nbands-1-i));
         le->SetTextAlign(12);
         le->SetTextFont(62);
         le->SetTextSize(0.03);
      }
      else if (bandtype == 1) {
         TLegendEntry* le = legend->AddEntry((TObject*)0, Form("smallest %.1f%% interval(s)", intervals[nbands-1-i]*100), "F");
         le->SetLineColor(GetColor(nbands-1-i));
         le->SetTextAlign(12);
         le->SetTextFont(62);
         le->SetTextSize(0.03);
      }
   }

   // mean, mode, median
   TMarker* marker_mode_global = new TMarker(fMode[0], fMode[1], 24);
   marker_mode_global->SetMarkerColor(GetColor(4));
   marker_mode_global->SetMarkerSize(1.5);


   int binx, biny, binz;
   fHistogram->GetBinXYZ( fHistogram->GetMaximumBin(), binx, biny, binz);
   TMarker* marker_mode_local = new TMarker(fHistogram->GetXaxis()->GetBinCenter(binx), fHistogram->GetYaxis()->GetBinCenter(biny), 25);
   marker_mode_local->SetMarkerColor(GetColor(4));
   marker_mode_local->SetMarkerSize(1.5);

   double xmean = fHistogram->GetMean(1);
   double ymean = fHistogram->GetMean(2);
   double xrms = fHistogram->GetRMS(1);
   double yrms = fHistogram->GetRMS(2);

   TMarker* marker_mean = new TMarker(xmean, ymean, 32);
   marker_mean->SetMarkerColor(GetColor(4));
   marker_mean->SetMarkerSize(1.5);

   // standard deviation
   TArrow* arrow_std1 = new TArrow(xmean-xrms, ymean,
         xmean+xrms, ymean,
         0.02, "<|>");
   arrow_std1->SetLineColor(GetColor(4));
   arrow_std1->SetFillColor(GetColor(4));

   TArrow* arrow_std2 = new TArrow(xmean, ymean-yrms,
         xmean, ymean+yrms,
         0.02, "<|>");
   arrow_std2->SetLineColor(GetColor(4));
   arrow_std2->SetFillColor(GetColor(4));

   // add marker_mean and arrow_std to list of ROOT objects
   fROOTObjects.push_back(marker_mean);
   fROOTObjects.push_back(marker_mode_global);
   fROOTObjects.push_back(marker_mode_local);
   fROOTObjects.push_back(arrow_std1);
   fROOTObjects.push_back(arrow_std2);

   if (flag_mode_global) {
      TLegendEntry* le = legend->AddEntry(marker_mode_global, "global mode", "P");
      le->SetMarkerStyle(24);
      le->SetMarkerSize(1.5);
      le->SetMarkerColor(GetColor(4));
   }

   if (flag_mode_local) {
      TLegendEntry* le = legend->AddEntry(marker_mode_local, "local mode", "P");
      le->SetMarkerStyle(25);
      le->SetMarkerSize(1.5);
      le->SetMarkerColor(GetColor(4));
   }

   if (flag_mean) {
      TLegendEntry* le = legend->AddEntry(arrow_std1, "mean and standard deviation", "PL");
      le->SetLineColor(GetColor(4));
      le->SetMarkerStyle(32);
      le->SetMarkerSize(1.5);
      le->SetMarkerColor(GetColor(4));
   }

   // calculate legend height in NDC coordinates
   double height = 0.03*legend->GetNRows();

   // make room for legend
   if (flag_legend)
      ymax+=(ymax-ymin)*(0.1+height);

   double deltax = 0.0015*(xmax - xmin);
   double deltay = 0.0015*(ymaxhist - ymin);
   TH2D* hist_axes = new TH2D("", "", 1, xmin-deltax, xmax+deltax, 1, ymin-deltay, ymaxhist+deltay);
   hist_axes->SetXTitle(fHistogram->GetXaxis()->GetTitle());
   hist_axes->SetYTitle(fHistogram->GetYaxis()->GetTitle());
   hist_axes->SetLineWidth(fHistogram->GetLineWidth());
   hist_axes->SetStats(kFALSE);
   fROOTObjects.push_back(hist_axes);

   // draw axes
   hist_axes->Draw("COL");

   // smooth
   if (flag_smooth1) {
      fHistogram->Smooth(1);
      fHistogram->Scale(1.0/fHistogram->Integral("width"));
   }
   if (flag_smooth3) {
      fHistogram->Smooth(3);
      fHistogram->Scale(1.0/fHistogram->Integral("width"));
   }
   if (flag_smooth5) {
      fHistogram->Smooth(5);
      fHistogram->Scale(1.0/fHistogram->Integral("width"));
   }
   if (flag_smooth10) {
      fHistogram->Smooth(10);
      fHistogram->Scale(1.0/fHistogram->Integral("width"));
   }

   // draw histogram
   if (bandtype == 0)
      hist_band->Draw("COL SAME");
   else if (bandtype == 1)
      hist_band->Draw("CONT1 SAME");

   // mean, mode, median
   if (flag_mode_global) {
      marker_mode_global->Draw();
   }

   if (flag_mode_local) {
      marker_mode_local->Draw();
   }

   if (flag_mean) {
      arrow_std1->Draw();
      arrow_std2->Draw();
      marker_mean->Draw();
   }

   if (flag_legend)
      gPad->SetTopMargin(0.02);

   double xlegend1 = gPad->GetLeftMargin();
   double xlegend2 = 1.0-gPad->GetRightMargin();
   double ylegend1 = 1.-gPad->GetTopMargin()-height;
   double ylegend2 = 1.-gPad->GetTopMargin();

   // place legend on top of histogram
   legend->SetX1NDC(xlegend1);
   legend->SetX2NDC(xlegend2);
   legend->SetY1NDC(ylegend1);
   legend->SetY2NDC(ylegend2);

   // draw legend
   if (flag_legend) {
      legend->Draw();
   }

   hist_axes->Draw("SAME");

   gPad->SetTopMargin(1.-ylegend1+0.01);

   gPad->RedrawAxis();

   return;
}

// ---------------------------------------------------------
void BCH2D::Draw(std::string options, double interval)
{
   std::vector<double> tempvec;
   tempvec.push_back(interval);
   Draw(options, tempvec);
}

// ---------------------------------------------------------

void BCH2D::DrawOld(int options, bool drawmode)
{
   // draw histogram
   fHistogram->SetLineColor(kBlack);
   //   fHistogram->SetLineWidth(4);

   fHistogram->Scale(1./fHistogram->Integral("width"));

   double modex,modey;

   // set mode to display
   if(fModeFlag)
   {
      modex = fMode[0];
      modey = fMode[1];
   }
   else
   {
      int maximumbin = fHistogram->GetMaximumBin();

      int binx = maximumbin % (fHistogram->GetNbinsX() + 2);
      int biny = maximumbin / (fHistogram->GetNbinsX() + 2);

      modex = fHistogram->GetXaxis()->GetBinCenter(binx);
      modey = fHistogram->GetYaxis()->GetBinCenter(biny);
   }

   // normalize histogram
   fHistogram->Scale(1./fHistogram->Integral("width"));

   // draw according to selected option
   if (options == 0)
      fHistogram->Draw("CONT0");
   else if (options == 1)
   {
      fHistogram->Draw("CONT3");

      // set contours
      CalculateIntegratedHistogram();

      double levels[4];
      levels[0] = 0.;
      levels[1] = GetLevel(1.0 - 0.6827);
      levels[2] = GetLevel(1.0 - 0.9545);
      levels[3] = GetLevel(1.0 - 0.9973);

      fHistogram->SetContour(4, levels);

      // best fit value
      TMarker* marker = new TMarker(modex, modey, 24);
      marker->Draw();

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->AddEntry(fHistogram, "68% prob. region", "L");
      legend->AddEntry(marker, "Best fit", "P");
      legend->Draw();
   }
   else if (options == 2)
   {
      fHistogram->Draw("CONT3");

      // set contours
      CalculateIntegratedHistogram();

      double levels[2];
      double level32 = GetLevel(0.32);
      levels[0] = 0.;
      levels[1] = level32;

      fHistogram->SetContour(2, levels);

      // best fit value
      TMarker* marker = new TMarker(modex, modey, 24);
      marker->Draw();

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->AddEntry(fHistogram, "68% prob. region", "L");
      legend->AddEntry(marker, "Best fit", "P");
      legend->Draw();

   }
   else if (options == 3)
   {
      fHistogram->Draw("CONT3");

      // set contours
      CalculateIntegratedHistogram();

      double levels[2];
      double level10 = GetLevel(0.10);
      levels[0] = 0.;
      levels[1] = level10;

      fHistogram->SetContour(2, levels);

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->AddEntry(fHistogram, "90% prob. region", "L");
      legend->Draw();
   }
   else if (options == 4)
   {
      fHistogram->Draw("CONT3");

      // set contours
      CalculateIntegratedHistogram();

      double levels[2];
      double level5 = GetLevel(0.05);
      levels[0] = 0.;
      levels[1] = level5;

      fHistogram->SetContour(2, levels);

      TLegend* legend = new TLegend(0.65, 0.80, 0.95, 0.95);
      legend->SetBorderSize(0);
      legend->SetFillColor(kWhite);
      legend->AddEntry(fHistogram, "95% prob. region", "L");
      legend->Draw();
   }
   else if (options == 5)
      fHistogram->Draw("COL");
   else if (options == 52 || options == 521)
   {
      // create new empty histogram
      int nx = fHistogram->GetNbinsX();
      int ny = fHistogram->GetNbinsY();
      double xmin = fHistogram->GetXaxis()->GetXmin();
      double xmax = fHistogram->GetXaxis()->GetXmax();
      TH2D * h = new TH2D(
            TString::Format("htemp52_%d",BCLog::GetHIndex()),fHistogram->GetTitle(),
            nx,xmin,xmax,
            ny,fHistogram->GetYaxis()->GetXmin(),fHistogram->GetYaxis()->GetXmax());
      h->SetXTitle(fHistogram->GetXaxis()->GetTitle());
      h->SetYTitle(fHistogram->GetYaxis()->GetTitle());

      // copy contents of the main histogram
      //      double min = fHistogram->GetMinimum(0.);
      for(int i=1;i<=nx;i++)
         for(int j=1;j<=ny;j++)
         {
            double val = fHistogram->GetBinContent(i,j);
            // if requested, change contents of bins to log scale
            if(options == 521)
            {
               val = log10(val);
            }
            h->SetBinContent(i,j,val);
         }

      // draw
      h->SetStats(0);
      h->Draw("col");

      // draw contour
      fHistogram->Draw("cont3 same");
      fHistogram->SetLineWidth(2);

      // set contours
      CalculateIntegratedHistogram();

      double levels[] = { GetLevel(0.32) };
      fHistogram->SetContour(1, levels);

      // best fit value
      if(drawmode)
      {
         TMarker * marker0 = new TMarker(modex, modey, 8);
         marker0->SetMarkerColor(0);
         marker0->SetMarkerSize(.7);
         marker0->Draw();
         TMarker * marker1 = new TMarker(modex, modey, 4);
         marker1->SetMarkerColor(1);
         marker1->SetMarkerSize(.7);
         marker1->Draw();
      }
   }
   else if (options == 53 || options == 531)
   {
      gPad->SetLogz();
      fHistogram->Draw("colz");

      // best fit value
      TMarker * marker0 = new TMarker(modex, modey, 8);
      marker0->SetMarkerColor(0);
      marker0->SetMarkerSize(.7);
      marker0->Draw();
      TMarker * marker1 = new TMarker(modex, modey, 4);
      marker1->SetMarkerColor(1);
      marker1->SetMarkerSize(.7);
      marker1->Draw();
      //      TMarker * marker = new TMarker(modex, modey, 5);
      //      marker->SetMarkerColor(0);
      //      marker->Draw();
   }
}

// ---------------------------------------------------------
void BCH2D::PrintIntegratedHistogram(const char* filename)
{
   TCanvas* c = new TCanvas();
   c->cd();
   c->SetLogy();
   CalculateIntegratedHistogram();
   fIntegratedHistogram->Draw();
   c->Print(filename);

   delete c;
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
	TH2D hist_temp(*fHistogram);
	
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

TGraph ** BCH2D::GetBandGraphs(TH2D * h, int &n)
{
   n=0;

   int nbands=0;
   TH2D * hcopy = (TH2D*)h->Clone(TString::Format("%s_copy_%d",h->GetName(),BCLog::GetHIndex()));

   std::vector<int> nint=GetNIntervalsY(hcopy,nbands);

   if(nbands>2)
   {
      BCLog::OutError(
            Form("BCH2D::GetBandGraphs : Detected %d bands. Maximum allowed is 2 (sorry).",nbands));
      return 0;
   }
   else if(nbands==0)
   {
      BCLog::OutError("BCH2D::GetBandGraphs : No bands detected.");
      return 0;
   }

   TGraph ** gxx = new TGraph*[nbands];

   TH2D * h0 = (TH2D*)h->Clone();

   if (nbands>0)
      gxx[0] = GetLowestBandGraph(h0,nint);

   if (nbands==2)
   {
      gxx[1] = GetLowestBandGraph(h0);
      n=2;
   }
   else
      n=1;

   return gxx;
}

// ---------------------------------------------------------

std::vector<int> BCH2D::GetNIntervalsY(TH2D * h, int &nfoundmax)
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

TGraph * BCH2D::GetLowestBandGraph(TH2D * h)
{
   int n;
   return GetLowestBandGraph(h,GetNIntervalsY(h,n));
}

// ---------------------------------------------------------

TGraph * BCH2D::GetLowestBandGraph(TH2D * h, std::vector<int> nint)
{
   int nx = h->GetNbinsX() - 1;
   int ny = h->GetNbinsY();

   TGraph * g = new TGraph(2*nx);

   for (int ix=1; ix<=nx; ix++)
   {
      // get x for the bin
      double x;
      if(ix==1)
         x = h->GetXaxis()->GetBinLowEdge(1);
      else if(ix==nx)
         x = h->GetXaxis()->GetBinLowEdge(nx+1);
      else
         x = h->GetXaxis()->GetBinCenter(ix);

      for(int iy=1; iy<=ny; iy++)
         if(h->GetBinContent(ix,iy)>0.)
         {
            // get low edge of the first not empty bin in y
            g->SetPoint(ix-1, x, h->GetYaxis()->GetBinLowEdge(iy));

            // delete content of all subsequent not empty bins
            if(nint[ix-1]==2)
               h->SetBinContent(ix,iy,0.);

            while(h->GetBinContent(ix,++iy)>0.)
               if(nint[ix-1]==2)
                  h->SetBinContent(ix,iy,0.);

            // get low edge of the first empty bin in y
            g->SetPoint(2*nx-ix, x, h->GetYaxis()->GetBinLowEdge(iy));

            break;
         }
   }

   return g;
}

// ---------------------------------------------------------

std::vector<double> BCH2D::GetLevelBoundary(double level)
{
   return GetLevelBoundary(fHistogram, level);
}

// ---------------------------------------------------------

std::vector<double> BCH2D::GetLevelBoundary(TH2D * h, double level)
{
   std::vector<double> b;

   int nx = h->GetNbinsX();

   b.assign(nx - 1, 0.0);

   // loop over x and y bins.
   for (int ix = 1; ix < nx; ix++)
   {
      TH1D * h1 = h->ProjectionY(TString::Format("temphist_%d",BCLog::GetHIndex()), ix, ix);

      int nprobSum = 1;
      double q[1];
      double probSum[] = { level };

      h1->GetQuantiles(nprobSum, q, probSum);

      b[ix-1] = q[0];
   }

   return b;
}

// ---------------------------------------------------------

TGraph * BCH2D::GetBandGraph(double l1, double l2)
{
   return GetBandGraph(fHistogram , l1, l2);
}

// ---------------------------------------------------------

TGraph * BCH2D::GetBandGraph(TH2D * h , double l1, double l2)
{
   // define new graph
   int nx = h->GetNbinsX() - 1;

   TGraph * g = new TGraph(2*nx);
   //   g->SetFillStyle(1001);
   //   g->SetFillColor(kYellow);

   // get error bands
   std::vector<double> ymin = GetLevelBoundary(h,l1);
   std::vector<double> ymax = GetLevelBoundary(h,l2);

   for (int i = 0; i < nx; i++)
   {
      g->SetPoint(i, h->GetXaxis()->GetBinCenter(i+1), ymin[i]);
      g->SetPoint(nx+i, h->GetXaxis()->GetBinCenter(nx-i), ymax[nx-i-1]);
   }

   return g;
}
