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

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TString.h>

#include <math.h>

// ---------------------------------------------------------

BCH2D::BCH2D()
{
   fHistogram = 0;
   fIntegratedHistogram = 0;

   fModeFlag = 0;
}

// ---------------------------------------------------------

BCH2D::BCH2D(TH2D * h)
{
   fHistogram = h;
   fIntegratedHistogram = 0;

   fModeFlag = 0;
}

// ---------------------------------------------------------

BCH2D::~BCH2D()
{
   if (fHistogram) delete fHistogram;
   if (fIntegratedHistogram) delete fIntegratedHistogram;
}

// ---------------------------------------------------------

void BCH2D::Print(const char * filename, int options, int ww, int wh)
{
   // create temporary canvas
   TCanvas * c;
   if(ww >0 && wh > 0)
      c = new TCanvas("c","c",ww,wh);
   else
      c = new TCanvas("c","c",700,700);
   c->cd();

   // draw histogram
   Draw(options);

   gPad->RedrawAxis();

   // print to file
   c->Print(filename);
}

// ---------------------------------------------------------

void BCH2D::Draw(int options, bool drawmode)
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
//               if(val == 0.)
//                  val = log(min)-1.;
//               else
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
//         TMarker * marker = new TMarker(modex, modey, 5);
//         marker->SetMarkerColor(0);
//         marker->Draw();
      }
   }
   else if (options == 53 || options == 531)
   {
      // create new empty histogram
//      int nx = fHistogram->GetNbinsX();
//      int ny = fHistogram->GetNbinsY();
//      double xmin = fHistogram->GetXaxis()->GetXmin();
//      double xmax = fHistogram->GetXaxis()->GetXmax();
//      TH2D * h = new TH2D(
//            TString::Format("htemp53_%d",BCLog::GetHIndex()),fHistogram->GetTitle(),
//            nx,xmin,xmax,
//            ny,fHistogram->GetYaxis()->GetXmin(),fHistogram->GetYaxis()->GetXmax());
//      h->SetXTitle(fHistogram->GetXaxis()->GetTitle());
//      h->SetYTitle(fHistogram->GetYaxis()->GetTitle());

      // copy contents of the main histogram
//      double min = fHistogram->GetMinimum(0.);
/*      for(int i=1;i<=nx;i++)
         for(int j=1;j<=ny;j++)
         {
            double val = fHistogram->GetBinContent(i,j);
            // if requested, change contents of bins to log scale
            if(options == 531)
            {
//               if(val == 0.)
//                  val = log(min)-1.;
//               else
                  val = log10(val);
            }
            h->SetBinContent(i,j,val);
         }

      // draw
      h->SetStats(0);
      h->Draw("colz");
*/
      gPad->SetLogz();
      fHistogram->Draw("colz");

      // draw contour
//      fHistogram->Draw("cont3 same");
//      fHistogram->SetLineWidth(2);

      // set contours
//      CalculateIntegratedHistogram();

//      double levels[] = { GetLevel(0.32) };
//      fHistogram->SetContour(1, levels);

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

void BCH2D::CalculateIntegratedHistogram()
{
   int nz = 100;

   double zmax = fHistogram->GetMaximum();
   double dz   = zmax / double(nz);

   double nx = fHistogram->GetNbinsX();
   double ny = fHistogram->GetNbinsY();

   // create histogram
   if (fIntegratedHistogram)
      delete fIntegratedHistogram;

   fIntegratedHistogram = new TH1D(
         TString::Format("%s_int_prob_%d",fHistogram->GetName(),BCLog::GetHIndex()), "", nz, 0.0, zmax);
   fIntegratedHistogram->SetXTitle("z");
   fIntegratedHistogram->SetYTitle("Integrated probability");
   fIntegratedHistogram->SetStats(kFALSE);

   // loop over histogram
   for (int ix = 1; ix <= nx; ix++)
      for (int iy = 1; iy <= ny; iy++)
      {
         int binmin = BCMath::Nint(fHistogram->GetBinContent(ix, iy) / dz);
         for (int i = binmin; i <= nz; i++)
            fIntegratedHistogram->SetBinContent(i,
                  fIntegratedHistogram->GetBinContent(i) +
                  fHistogram->GetBinContent(ix, iy));
      }
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

// ---------------------------------------------------------
