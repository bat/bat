/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCSummaryTool.h"
#include <string>

#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCModel.h"
#include "BCParameter.h"
#include "BCSummaryPriorModel.h"

#include <TArrow.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPostScript.h>
#include <TStyle.h>

#include <fstream>
#include <iostream>

unsigned int BCSummaryTool::fHCounter=0;

// ---------------------------------------------------------
BCSummaryTool::BCSummaryTool()
   : fModel(0)
   , fPriorModel(0)
   , fFlagInfoMarg(false)
   , fFlagInfoOpt(false)
{
   // define sum of probabilities for quantiles
   fSumProb.push_back(0.05);
   fSumProb.push_back(0.10);
   fSumProb.push_back(0.1587);
   fSumProb.push_back(0.50);
   fSumProb.push_back(0.8413);
   fSumProb.push_back(0.90);
   fSumProb.push_back(0.95);

   // set text style
   gStyle->SetPaintTextFormat(".2g");
}

// ---------------------------------------------------------
BCSummaryTool::BCSummaryTool(BCModel * model)
   : fModel(model)
   , fPriorModel(0)
   , fFlagInfoMarg(false)
   , fFlagInfoOpt(false)
{
   // define sum of probabilities for quantiles
   fSumProb.push_back(0.05);
   fSumProb.push_back(0.10);
   fSumProb.push_back(0.1587);
   fSumProb.push_back(0.50);
   fSumProb.push_back(0.8413);
   fSumProb.push_back(0.90);
   fSumProb.push_back(0.95);

   // set text style
   gStyle->SetPaintTextFormat(".2g");
}

// ---------------------------------------------------------
BCSummaryTool::~BCSummaryTool()
{
   delete fPriorModel;
}

// ---------------------------------------------------------
int BCSummaryTool::CopySummaryData()
{
   // check if model exists
   if (!fModel)
      return 0;

   // clear information
   fParName.clear();
   fParMin.clear();
   fParMax.clear();
   fMean.clear();
   fMargMode.clear();
   fGlobalMode.clear();
   fQuantiles.clear();
   fSmallInt.clear();
   fRMS.clear();
   fCorrCoeff.clear();

   // copy information from marginalized distributions
   for (unsigned i = 0; i < fModel->GetNParameters(); ++i) {

      // copy parameter information
      fParName.push_back( (fModel->GetParameter(i)->GetLatexName()) );
      fParMin.push_back( fModel->GetParameter(i)->GetLowerLimit() );
      fParMax.push_back( fModel->GetParameter(i)->GetUpperLimit() );

      // copy 1D marginalized information
      BCH1D * bch1d_temp = fModel->GetMarginalized(i);
      if (bch1d_temp) {
         fFlagInfoMarg = true;
         fMean.push_back( bch1d_temp->GetMean() );
         fRMS.push_back( bch1d_temp->GetRMS() );
         fMargMode.push_back( bch1d_temp->GetMode() );
         for (unsigned j = 0; j < fSumProb.size(); ++j)
            fQuantiles.push_back( bch1d_temp->GetQuantile( fSumProb.at(j) ) );
         std::vector<double> intervals = bch1d_temp->GetSmallestIntervals();
         int nintervals = int(intervals.size() / 5);
         fSmallInt.push_back(nintervals);
         fSmallInt.insert( fSmallInt.end(), intervals.begin(), intervals.end() );
      }
      else {
         double tmpval = fModel->GetParameter(i)->GetUpperLimit() - fModel->GetParameter(i)->GetLowerLimit();
         tmpval = fModel->GetParameter(i)->GetLowerLimit() - 2. * tmpval;
         fMean.push_back( tmpval );
         fRMS.push_back( tmpval );
         fMargMode.push_back( tmpval );
         for (unsigned j = 0; j < fSumProb.size(); ++j)
            fQuantiles.push_back( tmpval );
         fSmallInt.push_back( 0 );
      }

      // copy 2D marginal information
      for (unsigned j = 0; j < fModel->GetNParameters(); ++j) {
         if (i == j)
            fCorrCoeff.push_back(1.);
         else {
            BCH2D * bch2d_temp = fModel->GetMarginalized(i, j);
            if ( bch2d_temp ) {
               fFlagInfoMarg = true;
               fCorrCoeff.push_back( bch2d_temp->GetHistogram()->GetCorrelationFactor() );
            }
            else
               fCorrCoeff.push_back( 0. );
         }
      }

      // copy optimization information
      if ((fModel->GetBestFitParameters()).size() > 0) {
         fFlagInfoOpt = true;
         fGlobalMode.push_back ( (fModel->GetBestFitParameters()).at(i) );
      }
   }

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintParameterPlot(const char * filename)
{
   // copy summary data
   if (!CopySummaryData())
      return 0;

   // get number of parameters and quantiles
   int npar = fModel->GetNParameters();
   int nquantiles = int( fSumProb.size() );

   // create histogram
   TH1D * hist_axes = new TH1D(
         TString::Format("hist_axes_par_%d",getNextIndex()),
         ";;Scaled parameter range [a.u.]",npar, -0.5, npar-0.5);
   hist_axes->SetStats(kFALSE);
   for (int i = 0; i < npar; ++i)
      hist_axes->GetXaxis()->SetBinLabel( i+1, fParName.at(i).c_str() );
   hist_axes->GetXaxis()->SetLabelOffset(0.015);
   hist_axes->GetXaxis()->SetLabelSize(0.06);
   hist_axes->GetXaxis()->SetTickLength(0.0);
   hist_axes->GetYaxis()->SetRangeUser(-0.1, 1.1);
   hist_axes->GetYaxis()->SetTickLength(0.0);

   // create graphs
   TGraphErrors * graph_quantiles = new TGraphErrors(npar*nquantiles);
   graph_quantiles->SetMarkerSize(0);
   graph_quantiles->SetLineColor(38);
   graph_quantiles->SetLineStyle(2);

   TGraphErrors * graph_mean = new TGraphErrors(npar);
   graph_mean->SetMarkerColor(kBlack);
   graph_mean->SetMarkerStyle(20);

   TGraphErrors * graph_mode = new TGraphErrors(npar);
   graph_mode->SetMarkerColor(kRed);
   graph_mode->SetMarkerStyle(20);

   TGraphAsymmErrors * graph_intervals = new TGraphAsymmErrors(0);
   graph_intervals->SetFillColor(kYellow);
   graph_intervals->SetLineStyle(2);
   graph_intervals->SetLineColor(kRed);
   graph_intervals->SetMarkerSize(0);

   // fill graphs
   int indexintervals = 0;

   // fill graph quantiles
   if (fFlagInfoMarg) {
      for (int i = 0; i < npar; ++i) {
         for (int j = 0; j < nquantiles; ++j) {
            graph_quantiles->SetPoint(i*nquantiles+j,double(i),
                  (fQuantiles.at(i*nquantiles+j) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
            graph_quantiles->SetPointError(i*nquantiles+j, 0.5, 0.0);
         }
      }
   }

   // fill graph mean and rms
   if (fFlagInfoMarg) {
      for (int i = 0; i < npar; ++i) {
         // fill graph mean
         graph_mean->SetPoint(i, double(i), (fMean.at(i) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
         graph_mean->SetPointError(i, 0.0, fRMS.at(i)/(fParMax.at(i)-fParMin.at(i)));
      }
   }

   // fill graph mode
   if (fFlagInfoOpt) {
      for (int i = 0; i < npar; ++i)
         graph_mode->SetPoint(i, double(i), (fGlobalMode.at(i) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
   }

   // fill graph smallest intervals
   if (fFlagInfoMarg) {
      for (int i = 0; i < npar; ++i) {
         int nintervals = int(fSmallInt.at(indexintervals++));
         for (int j = 0; j < nintervals; ++j) {
            double xmin = fSmallInt.at(indexintervals++);
            double xmax = fSmallInt.at(indexintervals++);
            indexintervals++;
            double xlocalmaxpos = fSmallInt.at(indexintervals++);
            indexintervals++;
            int npoints = graph_intervals->GetN();
            graph_intervals->SetPoint(npoints,double(i),
                  (xlocalmaxpos - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
            graph_intervals->SetPointError(npoints,0.5, 0.5,
                  (xlocalmaxpos - xmin)/(fParMax.at(i)-fParMin.at(i)),
                  (xmax - xlocalmaxpos)/(fParMax.at(i)-fParMin.at(i)));
         }
      }
   }

   // create legend
   TLegend * legend = new TLegend(0.15, 0.88, 0.85, 0.99);
   legend->SetBorderSize(0);
   legend->SetFillColor(0);

   // create latex
   TLatex * latex = new TLatex();
   latex->SetTextSize(0.02);

   // create lines
   TLine * line_top = new TLine(-0.5, 1.0, npar-0.5, 1.0);
   line_top->SetLineColor(kBlack);
   line_top->SetLineStyle(1);
   line_top->SetLineWidth(2);

   TLine * line_bot = new TLine(-0.5, 0.0, npar-0.5, 0.0);
   line_bot->SetLineColor(kBlack);
   line_bot->SetLineStyle(1);
   line_bot->SetLineWidth(2);

   // print to file
   TCanvas * c_par = new TCanvas(TString::Format("c_par_%d",getNextIndex()));
   c_par->cd();
   hist_axes->Draw();
   line_top->Draw();
   line_bot->Draw();
   if (fFlagInfoMarg) {
      graph_intervals->DrawClone("SAME2");
      for (int i = 0; i < graph_intervals->GetN(); ++i)
         graph_intervals->SetPointError(i, 0.5, 0.5, 0.0, 0.0);
      graph_intervals->Draw("SAMEPZ");
      graph_quantiles->Draw("SAMEPZ");
      graph_mean->Draw("SAMEP");
      legend->AddEntry(graph_quantiles, "Quantiles (5%, 10%, 16%, 50%, 84%, 90%, 95%)", "L");
      legend->AddEntry(graph_mean,      "Mean and RMS", "LEP");
      legend->AddEntry(graph_intervals, "Smallest 68% intervals and local modes", "FL");
   }
   if (fFlagInfoOpt) {
      graph_mode->Draw("SAMEP");
      legend->AddEntry(graph_mode,      "Global mode", "P");
   }
   for (int i = 0; i < npar;++i) {
      //      latex->DrawLatex(double(i)-0.1, 0.010, Form("%+3.3f", fParMin.at(i)));
      //      latex->DrawLatex(double(i)-0.1, 0.965, Form("%+3.3f", fParMax.at(i)));
      latex->DrawLatex(double(i)-0.1, 0.010-0.07, Form("%+3.3g", fParMin.at(i)));
      latex->DrawLatex(double(i)-0.1, 0.965+0.07, Form("%+3.3g", fParMax.at(i)));
   }
   latex->SetNDC();
   latex->DrawLatex(0.9, 0.175, "Par. min.");
   latex->DrawLatex(0.9, 0.83, "Par. max.");
   legend->Draw("SAME");
   gPad->RedrawAxis();
   c_par->Print(filename);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintCorrelationMatrix(const char * filename)
{
   // copy summary data
   if (!CopySummaryData())
      return 0;

   // check if marginalized information is there
   if (!fFlagInfoMarg)
      return 0;

   // get number of parameters
   int npar = fModel->GetNParameters();

   // create histogram
   TH2D * hist_corr = new TH2D(
         TString::Format("hist_corr_%d",getNextIndex()),
         ";;",npar, -0.5, npar-0.5,npar, -0.5, npar-0.5);
   hist_corr->SetStats(kFALSE);
   hist_corr->GetXaxis()->SetTickLength(0.0);
   hist_corr->GetYaxis()->SetTickLength(0.0);
   hist_corr->GetZaxis()->SetRangeUser(-1.0, 1.0);

   for (int i = 0; i < npar; ++i) {
      hist_corr->GetXaxis()->SetLabelSize(0.06);
      hist_corr->GetYaxis()->SetLabelSize(0.06);
      if (npar < 5) {
    	 hist_corr->GetXaxis()->SetBinLabel( i+1, fParName.at(i).c_str() );
         hist_corr->GetYaxis()->SetBinLabel( npar-i, fParName.at(i).c_str() );
      }
      else {
    	 hist_corr->GetXaxis()->SetBinLabel( i+1, TString::Format("%d",i) );
         hist_corr->GetYaxis()->SetBinLabel( npar-i, TString::Format("%d",i) );
      }
   }

   // fill plot
   for (int i = 0; i < npar; ++i)
      for (int j = 0; j < npar; ++j) {
         int index = i * npar + j;
         double corr = fCorrCoeff.at(index);
         hist_corr->SetBinContent(i+1, npar-j, corr);
      }

   // print to file
   TCanvas * c_corr = new TCanvas(TString::Format("c_corr_matrix_%d",getNextIndex()));
   c_corr->cd();
   hist_corr->Draw("colz text");

   TF1 * f = new TF1("fUp","x",-0.5,npar-0.5);
   TGaxis * A1 = new TGaxis(-0.5,npar-0.5,npar-0.5,npar-0.5,"fUp",100,"-");
   A1->ImportAxisAttributes(hist_corr->GetXaxis());
   A1->Draw();

   // redraw the histogram to overlay thetop axis tick marks since
   // we don't know how to make them disappear
   hist_corr->GetXaxis()->SetLabelSize(0.);
   hist_corr->Draw("colz text same");

   for (int i = 0; i < npar; ++i)
      for (int j = 0; j < npar; ++j) {
         BCH2D * bch2d_temp = fModel->GetMarginalized(fModel->GetParameter(i),fModel->GetParameter(j));
         if ( bch2d_temp || i==j )
            continue;

         TBox * bempty = new TBox(
            hist_corr->GetXaxis()->GetBinLowEdge(i+1),
            hist_corr->GetYaxis()->GetBinLowEdge(npar-j),
            hist_corr->GetXaxis()->GetBinLowEdge(i+2),
            hist_corr->GetYaxis()->GetBinLowEdge(npar-j+1)
         );
         bempty->SetLineStyle(0);
         bempty->SetLineWidth(0);
         bempty->SetFillColor(kWhite);
         bempty->Draw();
      }

   // redraw top and right axes
   TLine * lA1 = new TLine(-0.5,npar-0.5,npar-0.5,npar-0.5);
   lA1->Draw("same");
   TLine * lA2 = new TLine(npar-0.5,npar-0.5,npar-0.5,-0.5);
   lA2->Draw("same");

   gPad->RedrawAxis();
   c_corr->Print(filename);

   delete f;
   delete A1;
   delete lA1;
   delete lA2;
   delete hist_corr;
   delete c_corr;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintCorrelationPlot(const char * filename)
{
   // get number of parameters
   int npar = fModel->GetNParameters();

   TCanvas * c = new TCanvas(TString::Format("c_corr_%d",getNextIndex()));
   c->cd();

   double margin = .1;
   double padsize = (1. - 2.*margin) / (double)npar;

   // array with pads holding the histograms
   TPad ** pad = new TPad*[npar*npar];

   // position of pads
   double xlow, xup, ylow, yup;
   double marginleft, marginright, margintop, marginbottom;
   marginleft=marginright=margintop=marginbottom=0.01;

   // rotate label to avoid overlap with long parameter names
   size_t maxlength = 0;
   for (int i = 0; i < npar ; ++i) {
      maxlength = std::max(maxlength, fModel->GetParameter(i)->GetName().length());
   }

   double rotation = 0;
   short xalignment = 22;
   short yalignment = 22;

   if (maxlength > 20) {
      rotation = 10;
      xalignment = 12;
      yalignment = 32;
   }

   // drawing all histograms
   for (int i=npar-1;i>=0;--i) {
      xlow = (double)i*padsize + margin;
      xup = xlow+padsize;

      for (int j=i;j<=npar-1;j++) {
         yup = 1. - (double)j*padsize - margin;
         ylow = yup-padsize;

         // preparing the pad
         int ipad=i*npar+j;
         pad[ipad]=new TPad(TString::Format("pad_%d_%d",ipad,getNextIndex()), "", xlow, ylow, xup, yup);
         pad[ipad]->SetTopMargin(margintop);
         pad[ipad]->SetBottomMargin(marginbottom);
         pad[ipad]->SetLeftMargin(marginleft);
         pad[ipad]->SetRightMargin(marginright);
         pad[ipad]->SetFillColor(kWhite);

         pad[ipad]->Draw();
         pad[ipad]->cd();

         // get the histogram
         TH1 * hh = 0;
         BCH1D * bh1 = 0;
         BCH2D * bh2 = 0;
         if(i==j) {
           bh1 = fModel->GetMarginalized(fModel->GetParameter(i));
           if (bh1)
             hh = (TH1D*)bh1->GetHistogram()->Clone();
         }
         else {
           bh2 = fModel->GetMarginalized(fModel->GetParameter(i),fModel->GetParameter(j));
           if (bh2)
             hh = (TH2D*)bh2->GetHistogram()->Clone();
         }

         // if the histogram is not available, draw N/A
         if (!hh) {
            pad[ipad]->SetFillColor(kGray);
            TBox * box = new TBox(marginleft,marginbottom,1.-marginright,1.-margintop);
            box->SetLineWidth(1);
            box->SetLineColor(kGray+1);
            box->SetFillColor(kWhite);
            box->Draw();

            TText * text_na = new TText(.5,.5,"N/A");
            text_na->SetTextFont(42);
            text_na->SetTextAlign(22);
            text_na->SetTextSize(.8/(double)npar);
            text_na->SetTextColor(kGray+1);
            text_na->Draw();
         }
         // otherwise draw the histogram
         else {
            if(i==j) {
              bh1->Draw("BTsiB3CS1D0");
            }
            else {
              bh2->Draw("BTfB3CS1nL");
            }

            hh->GetXaxis()->SetLabelOffset(5500);
            hh->GetYaxis()->SetLabelOffset(5500);
            hh->GetXaxis()->SetTitleSize(10.00);
            hh->GetYaxis()->SetTitleSize(10.00);
         }

         c->cd();

         // draw axis label
         double labelsize = .8/(double)npar/10.;
         double xtext, ytext;

         // y axis
         if(i==0) {
            TLatex * label = new TLatex();
            label->SetTextFont(62);
            label->SetTextSize(labelsize);
            label->SetTextAlign(yalignment);
            label->SetNDC();

            label->SetTextAngle(90 - rotation);

            xtext = margin * (1. - 8. * labelsize);
            ytext = yup - padsize / 2.;

            label->DrawLatex(xtext,ytext,fModel->GetParameter(j)->GetLatexName().c_str());
         }

         // x axis
         if(j==npar-1) {
            TLatex * label = new TLatex();
            label->SetTextFont(62);
            label->SetTextSize(labelsize);
            label->SetTextAlign(xalignment);
            label->SetNDC();

            label->SetTextAngle(360 - rotation);

            xtext = xlow + padsize / 2.;
            ytext = margin * (1. - 6. * labelsize);

            label->DrawLatex(xtext,ytext, fModel->GetParameter(i)->GetLatexName().c_str());
         }
      }
   }

   gPad->RedrawAxis();
   c->Print(filename);

   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintKnowledgeUpdatePlot1D(int index, const char * filename, std::string options_post, std::string options_prior)
{
   // perform analysis
   CalculatePriorModel();

   // create canvas
   TCanvas * c = new TCanvas();
   c->cd();

   // draw
   DrawKnowledgeUpdatePlot1D(index, options_post, options_prior);

   // print
   c->Print(filename);

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::DrawKnowledgeUpdatePlot1D(int index, std::string options_post, std::string options_prior)
{
   // option flags
   bool flag_slice_post = false;
   bool flag_slice_prior = false;

   // check content of options string
   if (options_post.find("slice") < options_post.size()) {
      flag_slice_post = true;
   }
   if (options_prior.find("slice") < options_prior.size()) {
      flag_slice_prior = true;
   }

   // prepare legend
   TLegend* legend = new TLegend();
   legend->SetBorderSize(0);
   legend->SetFillColor(kWhite);
   legend->SetTextAlign(12);
   legend->SetTextFont(62);
   legend->SetTextSize(0.03);

   // get histograms;
   const BCParameter * par = fModel->GetParameter(index);
   BCH1D* hist_prior = fPriorModel->GetMarginalized(par);
   BCH1D* hist_posterior = 0;

   if (flag_slice_prior && fPriorModel->GetNParameters()==2) {
      if (index == 0) {
         TH1D* hist = fPriorModel->GetSlice(fPriorModel->GetParameter(0),fPriorModel->GetParameter(1))->GetHistogram()->ProjectionX(Form("projx_%i",BCLog::GetHIndex()));
         hist->Scale(1.0/hist->Integral("width"));
         for (int i = 1; i <= hist_prior->GetHistogram()->GetNbinsX(); ++i)
            hist_prior->GetHistogram()->SetBinContent(i, hist->GetBinContent(i));
      }
      else {
         TH1D* hist = fPriorModel->GetSlice(fPriorModel->GetParameter(0),fPriorModel->GetParameter(1))->GetHistogram()->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
         hist->Scale(1.0/hist->Integral("width"));
         for (int i = 1; i <= hist_prior->GetHistogram()->GetNbinsX(); ++i)
            hist_prior->GetHistogram()->SetBinContent(i, hist->GetBinContent(i));
      }
      hist_prior->GetHistogram()->SetStats(kFALSE);
   }
   else if (flag_slice_prior && fPriorModel->GetNParameters()==1) {
      hist_prior = fPriorModel->GetSlice(par);
      hist_prior->GetHistogram()->SetStats(kFALSE);
   }
   // if marginal doesn't exist, skip ahead
   if ( !hist_prior)
     return 0;

   hist_prior->GetHistogram()->SetLineColor(kRed);

   hist_posterior = fModel->GetMarginalized(par);
   if (flag_slice_post && fModel->GetNParameters()==2) {
      if (index == 0) {
         TH1D* hist = fModel->GetSlice(fModel->GetParameter(0),fModel->GetParameter(1))->GetHistogram()->ProjectionX(Form("projx_%i",BCLog::GetHIndex()));
         hist->Scale(1.0/hist->Integral("width"));
         for (int i = 1; i <= hist_posterior->GetHistogram()->GetNbinsX(); ++i)
            hist_posterior->GetHistogram()->SetBinContent(i, hist->GetBinContent(i));
      }
      else {
         TH1D* hist = fModel->GetSlice(fModel->GetParameter(0),fModel->GetParameter(1))->GetHistogram()->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
         hist->Scale(1.0/hist->Integral("width"));
         for (int i = 1; i <= hist_posterior->GetHistogram()->GetNbinsX(); ++i)
            hist_posterior->GetHistogram()->SetBinContent(i, hist->GetBinContent(i));
      }
      hist_posterior->GetHistogram()->SetStats(kFALSE);
   }
   else if (flag_slice_post && fModel->GetNParameters()==1) {
      hist_posterior = fModel->GetSlice(par);
      hist_posterior->GetHistogram()->SetStats(kFALSE);
   }

   // if marginal doesn't exist, skip ahead
   if ( !hist_posterior)
     return 0;

   legend->AddEntry(hist_prior->GetHistogram(), "prior", "L");
   legend->AddEntry(hist_posterior->GetHistogram(), "posterior", "L");

   // scale histograms
   hist_posterior->GetHistogram()->Scale(1./hist_posterior->GetHistogram()->Integral("width"));
   hist_prior->GetHistogram()->Scale(1.0/hist_prior->GetHistogram()->Integral("width"));

   // get maximum
   double max_prior = hist_prior->GetHistogram()->GetMaximum();
   double max_posterior = hist_posterior->GetHistogram()->GetMaximum();
   double maxy = 1.1 * TMath::Max(max_prior, max_posterior);

   double height = 0.03*legend->GetNRows();

   // plot
   hist_prior->GetHistogram()->GetXaxis()->SetNdivisions(508);
   hist_posterior->GetHistogram()->GetXaxis()->SetNdivisions(508);

   hist_prior->Draw(options_prior);
   hist_posterior->Draw(std::string(options_post+"same").c_str());

   // scale axes
   hist_prior->GetHistogram()->GetYaxis()->SetRangeUser(0.0, maxy);
   hist_posterior->GetHistogram()->GetYaxis()->SetRangeUser(0.0, maxy);

   gPad->SetTopMargin(0.02);
   double xlegend1 = gPad->GetLeftMargin() + 0.10 * (1.0 - gPad->GetRightMargin() - gPad->GetLeftMargin());
   double xlegend2 = 1.0-gPad->GetRightMargin();
   double ylegend1 = 1.-gPad->GetTopMargin()-height;
   double ylegend2 = 1.-gPad->GetTopMargin();

   // place legend on top of histogram
   legend->SetX1NDC(xlegend1);
   legend->SetX2NDC(xlegend2);
   legend->SetY1NDC(ylegend1);
   legend->SetY2NDC(ylegend2);

   // draw legend
   legend->Draw();

   // rescale top margin
   gPad->SetTopMargin(1.-ylegend1+0.01);

   gPad->RedrawAxis();

   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintKnowledgeUpdatePlots(const char * filename, std::string options)
{
   // perform analysis
   CalculatePriorModel();

   // option flags
   bool flag_slice = false;

   // check content of options string
   if (options.find("slice") < options.size()) {
      flag_slice = true;
   }

   std::string file(filename);

   // check if file extension is pdf
   if ( (file.find_last_of(".") != std::string::npos) &&
         (file.substr(file.find_last_of(".")+1) == "pdf") ) {
      ; // it's a PDF file

   }
   else if ( (file.find_last_of(".") != std::string::npos) &&
         (file.substr(file.find_last_of(".")+1) == "ps") ) {
      ; // it's a PS file
   }
   else {
      ; // make it a PDF file
      file += ".pdf";
   }

   // create canvas and prepare postscript
   TCanvas * c = new TCanvas(TString::Format("c_%d",getNextIndex()));
   c->cd();

   // loop over all parameters and draw 1D plots
   int npar = fModel->GetNParameters();
   c->Print(std::string(file + "[").c_str());
   for (int i = 0; i < npar; ++i) {
      if ( !DrawKnowledgeUpdatePlot1D(i, options, options))
         continue;
      c->Print(file.c_str());
   }

   // create legend
   TLegend * legend2d = new TLegend();
   legend2d->SetBorderSize(0);
   legend2d->SetFillColor(0);
   legend2d->SetTextAlign(12);
   legend2d->SetTextFont(62);
   legend2d->SetTextSize(0.03);

  // create markers and arrows
   TMarker * marker_prior = new TMarker();
   marker_prior->SetMarkerStyle(24);
   marker_prior->SetMarkerColor(kRed);

   TMarker * marker_posterior = new TMarker();
   marker_posterior->SetMarkerStyle(24);
   marker_posterior->SetMarkerColor(kBlack);

   TArrow * arrow = new TArrow();
   arrow->SetArrowSize(0.02);
   arrow->SetLineColor(kBlue);
	 //   arrow->SetLineStyle(2);

   // loop over all parameters
   for (int i = 0; i < npar; ++i) {
      for (int j = 0; j < i; ++j) {
         c->cd();

         // get parameters
         const BCParameter * par1 = fModel->GetParameter(i);
         const BCParameter * par2 = fModel->GetParameter(j);

         // get 2-d histograms
         BCH2D* bch2d_2dprior = 0;
         BCH2D* bch2d_2dposterior = 0;
         if (flag_slice && npar == 2) {
            bch2d_2dprior = fPriorModel->GetSlice(par2, par1);
            bch2d_2dposterior = fModel->GetSlice(par2, par1);
         }
         else {
            bch2d_2dprior = fPriorModel->GetMarginalized(par1, par2);
            bch2d_2dposterior = fModel->GetMarginalized(par1, par2);
         }

         // can't draw anything
         if ( !bch2d_2dprior || !bch2d_2dposterior)
            continue;

         // get histograms
         TH2D* hist_2dprior = bch2d_2dprior->GetHistogram();
         hist_2dprior->SetLineColor(kRed);
         TH2D* hist_2dposterior = bch2d_2dposterior->GetHistogram();

         // scale histograms
         hist_2dprior->Scale(1.0/hist_2dprior->Integral("width"));
         hist_2dposterior->Scale(1.0/hist_2dposterior->Integral("width"));

         // calculate contours
         bch2d_2dprior->CalculateIntegratedHistogram();
         bch2d_2dposterior->CalculateIntegratedHistogram();

         double level[1] = {bch2d_2dprior->GetLevel(0.32)};
         hist_2dprior->SetContour(1, level);
				 hist_2dprior->Draw("CONT3");
         level[0] = bch2d_2dposterior->GetLevel(0.32);
         hist_2dposterior->SetContour(1, level);
				 hist_2dposterior->Draw("CONT3 SAME");

         std::vector<double> mode_prior = fPriorModel->GetBestFitParameters();
         std::vector<double> mode_posterior = fModel->GetBestFitParameters();

         marker_prior->DrawMarker(mode_prior.at(j), mode_prior.at(i));
         marker_posterior->DrawMarker(mode_posterior.at(j), mode_posterior.at(i));
         arrow->DrawArrow(mode_prior.at(j), mode_prior.at(i), mode_posterior.at(j), mode_posterior.at(i));

         if (i==1 && j == 0) {
            legend2d->AddEntry(hist_2dprior, "smallest 68% interval(s) of prior", "L");
            legend2d->AddEntry(hist_2dposterior, "smallest 68% interval(s) of posterior", "L");
            legend2d->AddEntry(marker_prior, "prior mode", "P");
            legend2d->AddEntry(marker_posterior, "posterior mode", "P");
            legend2d->AddEntry(arrow, "change in mode", "L");
         }
         gPad->SetTopMargin(0.02);

         double height = 0.03*legend2d->GetNRows();

         double xlegend1 = gPad->GetLeftMargin();
         double xlegend2 = 1.0-gPad->GetRightMargin();
         double ylegend1 = 1.-gPad->GetTopMargin()-height;
         double ylegend2 = 1.-gPad->GetTopMargin();

         // place legend on top of histogram
         legend2d->SetX1NDC(xlegend1);
         legend2d->SetX2NDC(xlegend2);
         legend2d->SetY1NDC(ylegend1);
         legend2d->SetY2NDC(ylegend2);

         legend2d->Draw();

         gPad->SetTopMargin(1.-ylegend1+0.01);

         gPad->RedrawAxis();
         c->Print(file.c_str());
      }
   }

   // close output
   c->Print(std::string(file + "]").c_str());
   c->Update();

   // free memory
   delete legend2d;
   delete marker_prior;
   delete marker_posterior;
   delete arrow;
   delete c;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintParameterLatex(const char * filename)
{
   // open file
   std::ofstream ofi(filename);
   ofi.precision(3);

   // check if file is open
   if(!ofi.is_open()) {
      std::cerr << "Couldn't open file " << filename <<std::endl;
      return 0;
   }

   // get number of parameters and quantiles
   int npar = fModel->GetNParameters();

   // print table
   ofi
      << "\\documentclass[11pt, a4paper]{article}" << std::endl
      << std::endl
      << "\\begin{document}" << std::endl
      << std::endl
      << "\\begin{table}[ht!]" << std::endl
      << "\\begin{center}" << std::endl
      <<"\\begin{tabular}{llllllll}" << std::endl
      << "\\hline" << std::endl
      << "Parameter & Mean & RMS & Gl. mode & Mode & Median & 16\\% quant. & 84\\% quant. \\\\" << std::endl
      << "\\hline" << std::endl;

   for (int i = 0; i < npar; ++i) {
      const BCParameter * par = fModel->GetParameter(i);
      BCH1D * bch1d = fModel->GetMarginalized(par);
      ofi
         << par->GetName() << " & "
         << bch1d->GetMean() << " & "
         << bch1d->GetRMS() << " & "
         << fModel->GetBestFitParameters().at(i) << " & "
         << bch1d->GetMode() << " & "
         << bch1d->GetMedian() << " & "
         << bch1d->GetQuantile(0.16) << " & "
         << bch1d->GetQuantile(0.84) << " \\\\" << std::endl;
   }
   ofi
      << "\\hline" << std::endl
      << "\\end{tabular}" << std::endl
      << "\\caption{Summary of the parameter estimates.}" << std::endl
      << "\\end{center}" << std::endl
      << "\\end{table}" << std::endl
      << std::endl
      << "\\end{document}" << std::endl;

   // close file
   ofi.close();

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::CalculatePriorModel()
{
   // create new prior model
   delete fPriorModel;

   fPriorModel = new BCSummaryPriorModel();

   // set model
   fPriorModel->SetModel(fModel);

   // perform marginalization
   fPriorModel->MarginalizeAll();

   // perform minimization
   fPriorModel->FindMode( fPriorModel->GetBestFitParameters() );

   // no error
   return 1;
}
