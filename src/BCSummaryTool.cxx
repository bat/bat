/*
 * Copyright (C) 2007-2013, the BAT core developer team
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
#include <TColor.h>
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
int BCSummaryTool::PrintParameterPlot(const char * filename, int npar) {
	TCanvas * c_par = new TCanvas(TString::Format("c_par_%d",getNextIndex()));
	c_par -> Print(Form("%s[",filename));
	int return_val = 1;
	if (npar<=0) {
		return_val *= PrintParameterPlot(0,fModel->GetNParameters()-1, filename);
		return_val *= PrintParameterPlot(fModel->GetNParameters()-1,fModel->GetNObservables(), filename);
	} else {
		for (unsigned i = 0; i<fModel->GetNParameters(); i += npar)
			return_val *= PrintParameterPlot(i,std::min(i+npar-1,fModel->GetNParameters()-1), filename);
		for (unsigned i = fModel->GetNParameters(); i<fModel->GetNObservables(); i += npar)
			return_val *= PrintParameterPlot(i,std::min(i+npar-1,fModel->GetNObservables()-1), filename);
	}
	c_par -> Print(Form("%s]",filename));
	return return_val;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintParameterPlot(unsigned i0, int i1, const char * filename) {

	unsigned i_1 = (i1>=0 && i1<(int)fModel->GetNObservables()) ? i1 : fModel->GetNObservables()-1;

	if (i_1 < i0) {
		BCLog::OutError(Form("BCSummaryTool::PrintParameterPlot : invalid parameter range [%d, %d]",i0,i_1));
		return 0;
	}

	/////////////////////////
	// Gather information

	double interval_content = 68e-2; // 68% credibility interval
	
	std::vector<double> x_quantiles;
	std::vector<double> quantiles;
	std::vector<double> x_i;
	std::vector<double> x_i_bf;
	std::vector<double> mean;
	std::vector<double> rms;
	std::vector<double> global_mode;
	std::vector<double> local_mode;
	std::vector<double> interval_lo;
	std::vector<double> interval_hi;
	
	for (unsigned i = i0; i <= i_1; ++i) {
		
		// Global Mode
		x_i_bf.push_back(i);
		global_mode.push_back(fModel->GetObservable(i)->PositionInRange(fModel->GetBestFitParameter(i)));

		if (!fModel->MarginalizedHistogramExists(i))
			continue;

		BCH1D * bch1d_temp = fModel -> GetMarginalized(i);
		if (!bch1d_temp)
			continue;

		x_i.push_back(i);

		// quantiles
		x_quantiles.insert(x_quantiles.end(),fSumProb.size(),i);
		for (unsigned j = 0; j < fSumProb.size(); ++j)
			quantiles.push_back(fModel->GetObservable(i)->PositionInRange(bch1d_temp->GetQuantile(fSumProb[j])));

		// mean
		mean.push_back(fModel->GetObservable(i)->PositionInRange(bch1d_temp->GetMean()));
		rms.push_back(bch1d_temp->GetRMS()/fModel->GetObservable(i)->GetRangeWidth());
		 
		// Local Mode
		local_mode.push_back(fModel->GetObservable(i)->PositionInRange(bch1d_temp->GetMode()));

		// smallest interval
		std::vector<double> intervals = bch1d_temp->GetSmallestIntervals(interval_content);
		if (intervals.size()>3) {
			interval_lo.push_back(fabs(intervals[3]-intervals[0])/fModel->GetObservable(i)->GetRangeWidth());
			interval_hi.push_back(fabs(intervals[3]-intervals[1])/fModel->GetObservable(i)->GetRangeWidth());
		} else {
			interval_lo.push_back(0);
			interval_hi.push_back(0);
		}
	}

	if (x_i.empty() and x_i_bf.empty())
		return 0;

	/////////////////////////
	// Draw it all

	TCanvas * c_par = new TCanvas(TString::Format("c_par_%d",getNextIndex()));
	c_par -> cd();

	// Create, label, and draw axes
	TH2D * hist_axes = new TH2D(TString::Format("hist_axes_par_%d",getNextIndex()), ";;Scaled parameter range [a.u.]",
															i_1-i0+1, i0-0.5, i_1+0.5, 10, -0.1, 1.1);
	hist_axes -> SetStats(kFALSE);
	hist_axes -> GetXaxis() -> SetLabelOffset(0.015);
	hist_axes -> GetXaxis() -> SetLabelSize(0.06);
	hist_axes -> GetXaxis() -> SetTickLength(0.0);
	// set bin labels
	for (int i=0; i<hist_axes->GetNbinsX(); ++i)
		hist_axes -> GetXaxis() -> SetBinLabel(i+1, fModel->GetObservable(i0+i)->GetLatexName().c_str());
	hist_axes->Draw();

	// Draw lines
	TLine * line = new TLine();
	line -> SetLineColor(kBlack);
	line -> SetLineStyle(1);
	line -> SetLineWidth(2);
	line -> DrawLine(hist_axes->GetXaxis()->GetXmin(), 0.0, hist_axes->GetXaxis()->GetXmax(), 0.0);
	line -> DrawLine(hist_axes->GetXaxis()->GetXmin(), 1.0, hist_axes->GetXaxis()->GetXmax(), 1.0);
	// Mark parameter ranges
	TLatex * latex = new TLatex();
	latex -> SetTextSize(0.02);
	latex -> SetTextAlign(11);
	latex -> DrawLatex(hist_axes->GetXaxis()->GetXmax(),  1.03, "  Par. max.");
	latex -> SetTextAlign(13);
	latex -> DrawLatex(hist_axes->GetXaxis()->GetXmax(), -0.03, "  Par. min.");
	latex -> SetTextAlign(21);
	for (unsigned i = 0; i < fModel->GetNObservables(); ++i) {
		latex -> SetTextAlign(21);
		latex->DrawLatex((double)i,  1.03, Form("%+.*g", fModel->GetObservable(i)->GetPrecision(),fModel->GetObservable(i)->GetUpperLimit()));
		latex -> SetTextAlign(23);
		latex->DrawLatex((double)i, -0.03, Form("%+.*g", fModel->GetObservable(i)->GetPrecision(),fModel->GetObservable(i)->GetLowerLimit()));
	}

	// create legend
	TLegend * legend = new TLegend(0.1, 0.91, 0.9, 0.99);
	legend -> SetBorderSize(0);
	legend -> SetFillColor(0);
	legend -> SetNColumns(2);

	if (!x_i.empty()) {

		// Smallest Interval
		std::vector<double> x_i_err(x_i.size(),0.5);
		TGraphAsymmErrors * graph_intervals = new TGraphAsymmErrors(x_i.size(), x_i.data(), local_mode.data(), x_i_err.data(), x_i_err.data(), interval_lo.data(), interval_hi.data());
		graph_intervals->SetFillColor(kYellow);
		graph_intervals->SetLineStyle(2);
		graph_intervals->SetLineColor(kRed);
		graph_intervals->SetMarkerSize(0);
		graph_intervals->DrawClone("SAME2"); // draw area
		//set y-error zero, to draw line at local mode
		for (int i = 0; i < graph_intervals->GetN(); ++i)
			graph_intervals->SetPointError(i, 0.5, 0.5, 0.0, 0.0);
		graph_intervals->Draw("SAMEZ"); // draw local mode

		// Quantiles graph
		if (!fSumProb.empty()) {
			std::vector<double> quantiles_err(x_quantiles.size(),0.5);
			TGraphErrors * graph_quantiles = new TGraphErrors(x_quantiles.size(), x_quantiles.data(), quantiles.data(), quantiles_err.data(), 0);
			graph_quantiles->SetMarkerSize(0);
			graph_quantiles->SetLineColor(38);
			graph_quantiles->SetLineStyle(2);
			graph_quantiles->Draw("SAMEZ");
			std::string quantiles_text = "Quantiles (";
			for (unsigned i=0; i<fSumProb.size()-1; ++i)
				quantiles_text += Form("%.0f%%, ",fSumProb[i]*100);
			quantiles_text += (fSumProb.size()>0) ? Form("%.0f%%)",fSumProb.back()) : "none)";
			legend -> AddEntry(graph_quantiles, quantiles_text.c_str(), "L");
		}

		// Means & RMSs
		TGraphErrors * graph_mean = new TGraphErrors(x_i.size(), x_i.data(), mean.data(), 0, rms.data());
		graph_mean->SetMarkerColor(kBlack);
		graph_mean->SetMarkerStyle(20);
		graph_mean->Draw("SAMEP");

		legend -> AddEntry(graph_mean, "Mean and RMS", "LEP");
		legend -> AddEntry(graph_intervals, "Smallest 68% interval and local mode", "FL");

	}

	// Global Modes
	if (!x_i_bf.empty()) {
		TGraph * graph_mode = new TGraph(x_i_bf.size(), x_i_bf.data(), global_mode.data());
		graph_mode->SetMarkerColor(kRed);
		graph_mode->SetMarkerStyle(20);
		graph_mode->Draw("SAMEP");
		legend->AddEntry(graph_mode, "Global mode", "P");
	}

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
int BCSummaryTool::PrintCorrelationPlot(const char * filename) {

	// Array of indices for which any maginalizations were stored
	std::vector<unsigned> I;
	for (unsigned i = 0; i < fModel->GetN1DMarginalizations(); ++i)
		if (fModel->MarginalizedHistogramExists(i))
			I.push_back(i);
		else 
			for (unsigned j = i+1; j < fModel->GetN2DMarginalizations(); ++j)
				if (fModel->MarginalizedHistogramExists(i,j))
					I.push_back(i);
	
	if (I.size() == 0)
		return 0;

	TCanvas * c = new TCanvas(TString::Format("c_corr_%d",getNextIndex()));
	c->cd();
	
	double margin = 0.1;
	double padsize = (1 - 2*margin) / I.size();

	// array with pads holding the histograms
	std::vector<std::vector<TPad*> > pad (I.size(), std::vector<TPad*>(I.size(),0));
	
	// position of pads
	double xlow, xup, ylow, yup;
	double marginleft   = 0.01;
	double marginright  = 0.01;
	double margintop    = 0.01;
	double marginbottom = 0.01;
	
	TLatex * ylabel = new TLatex();
	ylabel->SetTextFont(62);
	ylabel->SetTextSize(8e-2/I.size());
	ylabel->SetTextAlign(22);			// TODO: set to 32, if latex names too long
	ylabel->SetNDC();
	ylabel->SetTextAngle(90);			// TODO: set to 80, if latex names too long
	
	TLatex * xlabel = new TLatex();
	xlabel->SetTextFont(62);
	xlabel->SetTextSize(8e-2/I.size());
	xlabel->SetTextAlign(22);			// TODO: set to 12, if latex names too long
	xlabel->SetNDC();
	xlabel->SetTextAngle(0);			// TODO: set to 350, if latex names too long

	// Box + Text for empty squares:
	TBox * box_na = new TBox();
	box_na -> SetLineWidth(1);
	box_na -> SetLineColor(kGray+1);
	box_na -> SetFillColor(kWhite);
	TText * text_na = new TText();
	text_na -> SetTextFont(42);
	text_na -> SetTextAlign(22);
	text_na -> SetTextSize(8e-1/I.size());
	text_na -> SetTextColor(kGray+1);

	// drawing all histograms
	for (unsigned i = 0; i < I.size(); ++i) {
		xlow = i*padsize + margin;
		xup = xlow + padsize;

		for (unsigned j = i; j < I.size(); ++j) {
			yup = 1. - j*padsize - margin;
			ylow = yup - padsize;

			// preparing the pad
			pad[i][j] =  new TPad(TString::Format("pad_%d_%d_%d",i,j,getNextIndex()), "", xlow, ylow, xup, yup);
			pad[i][j] -> SetMargin(marginleft,marginright,marginbottom,margintop);
			pad[i][j] -> SetFillColor(kWhite);
			pad[i][j] -> Draw();
			pad[i][j] -> cd();

			// get the histogram
			TH1 * hh = 0;
			BCH1D * bh1 = 0;
			BCH2D * bh2 = 0;
			
			if (i==j) {
				bh1 = fModel->GetMarginalized(I[i]);
				hh = bh1 -> GetHistogram();
			}
			else {
				bh2 = fModel->GetMarginalized(I[i],I[j]);
				hh = bh2 -> GetHistogram();
			}
			
			if (!bh1 and !bh2) { // if the histogram is not available, draw N/A

				pad[i][j] -> SetFillColor(kGray);
				box_na -> DrawBox(marginleft,marginbottom,1.-marginright,1.-margintop);
				text_na -> DrawText(.5,.5,"N/A");

			}	else {									// otherwise draw the histogram

				if (bh1)
					bh1->Draw("BTsiB3CS1D0");
				else
					bh2->Draw("BTfB3CS1nL");

				hh->GetXaxis()->SetLabelOffset(5500);
				hh->GetYaxis()->SetLabelOffset(5500);
				hh->GetXaxis()->SetTitleSize(10.00);
				hh->GetYaxis()->SetTitleSize(10.00);
				
				c->cd();

				// y axis
				if(i==0) {
					if (I[j] < fModel->GetNParameters())
						ylabel -> DrawLatex(margin*(1-8*ylabel->GetTextSize()), yup-padsize/2., fModel->GetParameter(I[j])->GetLatexName().c_str());
					else
						ylabel -> DrawLatex(margin*(1-8*ylabel->GetTextSize()), yup-padsize/2., fModel->GetUserObservable(I[j]-fModel->GetNParameters())->GetLatexName().c_str());
				}
				
				// x axis
				if(j==I.size()-1) {
					if (I[i] < fModel->GetNParameters())
						xlabel -> DrawLatex(xlow+padsize/2., margin*(1-6*xlabel->GetTextSize()), fModel->GetParameter(I[i])->GetLatexName().c_str());
					else 
						xlabel -> DrawLatex(xlow+padsize/2., margin*(1-6*xlabel->GetTextSize()), fModel->GetUserObservable(I[i]-fModel->GetNParameters())->GetLatexName().c_str());
				}
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
int BCSummaryTool::DrawKnowledgeUpdatePlot1D(unsigned index, std::string options_post, std::string options_prior) {
	// option flags
	bool flag_slice_post  = (options_post.find("slice") < options_post.size());
	bool flag_slice_prior = (options_prior.find("slice") < options_prior.size());

	// Get Prior
	BCH1D * bch1d_prior = 0;
	TLine * const_prior = (fModel->IsPriorConstant(index)) ? new TLine() : 0;
	TF1   * f1_prior    = (const_prior) ? 0 : dynamic_cast<TF1*> (fModel->PriorContainer(index));
	TH1   * h1_prior    = (const_prior) ? 0 : dynamic_cast<TH1*> (fModel->PriorContainer(index));
	
	double max_prior = 0;

	if (const_prior) {
		max_prior = 1./fModel->GetParameter(index)->GetRangeWidth();
		const_prior -> SetLineColor(kRed);
	}
	else if (f1_prior) {
		max_prior = f1_prior -> GetMaximum(fModel->GetParameter(index)->GetLowerLimit(),fModel->GetParameter(index)->GetUpperLimit());
		f1_prior -> SetLineColor(kRed);
		f1_prior -> SetLineWidth(1);
	}
	else if (h1_prior) {
		max_prior = h1_prior -> GetMaximum();
		h1_prior -> SetLineColor(kRed);
		h1_prior -> SetStats(false);
		h1_prior -> GetXaxis() -> SetNdivisions(508);
	}
	else {
		if (flag_slice_prior and index<fPriorModel->GetNParameters()) {
			if (fPriorModel->GetNParameters()==2) {
				TH1D * hist = (index==0) ? fPriorModel->GetSlice(0,1)->ProjectionX(Form("projx_%i",BCLog::GetHIndex()))
					: fPriorModel->GetSlice(0,1)->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
				bch1d_prior = new BCH1D(hist);
			} else if (fPriorModel->GetNParameters()==1)
				bch1d_prior = new BCH1D(fPriorModel->GetSlice(index));
		}
		if (!bch1d_prior and fPriorModel->MarginalizedHistogramExists(index))
			bch1d_prior = fPriorModel -> GetMarginalized(index);
		if (!bch1d_prior)
			return 0;
		max_prior = bch1d_prior->GetHistogram()->GetMaximum();
		bch1d_prior -> GetHistogram() -> SetStats(false);
		bch1d_prior -> GetHistogram() -> SetLineColor(kRed);
		bch1d_prior -> GetHistogram() -> GetXaxis() -> SetNdivisions(508);
	}
	
	// if prior doesn't exist, exit
	if (!const_prior and !f1_prior and !h1_prior and !bch1d_prior)
		return 0;

	// Get Posterior
	BCH1D* bch1d_posterior = 0;
	if (flag_slice_post and index<fModel->GetNParameters()) {
		if (fModel->GetNParameters()==2) {
			TH1D * hist = (index==0) ? fModel->GetSlice(0,1)->ProjectionX(Form("projx_%i",BCLog::GetHIndex()))
				: fModel->GetSlice(0,1)->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
			hist -> Scale(1./hist->Integral("width"));
			bch1d_posterior = new BCH1D(hist);
		} else if (fModel->GetNParameters()==1)
			bch1d_posterior = new BCH1D(fModel->GetSlice(index));
	} else if (fModel->MarginalizedHistogramExists(index))
		bch1d_posterior = fModel->GetMarginalized(index);
	
	// if marginal doesn't exist, exit
	if (!bch1d_posterior)
		return 0;
	
	bch1d_posterior -> GetHistogram() -> Scale(1./bch1d_posterior->GetHistogram()->Integral("width"));
	bch1d_posterior -> GetHistogram() -> SetStats(kFALSE);
	
	// get maximum
	double maxy = 1.1 * TMath::Max(max_prior, bch1d_posterior->GetHistogram()->GetMaximum());

	// prepare legend
	TLegend* legend = new TLegend();
	legend->SetBorderSize(0);
	legend->SetFillColor(kWhite);
	legend->SetTextAlign(12);
	legend->SetTextFont(62);
	legend->SetTextSize(0.03);

	// draw axes
	TH2D * h2_axes = new TH2D(TString::Format("h2_axes_%d",getNextIndex()), TString::Format(";%s;P(%s|Data)",fModel->GetObservable(index)->GetLatexName().data(),fModel->GetObservable(index)->GetLatexName().data()),
														10, fModel->GetObservable(index)->GetLowerLimit(), fModel->GetObservable(index)->GetUpperLimit(),
														10, 0, maxy);
	h2_axes -> SetStats(false);
	h2_axes -> GetXaxis() -> SetNdivisions(508);
	h2_axes -> Draw();
	
	// draw prior
	if (const_prior) {
		legend -> AddEntry(const_prior,"prior","L");
		const_prior -> DrawLine(fModel->GetParameter(index)->GetLowerLimit(),max_prior,fModel->GetParameter(index)->GetUpperLimit(),max_prior);
	} else if (f1_prior) {
		legend -> AddEntry(f1_prior,"prior","L");
		f1_prior -> Draw("same");
	} else if (h1_prior) {
		legend -> AddEntry(h1_prior,"prior","L");
		h1_prior -> Draw("same");
	} else if (bch1d_prior) {
		legend -> AddEntry(bch1d_prior->GetHistogram(), "prior", "L");
		bch1d_prior -> Draw(std::string(options_prior+"same"));
	}
	
	// draw posterior
	legend -> AddEntry(bch1d_posterior->GetHistogram(), "posterior", "L");
	bch1d_posterior -> Draw(std::string(options_post+"same"));

	gPad->SetTopMargin(0.02);
   
	// Draw legend on top of histogram
	legend -> SetX1NDC(gPad->GetLeftMargin() + 0.10 * (1.0 - gPad->GetRightMargin() - gPad->GetLeftMargin()));
	legend -> SetX2NDC(1. - gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin() + legend->GetTextSize()*legend->GetNRows();
	legend -> SetY1NDC(1-y1);
	legend -> SetY2NDC(1. - gPad->GetTopMargin());
	legend -> Draw();

	// rescale top margin
	gPad -> SetTopMargin(y1+0.01);

	gPad -> RedrawAxis();

	return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice) {
	
	// Get Prior
	bool const_prior1 = fModel->IsPriorConstant(index1);
	TF1 * f1_prior1   = (const_prior1) ? 0 : dynamic_cast<TF1*> (fModel->PriorContainer(index1));
	TH1 * h1_prior1   = (const_prior1) ? 0 : dynamic_cast<TH1*> (fModel->PriorContainer(index1));
	bool auto_prior1 = const_prior1 or f1_prior1 or h1_prior1;

	bool const_prior2 = fModel->IsPriorConstant(index2);
	TF1 * f1_prior2   = (const_prior2) ? 0 : dynamic_cast<TF1*> (fModel->PriorContainer(index2));
	TH1 * h1_prior2   = (const_prior2) ? 0 : dynamic_cast<TH1*> (fModel->PriorContainer(index2));
	bool auto_prior2 = const_prior2 or f1_prior2 or h1_prior2;


	TH2D * h2d_2dprior = 0;
	TH2D * h2d_2dposterior = 0;

	double prior_mode_1 = fPriorModel->GetBestFitParameter(index1);
	double prior_mode_2 = fPriorModel->GetBestFitParameter(index2);

	if (!auto_prior1 and !auto_prior2) { // neither prior pre-defined
		if (flag_slice and fModel->GetNParameters()==2 and index1<fModel->GetNParameters() and index2<fModel->GetNParameters()) {
			h2d_2dprior = fPriorModel -> GetSlice(index1,index2);
			h2d_2dposterior = fModel  -> GetSlice(index1,index2);
		} else {
			h2d_2dprior = fPriorModel -> GetMarginalizedHistogram(index1,index2);
			h2d_2dposterior = fModel  -> GetMarginalizedHistogram(index1,index2);
		}
	}	else {
		
	}		

		if (!h2d_2dprior or !h2d_2dposterior) // no marginalizations to draw
		return 0;

	// TH1D drawing options
	h2d_2dprior -> SetLineColor(kRed);
	h2d_2dprior -> SetStats(false);
	h2d_2dposterior -> SetStats(false);

	// Create BCH2D's (these normalize the TH1D's)
	BCH2D * bch2d_2dprior     = new BCH2D(h2d_2dprior);
	BCH2D * bch2d_2dposterior = new BCH2D(h2d_2dposterior);

	// Calculate integrated histograms for getting contour line values
	bch2d_2dprior     -> CalculateIntegratedHistogram();
	bch2d_2dposterior -> CalculateIntegratedHistogram();

	// Set contour levels (0.32 = 1 - 68%)
	double level[1] = {bch2d_2dprior -> GetLevel(0.32)};
	h2d_2dprior -> SetContour(1, level);
	h2d_2dprior -> Draw("CONT3");
	level[0] = bch2d_2dposterior -> GetLevel(0.32);
	h2d_2dposterior -> SetContour(1, level);
	h2d_2dposterior -> Draw("CONT3 SAME");

	// create legend
	TLegend * legend2d = new TLegend();
	legend2d->SetBorderSize(0);
	legend2d->SetFillColor(0);
	legend2d->SetTextAlign(12);
	legend2d->SetTextFont(62);
	legend2d->SetTextSize(0.03);
	
  // create markers and arrows
	TMarker * marker_prior = new TMarker();
	marker_prior->SetMarkerStyle(20);
	marker_prior->SetMarkerColor(kRed);
	
	TMarker * marker_posterior = new TMarker();
	marker_posterior->SetMarkerStyle(20);
	marker_posterior->SetMarkerColor(kBlue);
	
	TArrow * arrow = new TArrow();
	arrow->SetArrowSize(0.02);
	arrow->SetLineColor(kBlack);
	//   arrow->SetLineStyle(2);

	double prior_mode_X = fPriorModel->GetBestFitParameter(index1);
	double prior_mode_Y = fPriorModel->GetBestFitParameter(index2);
	std::string marker_prior_text = "";
	if (fModel->IsPriorConstant(index1)) {
		prior_mode_X = fPriorModel -> GetParameter(index1) -> GetRangeCenter();
		marker_prior_text += Form("prior(%s) constant",fPriorModel->GetParameter(index1)->GetLatexName().data());
	}
	if (fModel->IsPriorConstant(index2)) {
		prior_mode_Y = fPriorModel -> GetParameter(index2) -> GetRangeCenter();
		if (!marker_prior_text.empty())
			marker_prior_text += ", ";
		marker_prior_text += Form("prior(%s) constant",fPriorModel->GetParameter(index2)->GetLatexName().data());
	}
	marker_prior_text = (marker_prior_text.empty()) ? "prior mode" : "prior mode* [" + marker_prior_text + "]";

	marker_prior     -> DrawMarker(prior_mode_X,prior_mode_Y);
	marker_posterior -> DrawMarker(fModel->GetBestFitParameter(index1),fModel->GetBestFitParameter(index2));
	arrow            -> DrawArrow(prior_mode_X,prior_mode_Y, fModel->GetBestFitParameter(index1), fModel->GetBestFitParameter(index2));
	
	if (index1==0 and index2==1) {
		legend2d->AddEntry(h2d_2dprior,      "smallest 68% interval(s) of prior", "L");
		legend2d->AddEntry(h2d_2dposterior,  "smallest 68% interval(s) of posterior", "L");
		legend2d->AddEntry(marker_prior,     marker_prior_text.data(), "P");
		legend2d->AddEntry(marker_posterior, "posterior mode", "P");
		legend2d->AddEntry(arrow,            "change in mode", "L");
	}

	gPad->SetTopMargin(0.02);

	// place legend on top of histogram
	legend2d->SetX1NDC(gPad->GetLeftMargin());
	legend2d->SetX2NDC(1.-gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin()-legend2d->GetTextSize()*legend2d->GetNRows();
	legend2d->SetY1NDC(1.-y1);
	legend2d->SetY2NDC(1.-gPad->GetTopMargin());

	legend2d->Draw();

	gPad->SetTopMargin(y1+0.01);

	gPad->RedrawAxis();
	return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintKnowledgeUpdatePlots(const char * filename, unsigned hdiv, unsigned vdiv, std::string options)
{
   // perform analysis
   CalculatePriorModel();

   // option flags
   bool flag_slice = false;

   // check content of options string
   if (options.find("slice") < options.size())
      flag_slice = true;

   std::string file(filename);

   // if file extension is neither .pdf nor .ps, force to .pdf
	 if ( file.rfind(".pdf") != file.size()-4 and file.rfind(".ps") != file.size()-3 )
		 file += ".pdf";

   // create canvas and prepare postscript
   TCanvas * c = new TCanvas(TString::Format("c_%d",getNextIndex()));
   c->cd();
   c->Print(std::string(file + "[").c_str());

	 if (hdiv<1) hdiv = 1;
	 if (vdiv<1) vdiv = 1;
	 int npads = hdiv * vdiv;

	 c -> Divide(hdiv,vdiv);

   // loop over all parameters and draw 1D plots
	 int ndrawn = 0;
	 int nprinted = -1;
	 c -> cd(1);
   for (unsigned i = 0; i < fModel->GetNObservables(); ++i)
		 if(DrawKnowledgeUpdatePlot1D(i, options, options)) {
			 ++ndrawn;
			 if (ndrawn!=0 and ndrawn%npads==0) {
				 c -> Print(file.c_str());
				 nprinted = ndrawn;
				 c -> Clear("D");
			 }
			 c -> cd(ndrawn%npads+1);
		 }
	 if (nprinted<ndrawn)
		 c -> Print(file.c_str());

	 c -> Clear("D");

   // loop over all parameter pairs
	 ndrawn = 0;
	 nprinted = -1;
	 c -> cd(1);
   for (unsigned i = 0; i < fModel->GetNObservables(); ++i)
		 for (unsigned j = i+1; j < fModel->GetNObservables(); ++j)
			 if (DrawKnowledgeUpdatePlot2D(i,j,flag_slice)) {
				 ++ndrawn;
				 if (ndrawn!=0 and ndrawn%npads==0) {
					 c -> Print(file.c_str());
					 nprinted = ndrawn;
					 c -> Clear("D");
				 }
				 c -> cd(ndrawn%npads+1);
			 }
	 if (nprinted<ndrawn)
		 c -> Print(file.c_str());

   // close output
   c->Print(std::string(file + "]").c_str());
   c->Update();

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

   fPriorModel = new BCSummaryPriorModel(fModel);

   // perform marginalization
   fPriorModel->MarginalizeAll();

   // perform minimization
   fPriorModel->FindMode( fPriorModel->GetBestFitParameters() );

   // no error
   return 1;
}
