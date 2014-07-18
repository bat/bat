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
#include "BCModel.h"
#include "BCVariable.h"
#include "BCParameter.h"
#include "BCObservable.h"

#include <TArrow.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>

#include <iostream>

// ---------------------------------------------------------
BCSummaryTool::BCSummaryTool(BCModel * model)
	: BCModel()
	, fModel(0)
{
   // set text style
   gStyle->SetPaintTextFormat(".2g");

	 // set model
	 if (model)
		 SetModel(model);
}

// ---------------------------------------------------------
BCSummaryTool::~BCSummaryTool()
{
}

// ---------------------------------------------------------
void BCSummaryTool::SetModel(BCModel * model) {
	fModel = model;
	SetName((fModel->GetName()+"_summary").data());
}

// ---------------------------------------------------------
int BCSummaryTool::CalculatePriorModel() {
	// Clear Parameters and User-defined Observables
	fParameters.Clear(true);
	fObservables.Clear(true);
	
	if (!fModel)
		return 0;
	
	// copy parameters from model
	for (unsigned i = 0; i < fModel->GetNParameters(); ++i) {
		BCParameter * par = const_cast<BCParameter *>(fModel->GetParameter(i));
		if (fModel->MarginalizedHistogramExists(i))
			// Set binning from marginalization rather than parameter copy
			par -> SetNbins(fModel->GetMarginalizedHistogram(i)->GetNbinsX());
		AddParameter(par);
	}	
	// copy user-observables from model
	for (unsigned i = 0; i < fModel->GetNObservables(); ++i) {
		BCObservable * obs = const_cast<BCObservable *>(fModel->GetObservable(i));
		if (fModel->MarginalizedHistogramExists(i+fModel->GetNParameters()))
			// Set binning from marginalization rather than observable copy
			obs -> SetNbins(fModel->GetMarginalizedHistogram(i+fModel->GetNParameters())->GetNbinsX());
		AddObservable(obs);
	}
	 
	// set default MCMC setup to the one of the original model
	MCMCSetPrecision(fModel);

	// perform marginalization
	MarginalizeAll();

	// perform minimization
	FindMode( GetBestFitParameters() );

	// no error
	return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintKnowledgeUpdatePlot1D(int index, const char * filename, std::string options_post, std::string options_prior)
{
   // perform analysis
	if (!CalculatePriorModel())
		return 0;

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
		max_prior = 1./fModel->GetVariable(index)->GetRangeWidth();
		const_prior -> SetLineColor(kRed);
	}
	else if (f1_prior) {
		max_prior = f1_prior -> GetMaximum(fModel->GetVariable(index)->GetLowerLimit(),fModel->GetVariable(index)->GetUpperLimit());
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
		if (flag_slice_prior and index<GetNParameters()) {
			if (GetNParameters()==2) {
				TH1D * hist = (index==0) ? GetSlice(0,1)->ProjectionX(Form("projx_%i",BCLog::GetHIndex())) : GetSlice(0,1)->ProjectionY(Form("projy_%i",BCLog::GetHIndex()));
				bch1d_prior = new BCH1D(hist);
			} else if (GetNParameters()==1)
				bch1d_prior = new BCH1D(GetSlice(index));
		}
		if (!bch1d_prior and MarginalizedHistogramExists(index))
			bch1d_prior = GetMarginalized(index);
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
	TH2D * h2_axes = new TH2D(TString::Format("h2_axes_%s_%d",GetName().data(),index), TString::Format(";%s;P(%s|Data)",fModel->GetVariable(index)->GetLatexName().data(),fModel->GetVariable(index)->GetLatexName().data()),
														10, fModel->GetVariable(index)->GetLowerLimit(), fModel->GetVariable(index)->GetUpperLimit(),
														10, 0, maxy);
	h2_axes -> SetStats(false);
	h2_axes -> GetXaxis() -> SetNdivisions(508);
	h2_axes -> Draw();
	
	// draw prior
	if (const_prior) {
		legend -> AddEntry(const_prior,"prior","L");
		const_prior -> DrawLine(fModel->GetVariable(index)->GetLowerLimit(),max_prior,fModel->GetVariable(index)->GetUpperLimit(),max_prior);
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
int BCSummaryTool::DrawKnowledgeUpdatePlot2D(unsigned index1, unsigned index2, bool flag_slice, double interval_content) {
	
	if (index1 == index2)
		return 0;
	if (index1 > index2)
		return DrawKnowledgeUpdatePlot2D(index2,index1,flag_slice);

	// Get Posterior
	TH2D * h2d_2dposterior = 0;
	if (flag_slice and fModel->GetNParameters()==2 and index1<fModel->GetNParameters() and index2<fModel->GetNParameters())
		h2d_2dposterior = fModel  -> GetSlice(index1,index2);
	else if (fModel->MarginalizedHistogramExists(index1,index2))
		h2d_2dposterior = fModel  -> GetMarginalizedHistogram(index1,index2);

	if (!h2d_2dposterior)
		return 0;

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

	if (!auto_prior1 or !auto_prior2) { // one or both prior pre-defined
		if (flag_slice and fModel->GetNParameters()==2 and index1<fModel->GetNParameters() and index2<fModel->GetNParameters())
			h2d_2dprior = GetSlice(index1,index2);
		else if (MarginalizedHistogramExists(index1, index2))
			h2d_2dprior = GetMarginalizedHistogram(index1,index2);
	}

	// if not predefined, use the projection of the marginalization
	if (!auto_prior1 and h2d_2dprior)
		h1_prior1 = h2d_2dprior -> ProjectionX(TString::Format("h1_prior1_%s_%d",fModel->GetName().data(),index1));
	if (!auto_prior2 and h2d_2dprior)
		h1_prior2 = h2d_2dprior -> ProjectionY(TString::Format("h1_prior2_%s_%d",fModel->GetName().data(),index2));

	if (!h2d_2dprior)
		h2d_2dprior = GetVariable(index1) -> CreateH2(TString::Format("h2d_2dprior_%s_%d_%d",GetName().data(),index1,index2).Data(),GetVariable(index2));
	
	for (int i = 1; i <= h2d_2dprior->GetNbinsX(); ++i) {
		// x prior
		double x = 1;
		if (f1_prior1)
			x = f1_prior1 -> Eval(h2d_2dprior->GetXaxis()->GetBinCenter(i));
		else if (h1_prior1)
			x = h1_prior1 -> GetBinContent(h1_prior1->FindFixBin(h2d_2dprior->GetXaxis()->GetBinCenter(i)));
		
		for (int j = 1; j <= h2d_2dprior->GetNbinsY(); ++j) {
			// y prior
			double y = 1;
			if (f1_prior2)
				y = f1_prior2 -> Eval(h2d_2dprior->GetYaxis()->GetBinCenter(j));
			else if (h1_prior2)
				y = h1_prior2 -> GetBinContent(h1_prior2->FindFixBin(h2d_2dprior->GetYaxis()->GetBinCenter(j)));
			
			h2d_2dprior -> SetBinContent(i,j,x*y);
		}
	}

	if (!h2d_2dprior)
		return 0;


	// TH2D drawing options
	h2d_2dprior -> SetLineColor(kRed);
	h2d_2dprior -> SetStats(false);
	h2d_2dposterior -> SetStats(false);

	// Create BCH2D's (these normalize the TH2D's)
	BCH2D * bch2d_2dprior     = new BCH2D(h2d_2dprior);
	BCH2D * bch2d_2dposterior = new BCH2D(h2d_2dposterior);

	// Calculate integrated histograms for getting contour line values
	bch2d_2dprior     -> CalculateIntegratedHistogram();
	bch2d_2dposterior -> CalculateIntegratedHistogram();

	// Set contour levels 
	if (interval_content <= 0 or interval_content >= 1)
		interval_content = 68e-2;
	double level[1] = {bch2d_2dprior -> GetLevel(1-interval_content)};
	h2d_2dprior -> SetContour(1, level);
	h2d_2dprior -> Draw("CONT3");
	level[0] = bch2d_2dposterior -> GetLevel(1-interval_content);
	h2d_2dposterior -> SetContour(1, level);
	h2d_2dposterior -> Draw("CONT3 SAME");

	// create legend
	TLegend * legend2d = new TLegend();
	legend2d -> SetBorderSize(0);
	legend2d -> SetFillColor(kWhite);
	legend2d -> SetTextAlign(12);
	legend2d -> SetTextFont(62);
	legend2d -> SetTextSize(0.03);
	
  // create markers and arrows
	TMarker * marker_prior = new TMarker();
	marker_prior -> SetMarkerStyle(20);
	marker_prior -> SetMarkerColor(kRed);
	marker_prior -> SetMarkerSize(1.5*gPad->GetWNDC());
	
	TMarker * marker_posterior = new TMarker();
	marker_posterior -> SetMarkerStyle(20);
	marker_posterior -> SetMarkerColor(kBlue);
	marker_posterior -> SetMarkerSize(1.5*gPad->GetWNDC());
	
	TArrow * arrow = new TArrow();
	arrow -> SetArrowSize(0.02*gPad->GetWNDC());
	arrow -> SetLineColor(kBlack);

	double prior_mode_X = GetBestFitParameter(index1);
	double prior_mode_Y = GetBestFitParameter(index2);

	std::string marker_prior_text = "";
	if (const_prior1) {
		prior_mode_X = GetParameter(index1) -> GetRangeCenter();
		marker_prior_text += Form("prior(%s) constant",GetParameter(index1)->GetLatexName().data());
	}
	if (const_prior2) {
		prior_mode_Y = GetParameter(index2) -> GetRangeCenter();
		if (!marker_prior_text.empty())
			marker_prior_text += ", ";
		marker_prior_text += Form("prior(%s) constant",GetParameter(index2)->GetLatexName().data());
	}
	marker_prior_text = (marker_prior_text.empty()) ? "prior mode" : "prior mode* [" + marker_prior_text + "]";

	marker_prior     -> DrawMarker(prior_mode_X,prior_mode_Y);
	marker_posterior -> DrawMarker(fModel->GetBestFitParameter(index1),fModel->GetBestFitParameter(index2));
	arrow            -> DrawArrow(prior_mode_X,prior_mode_Y, fModel->GetBestFitParameter(index1), fModel->GetBestFitParameter(index2));
	
	legend2d->AddEntry(h2d_2dprior,      TString::Format("smallest %.0f%% interval(s) of prior",     100*interval_content), "L");
	legend2d->AddEntry(h2d_2dposterior,  TString::Format("smallest %.0f%% interval(s) of posterior", 100*interval_content), "L");
	legend2d->AddEntry(marker_prior,     marker_prior_text.data(), "P");
	legend2d->AddEntry(marker_posterior, "posterior mode", "P");
	legend2d->AddEntry(arrow,            "change in mode", "L");
	
	gPad->SetTopMargin(0.02);

	// place legend on top of histogram
	legend2d->SetX1NDC(gPad->GetLeftMargin());
	legend2d->SetX2NDC(1. - gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin() + legend2d->GetTextSize()*legend2d->GetNRows();
	legend2d->SetY1NDC(1.-y1);
	legend2d->SetY2NDC(1. - gPad->GetTopMargin());

	legend2d->Draw();
	
	gPad->SetTopMargin(y1+0.01);
	
	gPad->RedrawAxis();
	return 1;
}

// ---------------------------------------------------------
int BCSummaryTool::PrintKnowledgeUpdatePlots(const char * filename, unsigned hdiv, unsigned vdiv, std::string options, double interval_content)
{
   // perform analysis
	if (!CalculatePriorModel())
		return 0;

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
   TCanvas * c = new TCanvas(TString::Format("c_%s_update",GetName().data()));
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
   for (unsigned i = 0; i < fModel->GetNVariables(); ++i)
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
   for (unsigned i = 0; i < fModel->GetNVariables(); ++i)
		 for (unsigned j = i+1; j < fModel->GetNVariables(); ++j)
			 if (DrawKnowledgeUpdatePlot2D(i,j,flag_slice,interval_content)) {
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
