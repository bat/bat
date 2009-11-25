#include "SummaryTool.h"
#include <BAT/BCH1D.h> 
#include <BAT/BCH2D.h> 
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <TStyle.h> 

// ---------------------------------------------------------
SummaryTool::SummaryTool()
{
	fModel = 0; 

	// define sum of probabilities for quantiles 
	fSumProb.push_back(0.05); 
	fSumProb.push_back(0.10); 
	fSumProb.push_back(0.1587); 
	fSumProb.push_back(0.50); 
	fSumProb.push_back(0.8413); 
	fSumProb.push_back(0.90); 
	fSumProb.push_back(0.95); 

	gStyle->SetPaintTextFormat("3.0g"); 
};

// ---------------------------------------------------------
SummaryTool::SummaryTool(BCModel* model) : 
	fModel(model)
{
};

// ---------------------------------------------------------
SummaryTool::~SummaryTool()
{}; 

// ---------------------------------------------------------
int SummaryTool::CopySummaryData()
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

	// get number of parameters and quantiles 
	int npar = fModel->GetNParameters(); 
	int nquantiles = int( fSumProb.size() ); 

	// copy information from marginalized distributions
	for (int i = 0; i < npar; ++i) {
		fParName.push_back( (fModel->GetParameter(i)->GetName()).c_str() ); 
		fParMin.push_back( fModel->GetParameter(i)->GetLowerLimit() ); 
		fParMax.push_back( fModel->GetParameter(i)->GetUpperLimit() ); 
		BCH1D* bch1d_temp = fModel->GetMarginalized( fModel->GetParameter(i) ); 
		fMean.push_back( bch1d_temp->GetMean() ); 
		fMargMode.push_back( bch1d_temp->GetMode() ); 
		fGlobalMode.push_back ( (fModel->GetBestFitParameters()).at(i) ); 
		for (int j = 0; j < nquantiles; ++j) 
			fQuantiles.push_back( bch1d_temp->GetQuantile( fSumProb.at(j) ) ); 
		std::vector <double> intervals = bch1d_temp->GetSmallestIntervals(); 
		int nintervals = int(intervals.size() / 5); 
		fSmallInt.push_back(nintervals); 
		fSmallInt.insert( fSmallInt.end(), intervals.begin(), intervals.end() ); 
		for (int j = 0; j < npar; ++j) { 
			if (i!=j) {
				BCH2D* bch2d_temp = fModel->GetMarginalized(fModel->GetParameter(i), 
																										fModel->GetParameter(j)); 
				fCorrCoeff.push_back( bch2d_temp->GetHistogram()->GetCorrelationFactor() ); 
			}
			else
				fCorrCoeff.push_back(1.0); 
		}
		fRMS.push_back( bch1d_temp->GetRMS() ); 
	}

	// no error
	return 1; 
}; 

// ---------------------------------------------------------
int SummaryTool::PrintParameterPlot(const char* filename)
{
	// copy summary data 
	if (!CopySummaryData())
		return 0; 

	// get number of parameters and quantiles 
	int npar = fModel->GetNParameters(); 
	int nquantiles = int( fSumProb.size() ); 
	
	// create histogram
	TH1D* hist_axes = new TH1D("hist_axes_par", ";;Scaled parameter range [a.u.]", 
														 npar, -0.5, npar-0.5); 
	hist_axes->SetStats(kFALSE); 
	for (int i = 0; i < npar; ++i) 
		hist_axes->GetXaxis()->SetBinLabel( i+1, fParName.at(i) ); 
	hist_axes->GetXaxis()->SetTickLength(0.0); 
	hist_axes->GetYaxis()->SetRangeUser(0.0, 1.0); 
	hist_axes->GetYaxis()->SetTickLength(0.0); 

	// create graphs
	TGraphErrors* graph_quantiles = new TGraphErrors(npar*nquantiles); 
	graph_quantiles->SetMarkerSize(0); 
	graph_quantiles->SetLineColor(38); 
	graph_quantiles->SetLineStyle(2); 

	TGraphErrors* graph_mean = new TGraphErrors(npar); 
	graph_mean->SetMarkerColor(kBlack);
	graph_mean->SetMarkerStyle(20); 

	TGraphErrors* graph_mode = new TGraphErrors(npar); 
	graph_mode->SetMarkerColor(kRed);
	graph_mode->SetMarkerStyle(20); 

	TGraphAsymmErrors* graph_intervals = new TGraphAsymmErrors(0); 
	graph_intervals->SetFillColor(kYellow); 
	graph_intervals->SetLineStyle(2); 
	graph_intervals->SetLineColor(kRed); 
	graph_intervals->SetMarkerSize(0); 

	// fill graphs
	int indexintervals = 0; 

	for (int i = 0; i < npar; ++i) {

		// fill graph quantiles 
		for (int j = 0; j < nquantiles; ++j) {
			graph_quantiles->SetPoint(i*nquantiles+j,
																double(i), 
																(fQuantiles.at(i*nquantiles+j) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i))); 
			graph_quantiles->SetPointError(i*nquantiles+j, 0.5, 0.0); 
		}

		// fill graph mean 
		graph_mean->SetPoint(i, double(i), (fMean.at(i) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
		graph_mean->SetPointError(i, 0.0, fRMS.at(i)/(fParMax.at(i)-fParMin.at(i)));

		// fill graph mode 
		graph_mode->SetPoint(i, double(i), (fGlobalMode.at(i) - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));

		// fill graph smallest intervals 
		int nintervals = fSmallInt.at(indexintervals++); 
		for (int j = 0; j < nintervals; ++j) {
			double xmin = fSmallInt.at(indexintervals++); 
			double xmax = fSmallInt.at(indexintervals++); 
			indexintervals++; 
			double xlocalmaxpos = fSmallInt.at(indexintervals++); 
			indexintervals++; 
			int npoints = graph_intervals->GetN(); 
			graph_intervals->SetPoint(npoints,
																double(i), 
																(xlocalmaxpos - fParMin.at(i))/(fParMax.at(i)-fParMin.at(i)));
			graph_intervals->SetPointError(npoints, 
																		 0.5, 0.5, 
																		 (xlocalmaxpos - xmin)/(fParMax.at(i)-fParMin.at(i)), 
																		 (xmax - xlocalmaxpos)/(fParMax.at(i)-fParMin.at(i))); 
			}
	}

	// print to file 
	TCanvas * c_par = new TCanvas("c_corr"); 
	c_par->cd(); 
	hist_axes->Draw(); 
	graph_intervals->DrawClone("SAME2"); 
	for (int i = 0; i < graph_intervals->GetN(); ++i)
		graph_intervals->SetPointError(i, 0.5, 0.5, 0.0, 0.0); 
	graph_intervals->Draw("SAMEPZ"); 
	graph_quantiles->Draw("SAMEPZ"); 
	graph_mean->Draw("SAMEP"); 
	graph_mode->Draw("SAMEP"); 
	c_par->Print(filename); 

	// no error 
	return 1; 
}

// ---------------------------------------------------------
int SummaryTool::PrintCorrelationPlot(const char* filename)
{
	// copy summary data 
	if (!CopySummaryData())
		return 0; 

	// get number of parameters 
	int npar = fModel->GetNParameters(); 

	// create histogram
	TH2D * hist_corr = new TH2D("hist_corr", ";;", 
															npar, -0.5, npar-0.5,
															npar, -0.5, npar-0.5); 
	hist_corr->SetStats(kFALSE); 
	hist_corr->GetXaxis()->SetTickLength(0.0); 
	hist_corr->GetYaxis()->SetTickLength(0.0); 
	hist_corr->GetZaxis()->SetRangeUser(-1.0, 1.0); 

	for (int i = 0; i < npar; ++i) {
		hist_corr->GetXaxis()->SetBinLabel( i+1, fParName.at(i) ); 
		hist_corr->GetYaxis()->SetBinLabel( i+1, fParName.at(i) ); 
	}

	// fill plot
	for (int i = 0; i < npar; ++i) {
		for (int j = 0; j < npar; ++j) {
			// debugKK
			// check if this is correct or if it has to be reversed 
			int index = i * npar + j; 
			double corr = fCorrCoeff.at(index); 
			hist_corr->SetBinContent(i+1, j+1, corr); 
		}
	}

	// print to file 
	TCanvas * c_corr = new TCanvas("c_corr"); 
	c_corr->cd(); 
	hist_corr->Draw("COLZTEXT"); 
	c_corr->Print(filename); 

	// no error 
	return 1; 
}

// // ---------------------------------------------------------
// int SummaryTool::Print2DOverviewPlots(const char* filename)
// {
// 	// copy summary data 
// 	if (!CopySummaryData())
// 		return 0; 

// 	// get number of parameters 
// 	int npar = fModel->GetNParameters(); 

// 	TCanvas * c_2doverview = new TCanvas("c_2doverview"); 
// 	c_2doverview->Divide(npar+1, npar+1, 0.005, 0.005); 

// 		for (int i = 0; i < npar; ++i) { 
// 		for (int j = 1; j < npar+1; ++j) {
// 			c_2doverview->cd(1+i*(npar+1)+j); 
// 			gPad->SetBottomMargin(0); 
// 			gPad->SetTopMargin(0); 
// 			gPad->SetLeftMargin(0); 
// 			gPad->SetRightMargin(0); 
// 			->DrawClone(); 
// 		}
// 	}
// 	for (int i = 1; i < npar+1; ++i) { 
// 		int index = (npar+1) * npar + i + 1; 
// 		c_2doverview->cd(index); 
// 		TPaveText * pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "NDC"); 
// 		pt->SetTextAlign(22); 
// 		pt->SetTextSize(0.1);
// 		pt->SetBorderSize(0); 
// 		pt->SetFillStyle(0);
// 		pt->AddText(fParName.at(i-1)); 
// 		pt->Draw(); 
// 	}


// 	// no error 
// 	return 1;
// }

// ---------------------------------------------------------
