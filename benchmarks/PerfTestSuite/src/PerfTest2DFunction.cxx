/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TCanvas.h>
#include <TF2.h>
#include <TH2D.h> 
#include <TMath.h>

#include <BAT/BCH2D.h>
#include <BAT/BCAux.h> 

#include <iostream>

#include "include/PerfTest2DFunction.h"

//______________________________________________________________________________
PerfTest2DFunction::PerfTest2DFunction(std::string name, TF2* func) : PerfTest(name)
																																		,	BCModel(name.c_str())
																																		, fFunction(func)
{
	// set test type 
	fTestType = PerfTest::kFunction2D; 

	// set style
	BCAux::SetStyle();

	// set options
	MCMCSetNLag(50);
	MCMCSetNIterationsRun(1000000);

	// manipulate function
	fFunction->SetNDF(100000);

	// get limits
	double xmin = fFunction->GetXmin();
	double xmax = fFunction->GetXmax();

	double ymin = fFunction->GetYmin();
	double ymax = fFunction->GetYmax();

	// add parameters
	AddParameter("x", xmin, xmax);
	AddParameter("y", ymin, ymax);

	DefineSubtests(); 
}
	
//______________________________________________________________________________
PerfTest2DFunction::~PerfTest2DFunction()
{
}
	
//______________________________________________________________________________
int PerfTest2DFunction::Run()
{
	// perform mcmc
	MarginalizeAll(); 
	
	// get histogram
	TH2D* hist_marg = (TH2D*) GetMarginalized( GetParameter(0), GetParameter(1) )->GetHistogram()->Clone(); 
	hist_marg->SetContour(20);
	hist_marg->Scale(hist_marg->GetEntries()/hist_marg->Integral());
	TH2D* hist_diff = (TH2D*) GetMarginalized( GetParameter(0), GetParameter(1) )->GetHistogram()->Clone(); 
	hist_diff->SetContour(20);
	TH2D* hist_func = (TH2D*) GetMarginalized( GetParameter(0), GetParameter(1) )->GetHistogram()->Clone(); 
	hist_func->SetContour(20);

	TH1D* hist_pull = new TH1D("hist_pull", ";#Deltaf/f;N", 11, -5.5, 5.5);

	// calculate ndf 
	int nbinsx = hist_marg->GetNbinsX();
	int nbinsy = hist_marg->GetNbinsY();
	int ndf = nbinsx*nbinsy;

	// calculate norms
	double norm_hist = hist_marg->Integral();

	// get bin widths
	double binwidthx = hist_marg->GetXaxis()->GetBinWidth(1);
	double binwidthy = hist_marg->GetYaxis()->GetBinWidth(1);

	// fill histograms
	for (int i = 1; i <= nbinsx; ++i) {
		for (int j = 1; j <= nbinsy; ++j) {
			double e = fFunction->Integral(hist_marg->GetXaxis()->GetBinCenter(i)-0.5*binwidthx, 
																		 hist_marg->GetXaxis()->GetBinCenter(i)+0.5*binwidthx, 
																		 hist_marg->GetYaxis()->GetBinCenter(j)-0.5*binwidthy,
																		 hist_marg->GetYaxis()->GetBinCenter(j)+0.5*binwidthy
																		 );
			hist_func->SetBinContent(i, j, e);
		}
	}
	hist_func->Scale(norm_hist/hist_func->Integral());

	// calculate chi2
	double chi2 = 0; 
	for (int i = 1; i <= nbinsx; ++i) {
		for (int j = 1; j <= nbinsy; ++j) {
			double n = hist_marg->GetBinContent(i, j);
			double e = hist_func->GetBinContent(i, j); 
			chi2 += (n-e)*(n-e)/e; 
			
			// fill histograms
			hist_func->SetBinContent(i, j, e);
			hist_diff->SetBinContent(i, j, n-e);
			hist_pull->Fill((n-e)/sqrt(e));
		}
	}

	// define test results
	GetSubtest("chi2")->SetTargetValue(ndf);
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kGood,   3.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kFlawed, 5.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kBad,    7.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetTestValue(chi2); 

	// add canvases
	TCanvas* c_func1 = new TCanvas();
	fFunction->Draw();
	AddCanvas(c_func1);

	TCanvas* c_func2 = new TCanvas();
	hist_func->Draw("COLZ");
	AddCanvas(c_func2);

	TCanvas* c_marg = new TCanvas();
	hist_marg->Draw("COLZ");
	AddCanvas(c_marg);

	TCanvas* c_func3 = new TCanvas();
	c_func3->SetLogz(kTRUE);
	hist_func->Draw("COLZ");
	AddCanvas(c_func3);

	TCanvas* c_marg_log = new TCanvas();
	c_marg_log->SetLogz(kTRUE);
	hist_marg->Draw("COLZ");
	AddCanvas(c_marg_log);

	TCanvas* c_diff = new TCanvas();
	c_diff->cd();
	hist_diff->Draw("COLZ");
	AddCanvas(c_diff); 

	TCanvas* c_pull = new TCanvas();
	c_pull->cd();
	hist_pull->Draw();
	AddCanvas(c_pull); 

	// debugKK
	// also add projections?

	// no error 
	return 1; 
}

//______________________________________________________________________________
void PerfTest2DFunction::DefineSubtests()
{
	PerfSubTest * subtest = new PerfSubTest("chi2"); 
	AddSubtest(subtest);
}

//______________________________________________________________________________
	
