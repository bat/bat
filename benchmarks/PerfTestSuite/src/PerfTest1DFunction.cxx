/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h> 
#include <TMath.h>

#include <BAT/BCH1D.h>
#include <BAT/BCAux.h> 

#include <iostream>

#include "include/PerfTest1DFunction.h"

//______________________________________________________________________________
PerfTest1DFunction::PerfTest1DFunction(std::string name, TF1* func) : PerfTest(name)
																																		,	BCModel(name.c_str())
																																		, fFunction(func)
{
	// set test type 
	fTestType = PerfTest::kFunction1D; 

	// set style
	BCAux::SetStyle();

	// set options
	MCMCSetNLag(10);
	MCMCSetNIterationsRun(1000000);

	// manipulate function
	fFunction->SetNDF(100000);

	// get limits
	double xmin = fFunction->GetXmin();
	double xmax = fFunction->GetXmax();

	// add parameter
	AddParameter("x", xmin, xmax);

	DefineSubtests(); 
}
	
//______________________________________________________________________________
PerfTest1DFunction::~PerfTest1DFunction()
{
}
	
//______________________________________________________________________________
int PerfTest1DFunction::Run()
{
	// perform mcmc
	MarginalizeAll(); 
	
	// get histogram
	TH1D* hist_marg = (TH1D*) GetMarginalized( GetParameter(0) )->GetHistogram()->Clone(); 
	TH1D* hist_diff = (TH1D*) GetMarginalized( GetParameter(0) )->GetHistogram()->Clone(); 
	TH1D* hist_line = new TH1D(*hist_marg);
	TH1D* hist_diff_1sigma = new TH1D(*(GetMarginalized( GetParameter(0) )->GetHistogram()));
	hist_diff_1sigma->SetFillColor(kGreen);
	hist_diff_1sigma->SetFillStyle(1001);
	hist_diff_1sigma->SetMarkerSize(0); 
	TH1D* hist_diff_2sigma = new TH1D(*hist_diff_1sigma);
	hist_diff_2sigma->SetFillColor(kYellow);
	TH1D* hist_diff_3sigma = new TH1D(*hist_diff_1sigma);
	hist_diff_3sigma->SetFillColor(kRed);

	// calculate ndf 
	int ndf = hist_marg->GetNbinsX();

	// calculate norms
	double norm_hist = hist_marg->Integral();
	double norm_func = fFunction->Integral(fFunction->GetXmin(), fFunction->GetXmax());

	// get bin width
	double binwidth = hist_marg->GetBinWidth(1);

	// calculate chi2
	double chi2 = 0; 
	for (int i = 1; i <= ndf; ++i) {
		double n = hist_marg->GetBinContent(i);
		//		double e = norm_hist * fFunction->Eval(hist_marg->GetBinCenter(i)) * binwidth / norm_func; 
		double e = norm_hist * fFunction->Integral(hist_marg->GetBinCenter(i)-0.5*binwidth, hist_marg->GetBinCenter(i)+0.5*binwidth) / norm_func; 
		chi2 += (n-e)*(n-e)/e; 

		// fill histograms
		hist_diff->SetBinContent(i, n-e);
		hist_line->SetBinContent(i, 0);
		hist_diff_1sigma->SetBinError(i, sqrt(e));		
		hist_diff_1sigma->SetBinContent(i, 0);		
		hist_diff_2sigma->SetBinError(i, 2.*sqrt(e));		
		hist_diff_2sigma->SetBinContent(i, 0);		
		hist_diff_3sigma->SetBinError(i, 3.*sqrt(e));		
		hist_diff_3sigma->SetBinContent(i, 0);		
	}
	
	// calculate quantiles
	double quantiles_hist[9]; 
	double quantiles_func[9]; 
	double probsum[9] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; 

	fFunction->GetQuantiles(9, quantiles_func, probsum);
	hist_marg->GetQuantiles(9, quantiles_hist, probsum);
	
	// define test results
	GetSubtest("chi2")->SetTargetValue(ndf);
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kGood, 3.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kFlawed, 5.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kBad, 7.0*sqrt(2.0*ndf)); 
	GetSubtest("chi2")->SetTestValue(chi2); 

	GetSubtest("mean")->SetTargetValue(fFunction->Mean( fFunction->GetXmin(), fFunction->GetXmax()) );
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kGood, 3.*hist_marg->GetRMS()/sqrt(ndf));
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kFlawed, 5.*hist_marg->GetRMS()/sqrt(ndf));
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kBad, 7.*hist_marg->GetRMS()/sqrt(ndf));
	GetSubtest("mean")->SetTestValue(hist_marg->GetMean()); 

	GetSubtest("mode")->SetTargetValue(fFunction->GetMaximumX());
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kGood, binwidth);
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kFlawed, 2.*binwidth);
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kBad, 3.*binwidth);
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kFatal, binwidth*ndf);
	GetSubtest("mode")->SetTestValue(hist_marg->GetBinCenter(hist_marg->GetMaximumBin())); 

	for (int i = 0; i < 9; ++i) {
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetTargetValue(quantiles_func[i]);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kGood, binwidth);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kFlawed, 2.*binwidth);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kBad, 3.*binwidth);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kFatal, binwidth*ndf);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetTestValue(quantiles_hist[i]); 
	} 

	// add canvases
	TCanvas* c_marg = new TCanvas();
	hist_marg->Scale(1./norm_hist/hist_marg->GetBinWidth(1) * norm_func);
	hist_marg->Draw();
	fFunction->Draw("SAMEC");
	AddCanvas(c_marg);

	TCanvas* c_marg_log = new TCanvas();
	c_marg_log->SetLogy(kTRUE);
	hist_marg->Draw();
	fFunction->Draw("SAMEC");
	AddCanvas(c_marg_log);

	TCanvas* c2 = new TCanvas();
	hist_diff_3sigma->Draw("E3");
	hist_diff_2sigma->Draw("SAMEE3");
	hist_diff_1sigma->Draw("SAMECE3");
	hist_diff->Draw("SAMEP");
	hist_line->Draw("SAMEC");
	hist_diff_3sigma->GetYaxis()->SetRangeUser(1.1*TMath::Min(hist_diff->GetMinimum(), 
																														- 3.0*sqrt(norm_hist*fFunction->GetMaximum()*binwidth/norm_func)),
																						 1.1*TMath::Max(hist_diff->GetMaximum(),
																														3.0*sqrt(norm_hist*fFunction->GetMaximum()*binwidth/norm_func)));
	AddCanvas(c2);

	// no error 
	return 1; 
}

//______________________________________________________________________________
void PerfTest1DFunction::DefineSubtests()
{
	PerfSubTest * subtest = new PerfSubTest("chi2"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("mean"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("mode"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile10"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile20"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile30"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile40"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile50"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile60"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile70"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile80"); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile90"); 
	AddSubtest(subtest);
}

//______________________________________________________________________________
	
