/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h> 
#include <TH2D.h> 
#include <TMath.h>
#include <TROOT.h> 
#include <TStyle.h> 
#include <TRandom3.h> 
#include <TGraph.h>

#include <BAT/BCH1D.h>
#include <BAT/BCAux.h> 
#include <BAT/BCMath.h> 

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
	//	gStyle->SetOptFit(1);
	//	gStyle->SetOptStat(111111);

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
	
	// define histograms
	TH1D* hist_marg = (TH1D*) GetMarginalized( GetParameter(0) )->GetHistogram()->Clone(); 
	TH1D* hist_func = (TH1D*) GetMarginalized( GetParameter(0) )->GetHistogram()->Clone(); 
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
	TH1D* hist_pull = new TH1D("", ";#Deltaf/sqrt(f);N", 50, -5.0, 5.0);

	// calculate nbins 
	int nbins = hist_marg->GetNbinsX();

	// calculate norms
	double norm_hist = hist_marg->Integral();
	double norm_func = fFunction->Integral(fFunction->GetXmin(), fFunction->GetXmax());

	// get bin width
	double binwidth = hist_marg->GetBinWidth(1);

	// calculate chi2
	double chi2 = 0; 
	for (int i = 1; i <= nbins; ++i) {
		double n = hist_marg->GetBinContent(i);
		//		double e = norm_hist * fFunction->Eval(hist_marg->GetBinCenter(i)) * binwidth / norm_func; 
		double e = norm_hist * fFunction->Integral(hist_marg->GetBinCenter(i)-0.5*binwidth, hist_marg->GetBinCenter(i)+0.5*binwidth) / norm_func; 
		chi2 += (n-e)*(n-e)/e; 

		// fill histograms
		hist_diff->SetBinContent(i, n-e);
		hist_func->SetBinContent(i, e);
		hist_pull->Fill( (n-e)/sqrt(e) );
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

	// calculate KS probability
	double KS = hist_marg->KolmogorovTest(hist_func);

	// calculate variance of variance
	double mu2 = fFunction->CentralMoment( 2, fFunction->GetXmin(), fFunction->GetXmax() );
	double mu4 = fFunction->CentralMoment( 4, fFunction->GetXmin(), fFunction->GetXmax() );
	double var_s2 = 1./double(nbins) * (mu4 - (double(nbins)-3)/(double(nbins)-1)*mu2*mu2);

	// Monte Carlo estimate of the uncertainties
	TRandom3* rnd = new TRandom3(0);

	double mode_mean = 0;
	double mode_variance = 0;
	double mean_mean = 0;
	double mean_variance = 0;
	double quantile_mean[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	double quantile_variance[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	int nhist = 10000;

	// loop over histograms
	for (int ihist = 1; ihist <= nhist; ++ihist) {
		TH1D* hist = (TH1D*) hist_func->Clone(); 
		// loop over bins
		for (int ibin = 1; ibin <= nbins; ++ibin) {
			double n = hist->GetBinContent(ibin);
			hist->SetBinContent(ibin, rnd->Poisson(n));
		}
		// calculate mode, etc.
		int maxbin = hist->GetMaximumBin();
		double mode = hist->GetBinCenter(maxbin);
		double mean = hist->GetMean();
		double q[9]; 
		double p[9] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; 
		hist->GetQuantiles(9, q, p);
		
		// update statistics: mode
		mode_mean += (mode-mode_mean)/double(ihist);
		if (ihist>1)
			mode_variance = (1.-1./double(ihist))*mode_variance
				+ (mode-mode_mean) * (mode-mode_mean)/double(ihist-1);											
		// update statistics: mean
		mean_mean += (mean-mean_mean)/double(ihist);
		if (ihist>1)
			mean_variance = (1.-1./double(ihist))*mean_variance
				+ (mean-mean_mean) * (mean-mean_mean)/double(ihist-1);											
		// update statistics: quantiles
		for (int iquant = 0; iquant<9; ++iquant) {
			quantile_mean[iquant] += (q[iquant]-quantile_mean[iquant])/double(ihist);
			if (ihist>1)
				quantile_variance[iquant] = (1.-1./double(ihist))*quantile_variance[iquant]
					+ (q[iquant]-quantile_mean[iquant]) * (q[iquant]-quantile_mean[iquant])/double(ihist-1);											
		}
		// free memory
		delete hist;
	}	

	// define test results
	GetSubtest("chi2")->SetTargetValue(nbins);
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kGood, 3.0*sqrt(2.0*nbins)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kFlawed, 5.0*sqrt(2.0*nbins)); 
	GetSubtest("chi2")->SetStatusRegion(PerfSubTest::kBad, 7.0*sqrt(2.0*nbins)); 
	GetSubtest("chi2")->SetTestValue(chi2); 

	GetSubtest("KS")->SetTargetValue(1);
	GetSubtest("KS")->SetStatusRegion(PerfSubTest::kGood, 1.-0.05); 
	GetSubtest("KS")->SetStatusRegion(PerfSubTest::kFlawed, 1.-0.01); 
	GetSubtest("KS")->SetStatusRegion(PerfSubTest::kBad, 1.-0.0001); 
	GetSubtest("KS")->SetTestValue(KS); 

	GetSubtest("mean")->SetTargetValue(fFunction->Mean( fFunction->GetXmin(), fFunction->GetXmax()) );
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kGood, 3.*sqrt(mean_variance));
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kFlawed, 5.*sqrt(mean_variance));
	GetSubtest("mean")->SetStatusRegion(PerfSubTest::kBad, 7.*sqrt(mean_variance));
	GetSubtest("mean")->SetTestValue(hist_marg->GetMean()); 

	GetSubtest("mode")->SetTargetValue(fFunction->GetMaximumX());
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kGood, 3.*sqrt(mode_variance));
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kFlawed, 5.*sqrt(mode_variance));
	GetSubtest("mode")->SetStatusRegion(PerfSubTest::kBad, 7.*sqrt(mode_variance));
	GetSubtest("mode")->SetTestValue(hist_marg->GetBinCenter(hist_marg->GetMaximumBin())); 

	GetSubtest("variance")->SetTargetValue(fFunction->Variance( fFunction->GetXmin(), fFunction->GetXmax()) );
	GetSubtest("variance")->SetStatusRegion(PerfSubTest::kGood, 3.*sqrt(var_s2));
	GetSubtest("variance")->SetStatusRegion(PerfSubTest::kFlawed, 5.*sqrt(var_s2));
	GetSubtest("variance")->SetStatusRegion(PerfSubTest::kBad, 7.*sqrt(var_s2));
	GetSubtest("variance")->SetTestValue(hist_marg->GetRMS()*hist_marg->GetRMS()*double(nbins*nbins)/(double(nbins)-1)/(double(nbins)-1)); 

	for (int i = 0; i < 9; ++i) {
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetTargetValue(quantiles_func[i]);
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kGood, 3.*sqrt(quantile_variance[i]));
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kFlawed, 5.*sqrt(quantile_variance[i]));
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetStatusRegion(PerfSubTest::kBad, 7.*sqrt(quantile_variance[i]));
		GetSubtest(Form("quantile%i", int(probsum[i]*100)))->SetTestValue(quantiles_hist[i]); 
	} 

	// define graph and axes histogram
	int nsubtests = GetNSubtests();
	TGraph* graph_summary = new TGraph(nsubtests);
	graph_summary->SetMarkerSize(2);
	graph_summary->SetMarkerColor(kBlack);
	graph_summary->SetMarkerStyle(20);
	TH2D* hist_axes = new TH2D("", ";#Deltat/#sigma_{t};Test index", 1, -5., 5., 1, -0.99, nsubtests-0.000001);
	hist_axes->SetStats(kFALSE);

	// loop over tests
	for (int i = 0; i < nsubtests; ++i){
		graph_summary->SetPoint(i, (GetSubtest(i)->GetTestValue() - GetSubtest(i)->GetTargetValue())/GetSubtest(i)->GetStatusRegion(PerfSubTest::kGood)*3., i);
	}

	// add canvases
	TCanvas* c_marg = new TCanvas();
	hist_marg->Scale(1./norm_hist/hist_marg->GetBinWidth(1) * norm_func);
	hist_marg->Draw();
	hist_marg->SetYTitle("f(x)");
	fFunction->Draw("SAMEC");
	AddCanvas(c_marg);
	AddCanvasDescription("Distribution from MCMC and analytic function.");

	TCanvas* c_marg_log = new TCanvas();
	c_marg_log->SetLogy(kTRUE);
	hist_marg->Draw();
	fFunction->Draw("SAMEC");
	AddCanvas(c_marg_log);
	AddCanvasDescription("Distribution from MCMC and analytic function in log-scale.");

	TCanvas* c2 = new TCanvas();
	hist_diff_3sigma->SetYTitle("f(x_{i})-f_{MCMC}(x_{i})");
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
	AddCanvasDescription("Difference between the distribution from MCMC and the analytic function. The one, two and three sigma uncertainty bands are colored green, yellow and red, respectively.");

	TCanvas* c_pull = new TCanvas();
	c_pull->cd();
	hist_pull->Draw();
	TF1* g = new TF1("g", "[0]/sqrt(2.0*TMath::Pi())*exp(-x*x/2.)", -5.0, 5.0);
	g->SetParameter(0, hist_pull->Integral("")*hist_pull->GetBinWidth(1));
	g->Draw("SAMEP");
	AddCanvas(c_pull); 
	AddCanvasDescription("Pull between the distribution from MCMC and the analytic function. The Gaussian has a mean value of 0 and a standard deviation of 1 (not fitted).");

	TCanvas* c_summary = new TCanvas();
	c_summary->cd();
	c_summary->SetGridx();
	c_summary->SetGridy();
	hist_axes->Draw();
	graph_summary->Draw("SAMEP");
	AddCanvas(c_summary);
	AddCanvasDescription("Summary of subtest values.");

	// free memory
	delete rnd;

	// no error 
	return 1; 
}

//______________________________________________________________________________
void PerfTest1DFunction::DefineSubtests()
{
	PerfSubTest * subtest = new PerfSubTest("chi2"); 
	subtest->SetDescription("Calculate &chi;<sup>2</sup> and compare with prediction for dof=number of bins. <br> Tolerance good: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 3 &middot; (2 dof)<sup>1/2</sup>, <br> Tolerance flawed: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 5 &middot; (2 dof)<sup>1/2</sup>, <br> Tolerance bad: |&chi;<sup>2</sup>-E[&chi;<sup>2</sup>]| < 7 &middot; (2 dof)<sup>1/2</sup>."); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("KS"); 
	subtest->SetDescription("Calculate the Kolmogorov-Smirnov probability based on the ROOT implemention. <br> Tolerance good: KS prob > 0.05, <br> Tolerance flawed:  KS prob > 0.01<br> Tolerance bad: KS prob > 0.0001."); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("mean"); 
	subtest->SetDescription("Compare sample mean, &lt;x&gt;, with expectation value of function, E[x].<br> Tolerance good: |&lt;x&gt; -E[x]| < 3 &middot; (V[x]/n)<sup>1/2</sup>,</br>Tolerance flawed: |&lt;x&gt; -E[x]| < 5 &middot; (V[x]/n)<sup>1/2</sup>,</br>Tolerance bad: |&lt;x&gt; -E[x]| < 7 &middot; (V[x]/n)<sup>1/2</sup>."); 
	AddSubtest(subtest);

	subtest = new PerfSubTest("mode"); 
	subtest->SetDescription("Compare mode of distribution with mode of the analytic function. </br> Tolerance good: |x<sup>*</sup>-mode| < 3 &middot V[mode]<sup>1/2</sup>, </br> Tolerance flawed: |x<sup>*</sup>-mode| < 5 &middot V[mode]<sup>1/2</sup> bin widths, <br> Tolerance bad: |x<sup>*</sup>-mode| < 7 &middot V[mode]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("variance"); 
	subtest->SetDescription("Compare sample variance s<sup>2</sup> of distribution with variance of function. </br> Tolerance good: 3 &middot; V[s<sup>2</sup>]<sup>1/2</sup>, </br> Tolerance flawed: 5 &middot; V[s<sup>2</sup>]<sup>1/2</sup>, <br> Tolerance bad: 7 &middot; V[s<sup>2</sup>]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile10"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile20"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile30"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile40"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile50"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile60"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile70"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile80"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);

	subtest = new PerfSubTest("quantile90"); 
	subtest->SetDescription("Compare quantile of distribution from MCMC with the quantile of analytic function. </br> Tolerance good: |q_{X}-E[q_{X}]|<3&middot;V[q]<sup>1/2</sup>, </br> Tolerance flawed: |q_{X}-E[q_{X}]|<5&middot;V[q]<sup>1/2</sup>, <br> Tolerance bad: |q_{X}-E[q_{X}]|<7&middot;V[q]<sup>1/2</sup>.");
	AddSubtest(subtest);


}

//______________________________________________________________________________
int PerfTest1DFunction::WriteResults()
{
	PerfTest::WriteResults(); 

	PrintResults( Form("%s.log", PerfTest::GetName().c_str()));

	return 1;
}

//______________________________________________________________________________
	
