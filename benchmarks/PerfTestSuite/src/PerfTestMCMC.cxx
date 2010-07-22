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

#include "include/PerfTestMCMC.h"

//______________________________________________________________________________
PerfTestMCMC::PerfTestMCMC(std::string name) : PerfTest(name)
																						 , BCModel(name.c_str())
																						 , fCorrelation(std::vector<TGraph*>(0))
																						 , fHistCorr(std::vector<TH2D*>(0))
																						 , fXOld(std::vector<double>(0))
{
	// define subtests
	DefineSubtests();
}
	
//______________________________________________________________________________
PerfTestMCMC::~PerfTestMCMC()
{
}
	
//______________________________________________________________________________
int PerfTestMCMC::SetVarPar(double value, std::string name)
{
	if (name == "lag") {
		int n = MCMCGetNIterationsRun();
		int lag = MCMCGetNLag();
		int iter = int( n * value/double(lag) );
		MCMCSetNLag(int(value));
		MCMCSetNIterationsRun(iter);
		return 1;
	}
	else if (name == "iterations") {
		MCMCSetNIterationsRun(int(value));
		return 1;
	}
	else 
		return 0;
}

//______________________________________________________________________________
int PerfTestMCMC::PreTest()
{
	//add histograms
	int npar = GetNParameters();

	// clear histograms
	fHistCorr.clear();
	fCorrelation.clear();

	// loop over parameters
	for (int i = 0; i < npar; ++i) {
		double xmin = GetParameter(i)->GetLowerLimit();
		double xmax = GetParameter(i)->GetLowerLimit();
		TH2D* hist = new TH2D("", "", 100, xmin, xmax, 100, xmin, xmax);
		fHistCorr.push_back(hist);
		fCorrelation.push_back(new TGraph(0));
	}

	return 1;
}

//______________________________________________________________________________
int PerfTestMCMC::PostTest()
{

	// loop over parameters
	int npar = GetNParameters();
	for (int i = 0; i < npar; ++i) {
		TCanvas* c_corr = new TCanvas();
		TH2D* hist = new TH2D("", ";Iteration; Auto-correlation]", 1, 0., MCMCGetNIterationsRun()+1, 1, -1., 1.0);
		hist->SetStats(kFALSE);
		hist->Draw();
		fCorrelation.at(i)->Draw("SAMEP");
		AddCanvas(c_corr);
		AddCanvasDescription("Auto-correlation for the parameter.");
		
		double x; 
		double y;
		fCorrelation.at(i)->GetPoint(fCorrelation.at(i)->GetN()-1, x, y);
		
		// define test results
		GetSubtest(Form("correlation par %i", i))->SetTargetValue(0.0);
		GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kGood,   0.3); 
		GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kFlawed, 0.5); 
		GetSubtest(Form("correlation par %i", i))->SetStatusRegion(PerfSubTest::kBad,    0.7); 
		GetSubtest(Form("correlation par %i", i))->SetTestValue(y); 
		GetSubtest(Form("correlation par %i", i))->SetTestUncertainty(fCorrelation.at(i)->GetRMS(2)); 
		GetSubtest(Form("correlation par %i", i))->SetStatusOff(true);
	}

	// no error
	return 1;
}

//______________________________________________________________________________
int PerfTestMCMC::RunTest()
{
	// define error code
	int err = 1;

	// perform mcmc
	err *= MarginalizeAll(); 

	// return error code
	return err;
}

//______________________________________________________________________________
void PerfTestMCMC::DefineSubtests()
{

	// loop over parameters
	int npar = GetNParameters();
	for (int i = 0; i < npar; ++i) {	
		PerfSubTest * subtest = new PerfSubTest(Form("correlation par %i", i));
		subtest->SetDescription("Calculate the auto-correlation among the points."); 
		AddSubtest(subtest);
	}
}

//______________________________________________________________________________
int PerfTestMCMC::WriteResults()
{
	PerfTest::WriteResults(); 

	PrintResults( Form("%s.log", PerfTest::GetName().c_str()));

	return 1;
}

//______________________________________________________________________________
void PerfTestMCMC::PrecisionSettings(PerfTest::Precision precision)
{
	if (precision == PerfTest::kCoarse) {
	  MCMCSetPrecision(BCEngineMCMC::kLow);
	}
	else if (precision == PerfTest::kMedium) {
	  MCMCSetPrecision(BCEngineMCMC::kMedium);
	}
	else if (precision == PerfTest::kDetail) {
	  MCMCSetPrecision(BCEngineMCMC::kVeryHigh);
	}
	else {
		MCMCSetNLag(5);
		MCMCSetNIterationsRun(100000);
	}
}

//______________________________________________________________________________
void PerfTestMCMC::MCMCUserIterationInterface()
{
	// loop over parameters
	int npar = GetNParameters();
	int nlag = MCMCGetNLag();
	int iteration = MCMCGetCurrentIteration();
	int nchains = MCMCGetNChains(); 

	for (int i = 0; i < npar; ++i) {
		TH2D* hist = fHistCorr.at(i); 

		if (iteration > nlag) {
			for (int j = 0; j < nchains; ++j) 
				hist->Fill(fXOld.at(j * npar + i), fMCMCx.at(j * npar + i));
		}

		if (iteration/nlag % (MCMCGetNIterationsRun()/100/nlag) == 0) {
			(fCorrelation[i])->SetPoint( (iteration/nlag) / (MCMCGetNIterationsRun()/100/nlag), 
																	iteration, 
																	hist->GetCorrelationFactor());
		}
	}

	// copy point
	fXOld = fMCMCx; 

}

//______________________________________________________________________________
	
