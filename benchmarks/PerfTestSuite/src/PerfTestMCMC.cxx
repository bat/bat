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
int PerfTestMCMC::PreTest()
{
	//add histograms
	int npar = GetNParameters();

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
		TH2D* hist = new TH2D("", ";Iteration; Auto-correlation]", 1, 0., fMCMCNIterations.at(0), 1, -1., 1.0);
		hist->SetStats(kFALSE);
		hist->Draw();
		fCorrelation.at(i)->Draw("SAMEP");
		AddCanvas(c_corr);
		AddCanvasDescription("Auto-correlation for the parameter.");
	}

	// no error
	return 1;
}

//______________________________________________________________________________
int PerfTestMCMC::RunTest()
{
	// define error code
	int err = 1;

	// call pre-test
	err *= PreTest();

	// perform mcmc
	err *= MarginalizeAll(); 

	// call post-test
	err *= PostTest();

	// return error code
	return err;
}

//______________________________________________________________________________
void PerfTestMCMC::DefineSubtests()
{
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
		MCMCSetNLag(1);
		MCMCSetNIterationsRun(10000);
	}
	else if (precision == PerfTest::kMedium) {
		MCMCSetNLag(5);
		MCMCSetNIterationsRun(100000);
	}
	else if (precision == PerfTest::kDetail) {
		MCMCSetNLag(10);
		MCMCSetNIterationsRun(1000000);
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
	int nchains = MCMCGetNChains(); 

	for (int i = 0; i < npar; ++i) {
		TH2D* hist = fHistCorr.at(i); 

		if (fMCMCNIterations.at(0)>npar)
			for (int j = 0; j < nchains; ++j)
				hist->Fill(fXOld.at(j * npar + i), fMCMCx.at(j * npar + i));

		int counter = (fMCMCNIterations.at(0)-npar)/nlag/npar;
		int countermax = fMCMCNIterationsRun/nlag;

		if (counter % (countermax/500) == 0) {
			(fCorrelation[i])->SetPoint(counter / (countermax/500),
																	fMCMCNIterations.at(0), 
																	hist->GetCorrelationFactor());
		}
	}

	// copy point
	fXOld = fMCMCx; 

}

//______________________________________________________________________________
	
