// ***************************************************************
// This file was created using the ./CreateFitModel.sh script
// ./CreateFitModel.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <TemplateFit.h>

#include <TFile.h>
#include <TH1D.h>

int main()
{
	// ----------------------------------------------------
	// configure BAT
	// ----------------------------------------------------

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// ----------------------------------------------------
	// read histograms from file
	// ----------------------------------------------------

	TFile * file = new TFile("histograms.root", "read");

	TH1D* hist_0 = (TH1D*) file->Get("hist_0");
	TH1D* hist_1 = (TH1D*) file->Get("hist_1");
	TH1D* hist_2 = (TH1D*) file->Get("hist_2");
	TH1D* hist_3 = (TH1D*) file->Get("hist_3");
	TH1D* hist_4 = (TH1D*) file->Get("hist_4");
	TH1D* hist_5 = (TH1D*) file->Get("hist_5");

	// ----------------------------------------------------
	// create new fit model and perform the fit
	// ----------------------------------------------------

	// create new fit model
	TemplateFit * m = new TemplateFit();

	// add histograms
	m->AddHistogram(hist_0, 150.0);
	m->AddHistogram(hist_1, 160.0);
	m->AddHistogram(hist_2, 170.0);
	m->AddHistogram(hist_3, 180.0);
	m->AddHistogram(hist_4, 190.0);
	m->AddHistogram(hist_5, 200.0);

	// set MCMC properties
	//	m->MCMCSetNIterationsRun(1000000);

	// perform fit
	m->MarginalizeAll();
	m->FindMode();

	// print results
	m->PrintAllMarginalized("plots.pdf");
	m->PrintResults("results.txt");
	m->PrintHistograms();
	m->PrintChi2Summary();

	// print goodness-of-fit
	std::cout << " chi2 / dof (chi2-prob) : "
						<< m->CalculateChi2() << "/"
						<< m->GetNDF() << " ("
						<< m->CalculateChi2Prob() << ")" << std::endl;

	// ----------------------------------------------------
	// clean up
	// ----------------------------------------------------

	// delete model
	delete m;

	// close file
  //	file->Close();
  //	delete file;

	// no errors
	return 0;

}

