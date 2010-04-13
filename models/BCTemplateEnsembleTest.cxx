#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom3.h>

#include <math.h>
#include <vector>
#include <iostream>

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include "BCTemplateFitter.h"
#include "BCTemplateEnsembleTest.h"

using namespace std;

// ---------------------------------------------------------
BCTemplateEnsembleTest::BCTemplateEnsembleTest()
	: fTemplateFitter(0)
	, fFile(0)
	, fTree(0)
	, fEnsembleCounter(0)
	, fEnsembleExpectation(0)
	, fNEnsembles(0)
	, fFlagMCMC(true)
{
	fRandom = new TRandom3(0);
}

// ---------------------------------------------------------
BCTemplateEnsembleTest::~BCTemplateEnsembleTest()
{
	if (fRandom)
		delete fRandom;
}

// ---------------------------------------------------------
int BCTemplateEnsembleTest::SetEnsembleTemplate(TH1D hist)
{
	// calculate integral
	double integral = hist.Integral();

	// check if integral is ok
	if (integral <= 0) {
		std::cout << "Template not valid. Integral is lower or equal to 0." << std::endl;
		return 0;
	}

	// set template
	fEnsembleTemplate = hist;

	// scale template
	fEnsembleTemplate.Scale(1.0/integral);

	// no error
	return 1;
};

// ---------------------------------------------------------
int BCTemplateEnsembleTest::PerformEnsembleTest()
{
	// set log level to nothing
	BCLog::LogLevel ll = BCLog::GetLogLevelScreen();
	BCLog::SetLogLevel(BCLog::nothing);

	// Prepare tree
	PrepareTree();

	// loop over ensembles
	for(int j = 0; j < fNEnsembles; j++){

		// print status
		if ((j+1) % 100 == 0 && j > 0)
			cout << "Fraction of ensembles analyzed: " << double(j+1) / double(fNEnsembles) * 100 << "%" << std::endl;

		// create new ensemble
		TH1D * ensemble = BuildEnsemble();

		// set ensemble as new data set
		fTemplateFitter->SetData(*ensemble);

		// find mode
		fTemplateFitter->FindMode();

		// perform MCMC
		if(fFlagMCMC) {
			fTemplateFitter->MarginalizeAll();

			// get number of parameters
			int npar = fTemplateFitter->GetNParameters();

			// loop over parameters and set tree variables
			for (int i = 0; i < npar; ++i) {
				BCH1D * hist = fTemplateFitter->GetMarginalized(fTemplateFitter->GetParameter(i));
				fOutParModeMarg[i]       = hist->GetMode();
				fOutParMedianMarg[i]     = hist->GetMedian();
				fOutParMeanMarg[i]       = hist->GetMean();
				fOutParRMSMarg[i]        = hist->GetRMS();
				fOutParErrorUpMarg[i]    = hist->GetQuantile(0.84)-hist->GetMode();
				fOutParErrorDownMarg[i]  = hist->GetMode()-hist->GetQuantile(0.16);
				fOutParQuantile5Marg[i]  = hist->GetQuantile(0.05);
				fOutParQuantile10Marg[i] = hist->GetQuantile(0.10);
				fOutParQuantile90Marg[i] = hist->GetQuantile(0.90);
				fOutParQuantile95Marg[i] = hist->GetQuantile(0.95);
			}
		}

		if (fFlagMCMC) {
			int nratios = fTemplateFitter->GetNRatios();
			for (int i = 0; i < nratios; ++i) {
				TH1D histtemp = fTemplateFitter->GetHistRatio1D(i);
 				BCH1D * hist = new BCH1D( &histtemp );
 				fOutRatioModeMarg[i]       = hist->GetMode();
 				fOutRatioMedianMarg[i]     = hist->GetMedian();
 				fOutRatioMeanMarg[i]       = hist->GetMean();
 				fOutRatioRMSMarg[i]        = hist->GetRMS();
 				fOutRatioErrorUpMarg[i]    = hist->GetQuantile(0.84)-hist->GetMode();
				fOutRatioErrorDownMarg[i]  = hist->GetMode()-hist->GetQuantile(0.16);
 				fOutRatioQuantile5Marg[i]  = hist->GetQuantile(0.05);
 				fOutRatioQuantile10Marg[i] = hist->GetQuantile(0.10);
 				fOutRatioQuantile90Marg[i] = hist->GetQuantile(0.90);
 				fOutRatioQuantile95Marg[i] = hist->GetQuantile(0.95);
			}
		}

		// set tree variables
		fOutParModeGlobal      = fTemplateFitter->GetBestFitParameters();
		fOutParErrorUpGlobal   = fTemplateFitter->GetBestFitParameterErrors();
		fOutParErrorDownGlobal = fTemplateFitter->GetBestFitParameterErrors();
		fOutChi2               = fTemplateFitter->CalculateChi2();
		fOutNDF                = fTemplateFitter->GetNDF();
		fOutChi2Prob           = fTemplateFitter->CalculateChi2Prob();
		fOutKSProb             = fTemplateFitter->CalculateKSProb();
		fOutPValue             = fTemplateFitter->CalculatePValue();
		fOutNEvents            = int(fTemplateFitter->GetData().Integral());

		// fill the tree
		fTree->Fill();
	}

	// reset log level
	BCLog::SetLogLevel(ll);

	// no error
	return 1;
}

//---------------------------------------------------------------------------------------------------------
TH1D* BCTemplateEnsembleTest::BuildEnsemble()
{
	// get histogram parameters
  int nbins   = fEnsembleTemplate.GetNbinsX();
  double xmin = fEnsembleTemplate.GetXaxis()->GetXmin();
  double xmax = fEnsembleTemplate.GetXaxis()->GetXmax();

	// create new ensemble
  TH1D* ensemble = new TH1D("", "", nbins, xmin, xmax);

	// increase ensemble counter
	fEnsembleCounter++;

	// loop over bins and fill them
  for(int i = 1; i <= nbins; ++i){
    double p      = fEnsembleTemplate.GetBinContent(i);
    double lambda = p * fEnsembleExpectation;
		double n      = gRandom -> Poisson(lambda);

		// set the bin content
    ensemble->SetBinContent(i, n);
  }

	// return the ensemble histogram
  return ensemble;
}

//---------------------------------------------------------------------------------------------------------
int BCTemplateEnsembleTest::Write(const char * filename)
{
	// open file
	fFile = new TFile(filename, "RECREATE");

	// write tree
	fTree->Write();

	// close file
	fFile->Close();

	// free memory
	delete fFile;

	// no error
	return 1;
}

//---------------------------------------------------------------------------------------------------------
int BCTemplateEnsembleTest::PrepareTree()
{
	// delete old tree if necessary
	if (fTree)
		delete fTree;

	// create new tree
	fTree = new TTree("fTree", "fTree");

	// get number of parameters and ratios
	int npar = fTemplateFitter->GetNParameters();
	int nratios = fTemplateFitter->GetNRatios();

	// initialize variables
	fOutParModeGlobal.assign(npar, 0);
	fOutParErrorUpGlobal.assign(npar, 0);
	fOutParErrorDownGlobal.assign(npar, 0);
	fOutParModeMarg.assign(npar, 0);
	fOutParMeanMarg.assign(npar, 0);
	fOutParMedianMarg.assign(npar, 0);
	fOutParRMSMarg.assign(npar, 0);
	fOutParErrorUpMarg.assign(npar, 0);
	fOutParErrorDownMarg.assign(npar, 0);
	fOutParQuantile5Marg.assign(npar, 0);
	fOutParQuantile10Marg.assign(npar, 0);
	fOutParQuantile90Marg.assign(npar, 0);
	fOutParQuantile95Marg.assign(npar, 0);
	fOutRatioModeMarg.assign(nratios, 0);
	fOutRatioMeanMarg.assign(nratios, 0);
	fOutRatioMedianMarg.assign(nratios, 0);
	fOutRatioRMSMarg.assign(nratios, 0);
	fOutRatioErrorUpMarg.assign(nratios, 0);
	fOutRatioErrorDownMarg.assign(nratios, 0);
	fOutRatioQuantile5Marg.assign(nratios, 0);
	fOutRatioQuantile10Marg.assign(nratios, 0);
	fOutRatioQuantile90Marg.assign(nratios, 0);
	fOutRatioQuantile95Marg.assign(nratios, 0);

	fTree->Branch("chi2",     &fOutChi2,     "chi2/D");
	fTree->Branch("ndf",      &fOutNDF,      "ndf/I");
	fTree->Branch("chi2prob", &fOutChi2Prob, "chi2 prob probability/D");
	fTree->Branch("KSprob",   &fOutKSProb,   "KS probability/D");
	fTree->Branch("pvalue",   &fOutPValue,   "p-value/D");
	fTree->Branch("nevents",  &fOutNEvents,  "n events/I");

	for (int i = 0; i < npar; ++i) {
		// add branches
		fTree->Branch(Form("par_global_mode_par_%i", i),       &fOutParModeGlobal[i],      Form("par_global_Mode_par_%i/D", i));
		fTree->Branch(Form("par_global_error_up_par_%i", i),   &fOutParErrorUpGlobal[i],   Form("par_global_error_up_par_%i/D", i));
		fTree->Branch(Form("par_global_error_down_par_%i", i), &fOutParErrorDownGlobal[i], Form("par_global_error_down_par_%i/D", i));

		if(fFlagMCMC) {
			fTree->Branch(Form("par_marg_mode_par_%i", i),       &fOutParModeMarg[i],       Form("par_marg_mode_par_%i/D", i));
			fTree->Branch(Form("par_marg_mean_par_%i", i),       &fOutParMeanMarg[i],       Form("par_marg_mean_par_%i/D", i));
			fTree->Branch(Form("par_marg_median_par_%i", i),     &fOutParMedianMarg[i],     Form("par_marg_median_par_%i/D", i));
			fTree->Branch(Form("par_marg_rms_par_%i", i),        &fOutParRMSMarg[i],        Form("par_marg_rms_par_%i/D", i));
			fTree->Branch(Form("par_marg_error_up_par_%i", i),   &fOutParErrorUpMarg[i],    Form("par_marg_ErrorUp_par_%i/D", i));
			fTree->Branch(Form("par_marg_error_down_par_%i", i), &fOutParErrorDownMarg[i],  Form("par_marg_error_down_par_%i/D", i));
			fTree->Branch(Form("par_marg_quantile5_par_%i", i),  &fOutParQuantile5Marg[i],  Form("par_marg_Quantile5_par_%i/D", i));
			fTree->Branch(Form("par_marg_quantile10_par_%i", i), &fOutParQuantile10Marg[i], Form("par_marg_Quantile10_par_%i/D", i));
			fTree->Branch(Form("par_marg_quantile90_par_%i", i), &fOutParQuantile90Marg[i], Form("par_marg_Quantile90_par_%i/D", i));
			fTree->Branch(Form("par_marg_quantile95_par_%i", i), &fOutParQuantile95Marg[i], Form("par_marg_Quantile95_par_%i/D", i));
		}
	}

	if (fFlagMCMC) {
		for (int i = 0; i < nratios; ++i) {
			fTree->Branch(Form("ratio_marg_mode_ratio_%i", i),       &fOutRatioModeMarg[i],       Form("ratio_marg_mode_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_mean_ratio_%i", i),       &fOutRatioMeanMarg[i],       Form("ratio_marg_mean_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_median_ratio_%i", i),     &fOutRatioMedianMarg[i],     Form("ratio_marg_median_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_rms_ratio_%i", i),        &fOutRatioRMSMarg[i],        Form("ratio_marg_rms_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_error_up_ratio_%i", i),   &fOutRatioErrorUpMarg[i],    Form("ratio_marg_ErrorUp_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_error_down_ratio_%i", i), &fOutRatioErrorDownMarg[i],  Form("ratio_marg_error_down_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_quantile5_ratio_%i", i),  &fOutRatioQuantile5Marg[i],  Form("ratio_marg_Quantile5_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_quantile10_ratio_%i", i), &fOutRatioQuantile10Marg[i], Form("ratio_marg_Quantile10_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_quantile90_ratio_%i", i), &fOutRatioQuantile90Marg[i], Form("ratio_marg_Quantile90_ratio_%i/D", i));
			fTree->Branch(Form("ratio_marg_quantile95_ratio_%i", i), &fOutRatioQuantile95Marg[i], Form("ratio_marg_Quantile95_ratio_%i/D", i));
		}
	}

	// no error
	return 1;
}

//---------------------------------------------------------------------------------------------------------
