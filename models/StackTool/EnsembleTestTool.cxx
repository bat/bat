#include <EnsembleTestTool.h>
#include <StackModel.h>

#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom3.h>

#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

// ---------------------------------------------------------
EnsembleTestTool::EnsembleTestTool() : fStackModel(0)
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
EnsembleTestTool::~EnsembleTestTool()
{
;
}
 
// ---------------------------------------------------------
int EnsembleTestTool::SetEnsembleTemplate(TH1D hist)
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
int EnsembleTestTool::PerformEnsembleTest()
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
      cout << "Fraction of ensembles analyzed: " << double(j-1) / double(fNEnsembles) * 100 << "%" << std::endl;
		
		// create new ensemble 
		TH1D ensemble = BuildEnsemble();

		// set ensemble as new data set
		fStackModel->SetDataHistogram(ensemble);

		// find mode
		fStackModel->FindMode();

		// perform MCMC
		if(fFlagMCMC) {
			fStackModel->MarginalizeAll();

			// get number of parameters
			int npar = fStackModel->GetNParameters();
			
			// loop over parameters and set tree variables
			for (int i = 0; i < npar; ++i) {
				BCH1D* hist = fStackModel->GetMarginalized(fStackModel->GetParameter(i)); 
				fOutModeMarg[i]       = hist->GetMode();
				fOutMedianMarg[i]     = hist->GetMedian();
				fOutMeanMarg[i]       = hist->GetMean();
				fOutRMSMarg[i]        = hist->GetRMS();
				fOutErrorUpMarg[i]    = hist->GetQuantile(0.84)-hist->GetMode();
				fOutErrorDownMarg[i]  = hist->GetMode()-hist->GetQuantile(0.16);
				fOutQuantile5Marg[i]  = hist->GetQuantile(0.05);
				fOutQuantile10Marg[i] = hist->GetQuantile(0.10);
				fOutQuantile90Marg[i] = hist->GetQuantile(0.90);
				fOutQuantile95Marg[i] = hist->GetQuantile(0.95);
			}
		}

		// set tree variables
		fOutModeGlobal      = fStackModel->GetBestFitParameters();
		fOutErrorUpGlobal   = fStackModel->GetBestFitParameterErrors();
		fOutErrorDownGlobal = fStackModel->GetBestFitParameterErrors();
		fOutChi2            = fStackModel->CalculateChi2(); 
		fOutNDF             = fStackModel->GetNDF();
		fOutChi2Prob        = fStackModel->CalculateChi2Prob(); 
		fOutKSProb          = fStackModel->CalculateKSProb(); 
		fOutPValue          = fStackModel->CalculatePValue(); 

		// fill the tree 
		fTree->Fill();
	}

	// reset log level
	BCLog::SetLogLevel(ll);

	// no error
	return 1;
}

//---------------------------------------------------------------------------------------------------------
TH1D EnsembleTestTool::BuildEnsemble()
{
	// get histogram parameters 
  int nbins   = fEnsembleTemplate.GetNbinsX();
  double xmin = fEnsembleTemplate.GetXaxis()->GetXmin();
  double xmax = fEnsembleTemplate.GetXaxis()->GetXmax();

	// create new ensemble 
  TH1D ensemble = TH1D(Form("ensemble_%i", fEnsembleCounter), "", nbins, xmin, xmax);

	// increase ensemble counter 
	fEnsembleCounter++; 

	// loop over bins and fill them 
  for(int i = 1; i <= nbins; ++i){
    double p      = fEnsembleTemplate.GetBinContent(i);
    double lambda = p * fEnsembleExpectation;
		double n      = gRandom -> Poisson(lambda);

		// set the bin content
    ensemble.SetBinContent(i, n);
  }

	// return the ensemble histogram
  return ensemble;
}

//---------------------------------------------------------------------------------------------------------
int EnsembleTestTool::Write(const char* filename)
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
int EnsembleTestTool::PrepareTree()
{
	// delete old tree if necessary
	if (fTree)
		delete fTree; 

	// create new tree 
	fTree = new TTree("fTree", "fTree"); 

	// get number of parameters
	int npar = fStackModel->GetNParameters(); 
	
	// initialize variables
	fOutModeGlobal.assign(npar, 0);
	fOutErrorUpGlobal.assign(npar, 0);
	fOutErrorDownGlobal.assign(npar, 0);
	fOutModeMarg.assign(npar, 0);
	fOutMeanMarg.assign(npar, 0);
	fOutMedianMarg.assign(npar, 0);
	fOutRMSMarg.assign(npar, 0);
	fOutErrorUpMarg.assign(npar, 0);
	fOutErrorDownMarg.assign(npar, 0);
	fOutQuantile5Marg.assign(npar, 0);
	fOutQuantile10Marg.assign(npar, 0);
	fOutQuantile90Marg.assign(npar, 0);
	fOutQuantile95Marg.assign(npar, 0);

	fTree->Branch("chi2",
								&fOutChi2, 
								"chi2/D");
	
	fTree->Branch("ndf",
								&fOutNDF, 
								"ndf/I");
	
	fTree->Branch("chi2prob",
								&fOutChi2Prob, 
								"chi2 prob probability/D");
	
	fTree->Branch("KSprob",
								&fOutKSProb, 
								"KS probability/D");
	
	fTree->Branch("pvalue",
								&fOutPValue, 
								"p-value/D");
	
	for (int i = 0; i < npar; ++i) {
		// add branches 
		fTree->Branch(Form("global_mode_par_%i", i), 
									&fOutModeGlobal[i], 
									Form("global_Mode_par_%i/D", i)); 
		
		fTree->Branch(Form("global_error_up_par_%i", i), 
									&fOutErrorUpGlobal[i], 
									Form("global_error_up_par_%i/D", i)); 
		
		fTree->Branch(Form("global_error_down_par_%i", i), 
									&fOutErrorDownGlobal[i], 
									Form("global_error_down_par_%i/D", i)); 

		if(fFlagMCMC) {
			fTree->Branch(Form("marg_mode_par_%i", i), 
										&fOutModeMarg[i], 
										Form("marg_mode_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_mean_par_%i", i), 
										&fOutMeanMarg[i], 
										Form("marg_mean_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_median_par_%i", i), 
										&fOutMedianMarg[i], 
										Form("marg_median_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_rms_par_%i", i), 
										&fOutRMSMarg[i], 
										Form("marg_rms_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_error_up_par_%i", i), 
										&fOutErrorUpMarg[i], 
										Form("marg_ErrorUp_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_error_down_par_%i", i), 
										&fOutErrorDownMarg[i], 
										Form("marg_error_down_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_quantile5_par_%i", i), 
										&fOutQuantile5Marg[i], 
										Form("marg_Quantile5_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_quantile10_par_%i", i), 
										&fOutQuantile10Marg[i], 
										Form("marg_Quantile10_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_quantile90_par_%i", i), 
										&fOutQuantile90Marg[i], 
										Form("marg_Quantile90_par_%i/D", i)); 
			
			fTree->Branch(Form("marg_quantile95_par_%i", i), 
										&fOutQuantile95Marg[i], 
										Form("marg_Quantile95_par_%i/D", i)); 
		}
		
	}
	
	// no error 
	return 1;
}

//---------------------------------------------------------------------------------------------------------
