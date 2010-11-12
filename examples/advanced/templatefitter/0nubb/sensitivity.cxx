#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCSummaryPriorModel.h>
#include <BAT/BCTemplateFitter.h>
#include <BAT/BCTemplateEnsembleTest.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>

int main()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("templates.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cerr << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	TH1D hist_signal     = *((TH1D*) file->Get("hist_sgn"));
	TH1D hist_background = *((TH1D*) file->Get("hist_bkg"));

	// close file
	file->Close();

	// ----------------------------------------------------
	// configure BAT
	// ----------------------------------------------------

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

 	// ----------------------------------------------------
	// Create new model
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set precision
	model->MCMCSetPrecision(BCEngineMCMC::kMedium);

	// set number of events
	double nbkg = 300.0; // prior assumption

	// add template histograms
	model->AddTemplate(hist_background, "Background", 100., 500.); 
	model->AddTemplate(hist_signal,     "Signal", 0., 200.); 

	// set efficiencies
	model->SetTemplateEfficiency("Signal",     1.0, 0.);
	model->SetTemplateEfficiency("Background", 1.0, 0.);

	// set priors 
	model->SetPriorGauss("Background", nbkg, nbkg/4.);
	model->SetPriorConstant("Signal");

	// set constraints
	// ... no constraints 

	// ----------------------------------------------------
	// create prior model
	// ----------------------------------------------------

	// create new prior model
	BCSummaryPriorModel* pmodel = new BCSummaryPriorModel();

	// set model
	pmodel->SetModel(model);

	// ----------------------------------------------------
	// create output for prior model
	// ----------------------------------------------------

  // create new output object
  BCModelOutput* mout = new BCModelOutput(pmodel, "prior.root");

  // switch writing of Markov Chains on
  mout->WriteMarkovChain(true);

	// ----------------------------------------------------
	// perform analysis
	// ----------------------------------------------------

	// run MCMC
	pmodel->MarginalizeAll(); 
	
	// find global mode
	pmodel->FindMode( pmodel->GetBestFitParameters() );

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest(); 

	// settings
	tet->SetTemplateFitter(model);
	tet->SetFlagMCMC(false);
	tet->SetNEnsembles(100); 
	tet->SetEnsembleExpectation(nbkg); 

	// get tree from output file
	TTree* tree = (TTree*) mout->GetFile()->Get("MarkovChainTree_0");

	// perform ensemble tests
	tet->PerformEnsembleTest(tree); 

	// print pulls
	tet->PrintPulls("pulls_prior.ps");

	// write results to file
	tet->Write("ensembles_prior.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

  // close output file
  mout->Close();

	// delete model
 	delete model;

	// delete prior model
	delete pmodel;

	// delete model output
	delete mout;
	
	// delete ensemble test tool
	delete tet;

	return 0;
}

