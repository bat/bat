#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <models/BCTemplateFitter.h>
#include <models/BCTemplateEnsembleTest.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

int main()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("templates.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cout << "Could not open file. Exit." << std::endl;
		return 1;
	}

	TH1D hist_signal     = *((TH1D*) file->Get("hist_sgn"));
	TH1D hist_background = *((TH1D*) file->Get("hist_bkg"));
	TH1D hist_sum        = *((TH1D*) file->Get("hist_sum"));

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
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest(); 

	// calculate number of events
	double nbkg = 100.0; // prior assumption

 	// ----------------------------------------------------
	// Create new model
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * model = new BCTemplateFitter("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 

	// set data histogram
	model->SetData(hist_sum);

	// add template histograms
	model->AddTemplate(hist_background, "Background", 0.0,  2.5*nbkg); 
	model->AddTemplate(hist_signal,     "Signal", 0.0,      nbkg+2.5*nbkg); 

	// set efficiencies
	model->SetTemplateEfficiency("Signal", 1.0, 0.);
	model->SetTemplateEfficiency("Background", 1.0, 0.);

	// set priors 
	model->SetTemplatePrior("Background", nbkg, nbkg/2.0);

	// set constraints
	// ... no constraints 

	// stetings
	tet->SetTemplateFitter(model);
	tet->SetEnsembleTemplate(hist_sum);
	tet->SetNEnsembles(1000); 
	tet->SetEnsembleExpectation(nbkg); 
	tet->SetFlagMCMC(false); 

	// perform ensemble tests
	tet->PerformEnsembleTest(); 

	// write results to file
	tet->Write("ensembles.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
 	delete model;

	return 0;
}

