#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <TemplateModel.h>
#include <EnsembleTestTool.h>

int main()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("templates_gerda.root", "READ");

	TH1D hist_signal     = *((TH1D*) file->Get("hist_signal"));
	TH1D hist_background = *((TH1D*) file->Get("hist_background"));
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
	EnsembleTestTool* ett = new EnsembleTestTool(); 

	// calculate number of events
	double exposure = 3000; 
	double bkgindex = 1e-3; 
	double nbkg = exposure * bkgindex * double(hist_sum.GetNbinsX());

 	// ----------------------------------------------------
	// Create new model
	// ----------------------------------------------------

	// create new TemplateModel object
	TemplateModel * model = new TemplateModel("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 

	// set data histogram
	model->SetData(hist_sum);

	// add template histograms
	model->AddTemplate(hist_signal,     "Signal", 0.0,      nbkg+5*sqrt(nbkg)); 
	model->AddTemplate(hist_background, "Background", 0.0,  nbkg+5*sqrt(nbkg)); 

	// set efficiencies
	model->SetTemplateEfficiency("Signal", 1.0, 0.);
	model->SetTemplateEfficiency("Background", 1.0, 0.);

	// set priors 
	model->SetTemplatePrior("Background", nbkg, nbkg/2.0);

	// set constraints
	// ... no constraints 

	// settings
	ett->SetTemplateModel(model);
	ett->SetEnsembleTemplate(hist_sum);
	ett->SetNEnsembles(1000); 
	ett->SetEnsembleExpectation(nbkg); 
	ett->SetFlagMCMC(false); 

	// perform ensemble tests
	ett->PerformEnsembleTest(); 

	// write results to file
	ett->Write("ensemble_gerda.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
 	delete model;

	return 0;
}

