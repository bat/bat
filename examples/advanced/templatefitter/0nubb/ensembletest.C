#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCTemplateFitter.h>
#include <BAT/BCTemplateEnsembleTest.h>

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
		std::cerr << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
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
	// Create new model
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set precision
	model->MCMCSetPrecision(BCEngineMCMC::kMedium);

	// set data histogram
	model->SetData(hist_sum);

	// set number of events
	double nbkg = 300.0; // prior assumption

	// add template histograms
	model->AddTemplate(hist_background, "Background", 100., 400.); 
	model->AddTemplate(hist_signal,     "Signal", 50., 150.); 

	// set efficiencies
	model->SetTemplateEfficiency("Signal",     1.0, 0.);
	model->SetTemplateEfficiency("Background", 1.0, 0.);

	// set priors 
	model->SetTemplatePrior("Background", nbkg, 2.*sqrt(nbkg));

	// set constraints
	// ... no constraints 

	// ----------------------------------------------------
	// perform analysis
	// ----------------------------------------------------

	// initialize model
	model->Initialize();

	// run MCMC
	model->MarginalizeAll(); 
	
	// find global mode
	model->FindMode( model->GetBestFitParameters() );

	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	// create summary tool
	BCSummaryTool* st = new BCSummaryTool(model); 

	// print results
	model->PrintAllMarginalized("model_marginalized.eps"); 
	model->PrintStack("model_stack.eps");
	model->PrintRatios("model_fraction.ps");
	model->PrintResults("model_results.txt"); 

	st->PrintParameterPlot("model_parameters.eps"); 
	st->PrintCorrelationPlot("model_correlation.eps"); 
	st->PrintKnowledgeUpdatePlots("model_update.eps"); 
	st->PrintCorrelationMatrix("model_matrix.eps");

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest(); 

	// settings
	tet->SetTemplateFitter(model);
	tet->SetTemplateParameters( model->GetBestFitParameters() );
	tet->SetFlagMCMC(true);
	tet->SetNEnsembles(100); 
	tet->SetEnsembleExpectation(nbkg); 

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

	// delete summary tool
	delete st; 

	// delete ensemble test tool
	delete tet;

	return 0;
}

