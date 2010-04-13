#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCSummaryTool.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <TemplateModel.h>

int main()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("templates_gerda.root", "READ");

	TH1D hist_signal     = *((TH1D*) file->Get("hist_signal"));
	TH1D hist_background = *((TH1D*) file->Get("hist_background"));
	TH1D hist_sum        = *((TH1D*) file->Get("hist_sum"));
	TH1D hist_data       = *((TH1D*) file->Get("hist_data"));

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

	// create new TemplateModel object
	TemplateModel * model = new TemplateModel("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 

	// set data histogram
	model->SetData(hist_data);

	// calculate number of events
	double exposure = 100; 
	double bkgindex = 1e-3; 
	double nbkg = exposure * bkgindex * double(hist_sum.GetNbinsX());

	// add template histograms
	model->AddTemplate(hist_signal,     "Signal", 0.0,      nbkg+5*sqrt(nbkg)); 
	model->AddTemplate(hist_background, "Background", 0.0,  nbkg+5*sqrt(nbkg)); 

	// set efficiencies
	model->SetTemplateEfficiency("Signal", 1., 0.);
	model->SetTemplateEfficiency("Background", 1., 0.);

	// set priors 
	model->SetTemplatePrior("Background", nbkg, nbkg/2.0);

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
	model->FindMode();


	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	// create summary tool
	BCSummaryTool* st = new BCSummaryTool(model); 

	// print data
 	TCanvas c1("c1");
	c1.cd();
	hist_data.Draw();
	c1.Print("data.ps");

	// print results
	model->PrintAllMarginalized("model_marginalized.eps"); 
	model->PrintStack("model_stack.eps");
	model->PrintRatios("model_fraction.ps");
	model->PrintResults("model_results.txt"); 

	st->PrintParameterPlot("model_parameters.eps"); 
	st->PrintCorrelationPlot("model_correlation.eps"); 
	st->PrintKnowlegdeUpdatePlot("model_update.eps"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete models
 	delete model;

	// delete summary tool
	delete st; 

	return 0;
}

