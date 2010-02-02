#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <../SummaryTool/SummaryTool.h>

#include <StackModel.h>
#include <StackModelManager.h>

int main()
{

	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("templates_gerda.root", "READ");

	TH1D hist_signal     = *((TH1D*) file->Get("hist_signal"));
	TH1D hist_background = *((TH1D*) file->Get("hist_background"));
	TH1D hist_data        = *((TH1D*) file->Get("hist_data"));
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

	// create new StackModel object
	StackModel * model = new StackModel("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 

	// set data histogram
	model->SetDataHistogram(hist_data);

	// calculate number of events
	double exposure = 100; 
	double bkgindex = 1e-3; 
	double nbkg = exposure * bkgindex * double(hist_sum.GetNbinsX());

	// add template histograms
	model->AddTemplateHistogram(hist_signal,     "Signal", 0.0,      nbkg+5*sqrt(nbkg)); 
	model->AddTemplateHistogram(hist_background, "Background", 0.0,  nbkg+5*sqrt(nbkg)); 

	// set efficiencies
	model->SetTemplateEfficiency(0, 1.0, -1.0);
	model->SetTemplateEfficiency(1, 1.0, -1.0);

	// set priors 
	model->SetTemplatePrior(1, nbkg, nbkg/2.0);

	// set constraints
	// ... no constraints 

	// ----------------------------------------------------
	// set up model manager
	// ----------------------------------------------------

 	// create new StackModelManager object
	StackModelManager * smm = new StackModelManager();

 	// add models
 	smm->AddStackModel(model);

 	// compare models
	smm->PerformAnalysis();

	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	// create summary tool
	SummaryTool* st = new SummaryTool(model); 

	// print results 
	smm->PrintResults("comparison.txt");

	// print data
 	TCanvas c1("c1");
	c1.cd();
	hist_data.Draw();
	c1.Print("data.ps");

	// print results
	model->PrintAllMarginalized("model_marginalized.eps"); 
	model->PrintStack("model_stack.eps");
	model->PrintFraction("model_fraction.ps");
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

	// delete model manager
	delete smm;

	// delete summary tool
	delete st; 

	return 0;
}

