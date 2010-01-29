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

	TFile * file = new TFile("templates.root", "READ");

	TH1D hist_process1 = *((TH1D*) file->Get("hist_process1"));
	TH1D hist_process2 = *((TH1D*) file->Get("hist_process2"));
	TH1D hist_process3 = *((TH1D*) file->Get("hist_process3"));
	TH1D hist_process4 = *((TH1D*) file->Get("hist_process4"));
	TH1D hist_process5 = *((TH1D*) file->Get("hist_process5"));
	TH1D hist_prior_process1 = *((TH1D*) file->Get("hist_prior_process1"));
	TH1D hist_prior_process2 = *((TH1D*) file->Get("hist_prior_process2"));
	TH1D hist_prior_process3 = *((TH1D*) file->Get("hist_prior_process3"));
	TH1D hist_prior_process4 = *((TH1D*) file->Get("hist_prior_process4"));
	TH1D hist_prior_process5 = *((TH1D*) file->Get("hist_prior_process5"));
	TH1D hist_data = *((TH1D*) file->Get("hist_data"));

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
	// model: 2 bkg 1 sgn
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * model = new StackModel("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(1000000); 

	// set data histogram
	model->SetDataHistogram(hist_data);

	// add template histograms
	model->AddTemplateHistogram(hist_process1, "process 1", hist_prior_process1);
	model->AddTemplateHistogram(hist_process2, "process 2", hist_prior_process2);
	model->AddTemplateHistogram(hist_process3, "process 3", hist_prior_process3);
	model->AddTemplateHistogram(hist_process4, "process 4", hist_prior_process4);
	model->AddTemplateHistogram(hist_process5, "process 5", hist_prior_process5);

	// set efficiencies
	model->SetTemplateEfficiency(0, 0.85, 0.03);
	model->SetTemplateEfficiency(1, 0.52, 0.05);
	model->SetTemplateEfficiency(2, 0.79, 0.05);
	model->SetTemplateEfficiency(3, 0.75, 0.05);
	model->SetTemplateEfficiency(4, 0.27, 0.05);

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

	return 0;
}

