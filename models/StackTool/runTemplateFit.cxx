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
	// create model
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * model = new StackModel("model");
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000000); 

	// set data histogram
	model->SetDataHistogram(hist_data);

	// add template histograms
	model->AddTemplateHistogram(hist_process1, "Background",    0.0, 2000000.0);
	model->AddTemplateHistogram(hist_process2, "Signal (h= 0)", 0.0,   20000.0);
	model->AddTemplateHistogram(hist_process3, "Signal (h=-1)", 0.0,   20000.0);
	model->AddTemplateHistogram(hist_process4, "Signal (h=+1)", 0.0,   20000.0);

	// set efficiencies
	model->SetTemplateEfficiency(0, 0.001, 0.0005, true);
	model->SetTemplateEfficiency(1, 0.20, 0.05, true);
	model->SetTemplateEfficiency(2, 0.20, 0.05, true);
	model->SetTemplateEfficiency(3, 0.20, 0.05, true);

	model->SetTemplatePrior(0, 1300000.0, 50000.0);

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

