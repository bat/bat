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

	TFile * file = new TFile("templates.root", "READ");

	TH1D hist_process1 = *((TH1D*) file->Get("hist_process1"));
	TH1D hist_process2 = *((TH1D*) file->Get("hist_process2"));
	TH1D hist_process3 = *((TH1D*) file->Get("hist_process3"));
	TH1D hist_process4 = *((TH1D*) file->Get("hist_process4"));
	TH1D hist_sum = *((TH1D*) file->Get("hist_sum"));

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

	// set options
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 
	model->SetFlagPhysicalLimits(true);

	// set data histogram
	model->SetData(hist_sum);

	// add template histograms
	model->AddTemplate(hist_process1, "background",    0.0, 2000000.0);
	model->AddTemplate(hist_process2, "signal (h= 0)", 0.0,    8000.0);
	model->AddTemplate(hist_process3, "signal (h=-1)", 0.0,    8000.0);
	model->AddTemplate(hist_process4, "signal (h=+1)", 0.0,    8000.0);

	// set efficiencies
	model->SetTemplateEfficiency("background", 0.001, 0.0005);
	model->SetTemplateEfficiency("signal (h= 0)", 0.20, 0.05);
	model->SetTemplateEfficiency("signal (h=-1)", 0.20, 0.05);
	model->SetTemplateEfficiency("signal (h=+1)", 0.20, 0.05);

	// set priors 
	model->SetTemplatePrior("background", 1300000.0, 2000.0);

	// set constraints
	std::vector <int> indices; 
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(3);
	model->ConstrainSum(indices, 7000.0, 100); 

	// set up printing of fractions
	model->CalculateRatio(1, indices); 
	model->CalculateRatio(2, indices); 
	model->CalculateRatio(3, indices); 

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	EnsembleTestTool* ett = new EnsembleTestTool(); 

	// settings
	ett->SetTemplateModel(model);
	ett->SetEnsembleTemplate(hist_sum);
	ett->SetNEnsembles(1000); 
	ett->SetEnsembleExpectation(2700); 
	ett->SetFlagMCMC(true); 

	// perform ensemble tests
	ett->PerformEnsembleTest(); 

	// write results to file
	ett->Write("ensemble.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
 	delete model;

	return 0;
}

