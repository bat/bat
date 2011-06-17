#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCTemplateFitter.h>

#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
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
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	TH1D hist_signal     = *((TH1D*) file->Get("hist_sgn"));
	TH1D hist_background = *((TH1D*) file->Get("hist_bkg"));
	TH1D hist_sum        = *((TH1D*) file->Get("hist_sum"));
	TH1D hist_data       = *((TH1D*) file->Get("hist_data"));

	TF1 func_signal = TF1("signal", "1./sqrt((2.*3.14159265358979312)*[1])*exp(- (x-[0])*(x-[0])/2./[1]/[1])", 2039. - 50., 2039. + 50.);
	func_signal.SetParName(0, "Signal mass");
	func_signal.SetParName(1, "Signal width");
	func_signal.SetParLimits(0, 2039.-5., 2039.+5.);
	func_signal.SetParLimits(1, 0.5, 10.0);

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
	model->SetData(hist_data);

	// calculate number of events
	double nbkg = 300.0; // background assumption

	// add template histograms
	model->AddTemplate(hist_background, "Background", 200., 400.);
	model->AddTemplate(func_signal, "Signal",     0., 200. );

	// set efficiencies
	model->SetTemplateEfficiency("Signal",     1., 0.);
	model->SetTemplateEfficiency("Background", 1., 0.);

	// set priors
	model->SetPriorConstant("Signal");
	model->SetPriorGauss("Background", nbkg, nbkg/4.);
	model->SetPriorConstant("Signal mass");
	model->SetPriorGauss("Signal width", 5.0, 1.0);

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
	model->PrintAllMarginalized("model_marginalized.ps");
	model->PrintStack("model_stack.eps");

	//	model->PrintRatios("model_fraction.ps");
	model->PrintResults("model_results.txt");

	st->PrintParameterPlot("model_parameters.eps");
	st->PrintCorrelationPlot("model_correlation.eps");
	st->PrintKnowledgeUpdatePlots("model_update.ps");
	st->PrintCorrelationMatrix("model_matrix.eps");

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

