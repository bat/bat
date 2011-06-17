#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCTemplateFitter.h>
#include <BAT/BCTemplateEnsembleTest.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TTree.h>

int main()
{
	// modify these settings
	bool flag_posterior = true; // create ensembles according to posterior (true)
                              // or according to the best fit parameter (false)

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

	// create new BCTemplateFitter object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set precision
	model->MCMCSetPrecision(BCEngineMCMC::kMedium);

	// set data histogram
	model->SetData(hist_data);

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
	// create output model
	// ----------------------------------------------------

	// create new output object
	BCModelOutput* mout = new BCModelOutput(model, "posterior.root");

	// switch writing of Markov Chains on
	mout->WriteMarkovChain(true);

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
	model->PrintAllMarginalized("model_marginalized.ps");
	model->PrintStack("model_stack.eps");
	model->PrintRatios("model_fraction.ps");
	model->PrintResults("model_results.txt");

	st->PrintParameterPlot("model_parameters.eps");
	st->PrintCorrelationPlot("model_correlation.eps");
	st->PrintKnowledgeUpdatePlots("model_update.ps");
	st->PrintCorrelationMatrix("model_matrix.eps");

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest();

	// settings
	tet->SetTemplateFitter(model);
	tet->SetTemplateParameters( model->GetBestFitParameters() );
	tet->SetFlagMCMC(false);
	tet->SetNEnsembles(100);
	tet->SetEnsembleExpectation(nbkg);

	// get tree from output file
	TTree* tree = (TTree*) mout->GetFile()->Get("MarkovChainTree_0");

	// perform ensemble tests with varying parameters for ensemble
	// generation. Parameters are varied according to posterior
	if (flag_posterior)
		tet->PerformEnsembleTest(tree);

	// perform ensemble tests with fixed parameter for ensemble
	// generation
	else
		tet->PerformEnsembleTest();

	// write results to file
	if (flag_posterior)
		tet->Write("ensembles_posterior.root");
	else
		tet->Write("ensembles_bestfit.root");

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// close output file
	mout->Close();

	// delete model
	delete model;

	// delete model output
	delete mout;

	// delete summary tool
	delete st;

	// delete ensemble test tool
	delete tet;

	return 0;
}

