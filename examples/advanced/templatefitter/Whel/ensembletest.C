#include <TDirectory.h>
#include <TFile.h> 
#include <TH1D.h> 
#include <TCanvas.h> 

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCTemplateFitter.h>
#include <BAT/BCTemplateEnsembleTest.h>
#include <BAT/BCSummaryTool.h>

int ensembletest()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	// remember old directory
	TDirectory* f = gDirectory;

	// open file
	TFile * file = new TFile("templates.root", "READ");

	// check if file is open
	if (!file->IsOpen()) {
		std::cout << "Could not open file. Exit." << std::endl;
		std::cerr << "Run macro CreateHistograms.C in Root to create the file." << std::endl;
		return 1;
	}

	// go back to old directory for memory handling
	f->cd();

	// get templates
	TH1D hist_bkg = *((TH1D*) file->Get("hist_bkg"));
	TH1D hist_sgn_h0 = *((TH1D*) file->Get("hist_sgn_h0"));
	TH1D hist_sgn_hL = *((TH1D*) file->Get("hist_sgn_hL"));
	TH1D hist_sgn_hR = *((TH1D*) file->Get("hist_sgn_hR"));

	// efficiency parameterization
	TH1D hist_efficiency_bkg = *((TH1D*) file->Get("hist_efficiency_bkg")); 
	TH1D hist_efficiency_sgn_h0 = *((TH1D*) file->Get("hist_efficiency_sgn_h0")); 
	TH1D hist_efficiency_sgn_hL = *((TH1D*) file->Get("hist_efficiency_sgn_hL")); 
	TH1D hist_efficiency_sgn_hR = *((TH1D*) file->Get("hist_efficiency_sgn_hR")); 

	// efficiency uncertainty parameterization
	TH1D hist_efferror_bkg = *((TH1D*) file->Get("hist_efferror_bkg")); 
	TH1D hist_efferror_sgn_h0 = *((TH1D*) file->Get("hist_efferror_sgn_h0")); 
	TH1D hist_efferror_sgn_hL = *((TH1D*) file->Get("hist_efferror_sgn_hL")); 
	TH1D hist_efferror_sgn_hR = *((TH1D*) file->Get("hist_efferror_sgn_hR")); 	

	// data template
	TH1D hist_sum = *((TH1D*) file->Get("hist_sum"));

	// close file
	file->Close();

	// delete file
	delete file;

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

	// set options
	model->MCMCSetNLag(10); 
	model->MCMCSetNIterationsRun(10000); 
	model->SetFlagPhysicalLimits(true);

	// set data histogram
	model->SetData(hist_sum);

	// add template histograms
	model->AddTemplate(hist_bkg, "background",       1300000.-10000., 1300000.+10000.);
	model->AddTemplate(hist_sgn_h0, "signal (h= 0)", 0.0,               10000.);
	model->AddTemplate(hist_sgn_hL, "signal (h=-1)", 0.0,               10000.);
	model->AddTemplate(hist_sgn_hR, "signal (h=+1)", 0.0,               10000.);

	// set efficiencies and uncertainties
	model->SetTemplateEfficiency("background",    hist_efficiency_bkg,    hist_efferror_bkg);
	model->SetTemplateEfficiency("signal (h= 0)", hist_efficiency_sgn_h0, hist_efferror_sgn_h0);
	model->SetTemplateEfficiency("signal (h=-1)", hist_efficiency_sgn_hL, hist_efferror_sgn_hL);
	model->SetTemplateEfficiency("signal (h=+1)", hist_efficiency_sgn_hR, hist_efferror_sgn_hR);

	// set prior on background
	model->SetTemplatePrior("background", 1300000., 2000.);

	// set prior on signal contribution
	std::vector<int> indices; 
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(3);
	model->ConstrainSum(indices, 10000., 100.); 

	// ----------------------------------------------------
	// create ensemble test tool
	// ----------------------------------------------------

	// create ensemble test tool
	BCTemplateEnsembleTest* tet = new BCTemplateEnsembleTest(); 

	// stetings
	tet->SetTemplateFitter(model);
	tet->SetEnsembleTemplate(hist_sum); // this histogram is used as a template for the pseudo data
	tet->SetNEnsembles(1000); 
	tet->SetEnsembleExpectation(3251); 

	// perform ensemble tests
	tet->PerformEnsembleTest(); 

	// write results to file
	tet->Write("ensemble.root"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete model
 	delete model;

	return 0;
}

