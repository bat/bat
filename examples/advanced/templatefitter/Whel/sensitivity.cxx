#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCTemplateFitter.h>
#include <BAT/BCTemplateEnsembleTest.h>
#include <BAT/BCSummaryPriorModel.h>

int main()
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

	// templates
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

	// systematic uncertainty 1 for all four contributions
	TH1D hist_systerror1_bkg = *((TH1D*) file->Get("hist_systerror1_bkg"));
	TH1D hist_systerror1_sgn_h0 = *((TH1D*) file->Get("hist_systerror1_sgn_h0"));
	TH1D hist_systerror1_sgn_hL = *((TH1D*) file->Get("hist_systerror1_sgn_hL"));
	TH1D hist_systerror1_sgn_hR = *((TH1D*) file->Get("hist_systerror1_sgn_hR"));

	// systematic uncertainty 2 for all four contributions
	TH1D hist_systerror2_bkg = *((TH1D*) file->Get("hist_systerror2_bkg"));
	TH1D hist_systerror2_sgn_h0 = *((TH1D*) file->Get("hist_systerror2_sgn_h0"));
	TH1D hist_systerror2_sgn_hL = *((TH1D*) file->Get("hist_systerror2_sgn_hL"));
	TH1D hist_systerror2_sgn_hR = *((TH1D*) file->Get("hist_systerror2_sgn_hR"));

	// systematic uncertainty 3 for all four contributions
	TH1D hist_systerror3_bkg = *((TH1D*) file->Get("hist_systerror3_bkg"));
	TH1D hist_systerror3_sgn_h0 = *((TH1D*) file->Get("hist_systerror3_sgn_h0"));
	TH1D hist_systerror3_sgn_hL = *((TH1D*) file->Get("hist_systerror3_sgn_hL"));
	TH1D hist_systerror3_sgn_hR = *((TH1D*) file->Get("hist_systerror3_sgn_hR"));

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
	// create model
	// ----------------------------------------------------

	// create new TemplateModel object
	BCTemplateFitter * model = new BCTemplateFitter("model");

	// set precision
	model->MCMCSetPrecision(BCEngineMCMC::kMedium);

	// set template fitting options
	model->SetFlagPhysicalLimits(true);

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

	// set priors
	model->SetPriorGauss("background", 1300000., 2000.);
	model->SetPriorConstant("signal (h= 0)");
	model->SetPriorConstant("signal (h=-1)");
	model->SetPriorConstant("signal (h=+1)");

	// add systematic uncertainty 1
	model->AddSystError("syst. error 1", "gauss");
	model->SetTemplateSystError("syst. error 1", "background",    hist_systerror1_bkg);
	model->SetTemplateSystError("syst. error 1", "signal (h= 0)", hist_systerror1_sgn_h0);
	model->SetTemplateSystError("syst. error 1", "signal (h=-1)", hist_systerror1_sgn_hL);
	model->SetTemplateSystError("syst. error 1", "signal (h=+1)", hist_systerror1_sgn_hR);

	// add systematic uncertainty 2
	model->AddSystError("syst. error 2", "gauss");
	model->SetTemplateSystError("syst. error 2", "background",    hist_systerror2_bkg);
	model->SetTemplateSystError("syst. error 2", "signal (h= 0)", hist_systerror2_sgn_h0);
	model->SetTemplateSystError("syst. error 2", "signal (h=-1)", hist_systerror2_sgn_hL);
	model->SetTemplateSystError("syst. error 2", "signal (h=+1)", hist_systerror2_sgn_hR);

	// add systematic uncertainty 3
	model->AddSystError("syst. error 3", "gauss");
	model->SetTemplateSystError("syst. error 3", "background",    hist_systerror3_bkg);
	model->SetTemplateSystError("syst. error 3", "signal (h= 0)", hist_systerror3_sgn_h0);
	model->SetTemplateSystError("syst. error 3", "signal (h=-1)", hist_systerror3_sgn_hL);
	model->SetTemplateSystError("syst. error 3", "signal (h=+1)", hist_systerror3_sgn_hR);

	// ----------------------------------------------------
	// create prior model
	// ----------------------------------------------------

	// create new prior model
	BCSummaryPriorModel* pmodel = new BCSummaryPriorModel();

	// set model
	pmodel->SetModel(model);

	// ----------------------------------------------------
	// create output model
	// ----------------------------------------------------

	// create new output object
	BCModelOutput* mout = new BCModelOutput(pmodel, "prior.root");

	// switch writing of Markov Chains on
	mout->WriteMarkovChain(true);

	// ----------------------------------------------------
	// perform analysis
	// ----------------------------------------------------

	// run MCMC
	pmodel->MarginalizeAll();

	// find global mode
	pmodel->FindMode( pmodel->GetBestFitParameters() );

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
	tet->SetEnsembleExpectation(3151);

	// get tree from output file
	TTree* tree = (TTree*) mout->GetFile()->Get("MarkovChainTree_0");

	// perform ensemble tests
	tet->PerformEnsembleTest(tree);

	// print pulls
	tet->PrintPulls("pulls_prior.ps");

	// write results to file
	tet->Write("ensembles_prior.root");

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// close output file
	mout->Close();

	// delete model
	delete model;

	// delete prior model
	delete pmodel;

	// delete model output
	delete mout;

	// delete ensemble test tool
	delete tet;

	// no error
	return 0;
}

