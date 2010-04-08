#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <StackModel.h>
#include <StackModelManager.h>

int main()
{
	// ----------------------------------------------------
	// open file with data and templates
	// ----------------------------------------------------

	TFile * file = new TFile("histograms.root", "READ");

	TH1D hist_bkg1 = *((TH1D*) file->Get("hist_bkg1"));
	TH1D hist_bkg2 = *((TH1D*) file->Get("hist_bkg2"));
	TH1D hist_sgn  = *((TH1D*) file->Get("hist_sgn"));
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
	// model: 1 bkg 0 sgn
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * m_1bkg0sgn = new StackModel("model 1 bkg 0 sgn");

	// set data histogram
	m_1bkg0sgn->SetDataHistogram(hist_data);

	// add template histograms
	m_1bkg0sgn->AddTemplateHistogram(hist_bkg1, "bkg 1");

	// ----------------------------------------------------
	// model: 1 bkg 1 sgn
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * m_1bkg1sgn = new StackModel("model 1 bkg 1 sgn");

	// set data histogram
	m_1bkg1sgn->SetDataHistogram(hist_data);

	// add template histograms
	m_1bkg1sgn->AddTemplateHistogram(hist_bkg1, "bkg 1");
	m_1bkg1sgn->AddTemplateHistogram(hist_sgn, "sgn", 0., 50.);

	// ----------------------------------------------------
	// model: 2 bkg 0 sgn
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * m_2bkg0sgn = new StackModel("model 2 bkg 0 sgn");

	// set data histogram
	m_2bkg0sgn->SetDataHistogram(hist_data);

	// add template histograms
	m_2bkg0sgn->AddTemplateHistogram(hist_bkg1, "bkg 1");
	m_2bkg0sgn->AddTemplateHistogram(hist_bkg2, "bkg 2");
 
	// ----------------------------------------------------
	// model: 2 bkg 1 sgn
	// ----------------------------------------------------

	// create new StackModel object
	StackModel * m_2bkg1sgn = new StackModel("model 2 bkg 1 sgn");
	m_2bkg1sgn->MCMCSetNLag(10); 

	// set data histogram
	m_2bkg1sgn->SetDataHistogram(hist_data);

	// add template histograms
	m_2bkg1sgn->AddTemplateHistogram(hist_bkg1, "bkg 1");
	m_2bkg1sgn->AddTemplateHistogram(hist_bkg2, "bkg 2");
	m_2bkg1sgn->AddTemplateHistogram(hist_sgn, "sgn", 0., 50.);

	// ----------------------------------------------------
	// set up model manager
	// ----------------------------------------------------

 	// create new StackModelManager object
	StackModelManager * smm = new StackModelManager();

 	// add models
 	smm->AddStackModel(m_1bkg0sgn);
 	smm->AddStackModel(m_1bkg1sgn);
 	smm->AddStackModel(m_2bkg0sgn);
 	smm->AddStackModel(m_2bkg1sgn);

	smm->SetFlagFixNorm(false); 

 	// compare models
	smm->PerformAnalysis();

	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	// print results 
	smm->PrintResults("comparison.txt");

	// print data
 	TCanvas c1("c1");
	c1.cd();
	hist_data.Draw();
	c1.Print("data.ps");

	// print marginalized distribution for signal
 	c1.cd();
 	m_1bkg1sgn->GetMarginalized("N_1")->Draw(0, -95);
 	c1.Print("1bkg1sgn_sgn.ps");

	c1.cd();
	m_2bkg1sgn->GetMarginalized("N_2")->Draw(0, -95);
	c1.Print("2bkg1sgn_sgn.ps");

	// print stack plots and results
 	m_1bkg0sgn->PrintStack("1bkg0sgn_stack.eps");
 	m_1bkg1sgn->PrintStack("1bkg1sgn_stack.eps");
 	m_2bkg0sgn->PrintStack("2bkg0sgn_stack.eps");
	m_2bkg1sgn->PrintStack("2bkg1sgn_stack.eps");

 	m_1bkg0sgn->PrintRatios("1bkg0sgn_fraction.ps");
 	m_1bkg1sgn->PrintRatios("1bkg1sgn_fraction.ps");
 	m_2bkg0sgn->PrintRatios("2bkg0sgn_fraction.ps");
	m_2bkg1sgn->PrintRatios("2bkg1sgn_fraction.ps");

 	m_1bkg0sgn->PrintResults("1bkg0sgn_results.txt"); 
 	m_1bkg1sgn->PrintResults("1bkg1sgn_results.txt"); 
 	m_2bkg0sgn->PrintResults("2bkg0sgn_results.txt"); 
	m_2bkg1sgn->PrintResults("2bkg1sgn_results.txt"); 

	// ----------------------------------------------------
	// clean-up and end
	// ----------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// delete models
 	delete m_1bkg0sgn;
 	delete m_1bkg1sgn;
 	delete m_2bkg0sgn;
 	delete m_2bkg1sgn;

	// delete model manager
	delete smm;

	return 0;
}

