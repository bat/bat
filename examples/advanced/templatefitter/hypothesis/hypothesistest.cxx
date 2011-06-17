#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCTemplateFitter.h>

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>

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

	// create new BCTemplateFitter object
	BCTemplateFitter * m_1bkg0sgn = new BCTemplateFitter("model 1 bkg 0 sgn");

	// set data histogram
	m_1bkg0sgn->SetData(hist_data);

	// add template histograms
	m_1bkg0sgn->AddTemplate(hist_bkg1, "bkg 1",0., 100.0);

	// set prior
	m_1bkg0sgn->SetPriorConstant("bkg 1");

	// ----------------------------------------------------
	// model: 1 bkg 1 sgn
	// ----------------------------------------------------

	// create new StackModel object
	BCTemplateFitter * m_1bkg1sgn = new BCTemplateFitter("model 1 bkg 1 sgn");

	// set data histogram
	m_1bkg1sgn->SetData(hist_data);

	// add template histograms
	m_1bkg1sgn->AddTemplate(hist_bkg1, "bkg 1", 0., 100.0);
	m_1bkg1sgn->AddTemplate(hist_sgn, "sgn", 0., 50.);

	// set priors
	m_1bkg1sgn->SetPriorConstant("bkg 1");
	m_1bkg1sgn->SetPriorConstant("sgn");

	// ----------------------------------------------------
	// model: 2 bkg 0 sgn
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * m_2bkg0sgn = new BCTemplateFitter("model 2 bkg 0 sgn");

	// set data histogram
	m_2bkg0sgn->SetData(hist_data);

	// add template histograms
	m_2bkg0sgn->AddTemplate(hist_bkg1, "bkg 1", 0., 100.0);
	m_2bkg0sgn->AddTemplate(hist_bkg2, "bkg 2", 0., 100.0);

	// set priors
	m_2bkg0sgn->SetPriorConstant("bkg 1");
	m_2bkg0sgn->SetPriorConstant("bkg 2");

	// ----------------------------------------------------
	// model: 2 bkg 1 sgn
	// ----------------------------------------------------

	// create new BCTemplateFitter object
	BCTemplateFitter * m_2bkg1sgn = new BCTemplateFitter("model 2 bkg 1 sgn");
	m_2bkg1sgn->MCMCSetNLag(10); 

	// set data histogram
	m_2bkg1sgn->SetData(hist_data);

	// add template histograms
	m_2bkg1sgn->AddTemplate(hist_bkg1, "bkg 1", 0., 100.);
	m_2bkg1sgn->AddTemplate(hist_bkg2, "bkg 2", 0., 100.);
	m_2bkg1sgn->AddTemplate(hist_sgn, "sgn", 0., 50.);

	// set priors
	m_2bkg1sgn->SetPriorConstant("bkg 1");
	m_2bkg1sgn->SetPriorConstant("bkg 2");
	m_2bkg1sgn->SetPriorConstant("sgn");

	// ----------------------------------------------------
	// perform analyses
	// ----------------------------------------------------

	// initialize models
	m_1bkg0sgn->Initialize();
	m_1bkg1sgn->Initialize();
	m_2bkg0sgn->Initialize();
	m_2bkg1sgn->Initialize();

	// run MCMC and find global mode
	m_1bkg0sgn->MarginalizeAll();
	m_1bkg0sgn->FindMode();
	m_1bkg1sgn->MarginalizeAll();
	m_1bkg1sgn->FindMode();
	m_2bkg0sgn->MarginalizeAll();
	m_2bkg0sgn->FindMode();
	m_2bkg1sgn->MarginalizeAll();
	m_2bkg1sgn->FindMode();

	// ----------------------------------------------------
	// print
	// ----------------------------------------------------

	std::cout << " Model 1 bkg 0 sgn : " << std::endl;
	std::cout << " Chi2 / ndf = " << m_1bkg0sgn->CalculateChi2( m_1bkg0sgn->GetBestFitParameters() ) << " / ";
	std::cout << m_1bkg0sgn->GetNDF() << " (";
	std::cout << m_1bkg0sgn->CalculateChi2Prob( m_1bkg0sgn->GetBestFitParameters() ) << ")" << std::endl;
	std::cout << " KS probability = " << m_1bkg0sgn->CalculateKSProb() << std::endl;
	std::cout << std::endl;

	std::cout << " Model 1 bkg 1 sgn : " << std::endl;
	std::cout << " Chi2 / ndf = " << m_1bkg1sgn->CalculateChi2( m_1bkg1sgn->GetBestFitParameters() ) << " / ";
	std::cout << m_1bkg1sgn->GetNDF() << " (";
	std::cout << m_1bkg1sgn->CalculateChi2Prob( m_1bkg1sgn->GetBestFitParameters() ) << ")" << std::endl;
	std::cout << " KS probability = " << m_1bkg1sgn->CalculateKSProb() << std::endl;
	std::cout << std::endl;

	std::cout << " Model 2 bkg 0 sgn : " << std::endl;
	std::cout << " Chi2 / ndf = " << m_2bkg0sgn->CalculateChi2( m_2bkg0sgn->GetBestFitParameters() ) << " / ";
	std::cout << m_2bkg0sgn->GetNDF() << " (";
	std::cout << m_2bkg0sgn->CalculateChi2Prob( m_2bkg0sgn->GetBestFitParameters() ) << ")" << std::endl;
	std::cout << " KS probability = " << m_2bkg0sgn->CalculateKSProb() << std::endl;
	std::cout << std::endl;

	std::cout << " Model 2 bkg 1 sgn : " << std::endl;
	std::cout << " Chi2 / ndf = " << m_2bkg1sgn->CalculateChi2( m_2bkg1sgn->GetBestFitParameters() ) << " / ";
	std::cout << m_2bkg1sgn->GetNDF() << " (";
	std::cout << m_2bkg1sgn->CalculateChi2Prob( m_2bkg1sgn->GetBestFitParameters() ) << ")" << std::endl;
	std::cout << " KS probability = " << m_2bkg1sgn->CalculateKSProb() << std::endl;
	std::cout << std::endl;


	// print data
	TCanvas c1("c1");
	c1.cd();
	hist_data.Draw();
	c1.Print("data.ps");

	// print marginalized distribution for signal
	c1.cd();
	m_1bkg1sgn->GetMarginalized("sgn")->Draw(0, -95);
	c1.Print("1bkg1sgn_sgn.ps");

	c1.cd();
	m_2bkg1sgn->GetMarginalized("sgn")->Draw(0, -95);
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

	return 0;
}

