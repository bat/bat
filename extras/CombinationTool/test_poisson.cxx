#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h> 

#include "CombinationXSec.h"
#include <iostream>
int main()
{
	// ----------------------------------------------------------
	// setup BAT infrastructure
	// ----------------------------------------------------------

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// create new CombinationModel object
	// and define the parameter region
	CombinationXSec * model = new CombinationXSec("#sigma [pb]", 0.0, 20.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNChains(10);
	model->MCMCSetNIterationsRun(1000000);

	// ----------------------------------------------------------
	// define cross-section contributions, background sources and
	// systematics here 
	// ----------------------------------------------------------

	//
	// set fitting options
	//
	model->SetFlagSystErrors(true);

	//
	// add channels 
	// 

	// add channel
	model->AddChannel("ee");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("ee",  55); 
	model->SetChannelEfficiency("ee", 0.011);
	model->SetChannelLuminosity("ee", 4270.);
	model->SetChannelBR("ee", 0.10498);

	// add backgrounds for this channel 
	model->AddChannelBackground("ee", "Z->ll",   8.5);
	model->AddChannelBackground("ee", "Diboson", 2.1); 
	model->AddChannelBackground("ee", "fake e",  0.1); 

	/*
	// add channel
	model->AddChannel("emu");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("emu", 204 ); 
	model->SetChannelEfficiency("emu", 0.0427);
	model->SetChannelLuminosity("emu", 4280.);
	model->SetChannelBR("emu", 0.10498);

	// add backgrounds for this channel 
	model->AddChannelBackground("emu", "Z->tautau", 11.9);
	model->AddChannelBackground("emu", "Diboson",    6.5); 
	model->AddChannelBackground("emu", "fake e",     8.1); 
	model->AddChannelBackground("emu", "fake mu",    2.6); 
	*/

	// ----------------------------------------------------------
	// run analysis and plot
	// ----------------------------------------------------------

	// perform analysis
	model->PerformAnalysis();

	// get histogram
	BCH1D* hist_signal = model->GetMarginalized( model->GetParameter(0)->GetName().c_str() );

	// define histogram
	double xmin = model->GetParameter(0)->GetLowerLimit();
	double xmax = model->GetParameter(0)->GetUpperLimit(); 
	TH1D* hist_calc = new TH1D("hist_calc", "", 100, xmin, xmax);
	hist_calc->SetStats(kFALSE); 
	hist_calc->SetLineColor(kRed); 
	hist_calc->SetLineStyle(1);
	hist_calc->SetLineWidth(1); 

	// calculate function
	for (int i = 1; i <= 100; ++i) {
		double x1 = hist_calc->GetBinCenter(i) * 4280. * 0.011 * 0.10498;
		double x2 = hist_calc->GetBinCenter(i) * 4280. * 0.0427 * 0.10498;
		double p = 1.0 
			* TMath::PoissonI(55, x1+8.5+2.1+0.1);
			//			* TMath::PoissonI(204, x2+11.9+6.5+8.1+2.6);
		hist_calc->SetBinContent(i, p); 
	}

	// create new canvas
	TCanvas* c1 = new TCanvas();
	c1->cd();
	TH1D* hist = hist_signal->GetHistogram(); 
	hist->Scale( 1.0/hist->Integral() );
	hist->Draw("");
	hist_calc->Scale( 1.0/hist_calc->Integral() );
	hist_calc->Draw("SAME");
	c1->Print("test_poisson_plots.ps"); 

	std::cout << hist_calc->GetBinCenter(hist_calc->GetMinimumBin()) << std::endl;

	// print results
	model->PrintResults("test_poisson_results.txt");
	
	// clean up memory
	delete model;

	// no error
	return 1;
}
