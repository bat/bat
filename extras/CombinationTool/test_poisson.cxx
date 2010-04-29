#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h> 

#include "CombinationXSec.h"

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
	model->AddChannel("e+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("e+jets",  3); 
	model->SetChannelEfficiencyPriorGauss("e+jets", 1., 0.);
	model->SetChannelLuminosityPriorGauss("e+jets", 1., 0.);
	model->SetChannelBR("e+jets", 1.0);

	// add backgrounds for this channel 
	model->AddChannelBackground("e+jets", "W+jets", 0.4, 0.4);
	model->AddChannelBackground("e+jets", "Z+jets", 0.2, 0.2); 
	model->AddChannelBackground("e+jets", "QCD",    0.5, 0.5); 
	model->AddChannelBackground("e+jets", "other",  0.9, 0.9); 


	// add channel
	model->AddChannel("mu+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("mu+jets",  9); 
	model->SetChannelEfficiencyPriorGauss("mu+jets", 1., 0.);
	model->SetChannelLuminosityPriorGauss("mu+jets", 1., 0.);
	model->SetChannelBR("mu+jets", 1.0);

	// add backgrounds for this channel 
	model->AddChannelBackground("mu+jets", "W+jets", 1.4, 1.4);
	model->AddChannelBackground("mu+jets", "Z+jets", 0.2, 0.2); 
	model->AddChannelBackground("mu+jets", "QCD",    1.5, 1.5); 
	model->AddChannelBackground("mu+jets", "other",  0.9, 0.9); 

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
		double x = hist_calc->GetBinCenter(i); 
		double p = TMath::PoissonI(3, x+2)
			* TMath::PoissonI(9, x+4);
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

	// print results
	model->PrintResults("test_poisson_results.txt");
	
	// clean up memory
	delete model;

	// no error
	return 1;
}
