// BAT 
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h> 

#include "CombinationXSec.h"

// ------------------------------------------------------------
int main()
{

	// ----------------------------------------------------------
	// setup BAT infrastructure
	// ----------------------------------------------------------

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::summary);

	// create new CombinationModel object
	// and define the parameter region
	CombinationXSec * model = new CombinationXSec("#sigma [pb]", 0.0, 20.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	//	model->MCMCSetNIterationsRun(1000);
	//	model->MCMCSetNIterationsRun(100000);

	// ----------------------------------------------------------
	// define cross-section contributions, background sources and
	// systematics here 
	// ----------------------------------------------------------

	// add channel
	model->AddChannel("e+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("e+jets",  200.0); 
	model->SetChannelEfficiencyPriorGauss("e+jets", 0.01, 0.001);
	model->SetChannelLuminosity("e+jets",  1000.0);
	model->SetChannelBR("e+jets", 0.9);

	// add backgrounds for this channel 
	model->AddChannelBackground("e+jets", "W+jets", 0.0, 100.0); 
	model->AddChannelBackground("e+jets", "Z+jets", 0.0,  12.0);
	model->AddChannelBackground("e+jets", "QCD",    0.0,  50.0);
	model->AddChannelBackground("e+jets", "other",  5.0,  15.0);

	// set channel priors
 	model->SetChannelBackgroundPriorGauss("e+jets", "W+jets", 55.0,  5.7, 6.0); 
 	model->SetChannelBackgroundPriorGauss("e+jets", "Z+jets",  5.0,  1.5);
 	model->SetChannelBackgroundPriorGauss("e+jets", "QCD",    30.0,  3.7, 3.5);
 	model->SetChannelBackgroundPriorGauss("e+jets", "other",  10.0,  1.1, 1.0); 

	// add channel
	model->AddChannel("mu+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("mu+jets", 200.0); 
	model->SetChannelEfficiencyPriorGauss("mu+jets", 0.01, 0.001);
	model->SetChannelLuminosity("mu+jets", 1000.0);
	model->SetChannelBR("mu+jets", 0.9);

	// add backgrounds
	model->AddChannelBackground("mu+jets", "W+jets", 0.0, 100.0);
	model->AddChannelBackground("mu+jets", "Z+jets", 5.0,  15.0); 
	model->AddChannelBackground("mu+jets", "QCD",    0.0,  20.0); 
	model->AddChannelBackground("mu+jets", "other",  0.0,  20.0); 

	// set channel priors
	model->SetChannelBackgroundPriorGauss("mu+jets", "W+jets", 50.0, 5.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "Z+jets", 10.0, 1.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "QCD",    10.0, 1.9, 2.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "other",  10.0, 1.5, 1.4); 

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// run MCMC and Minuit
	model->MarginalizeAll();
	model->FindMode( model->GetBestFitParameters() );

	// print results
	model->PrintAllMarginalized("model_plots.ps");
	model->PrintResults("model_results.txt");
	model->PrintChannelOverview("channels.ps");

	// print summary plots
	BCSummaryTool* summary = new BCSummaryTool(model);
	summary->PrintParameterPlot("summary_parameters.ps");
	summary->PrintCorrelationPlot("summary_correlation.ps");
	summary->PrintKnowlegdeUpdatePlot("summary_update.ps"); 

	// ----------------------------------------------------------
	// clean-up and return
	// ----------------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// clean up memory
	delete model;
	delete summary;

	return 0;

}

