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
	CombinationXSec * model = new CombinationXSec("#sigma [pb]", 0.0, 25.0);

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
	model->SetChannelObservation("e+jets",  200.0); 
	model->SetChannelEfficiencyPriorGauss("e+jets", 0.01, 0.001);
	model->SetChannelLuminosityPriorGauss("e+jets", 1000.0, 50.0);
	model->SetChannelBR("e+jets", 0.9);

	// add backgrounds for this channel 
	model->AddChannelBackground("e+jets", "W+jets", 0.0, 100.0); 
	model->AddChannelBackground("e+jets", "Z+jets", 0.0,  12.0);
	model->AddChannelBackground("e+jets", "QCD",    0.0,  50.0);
	model->AddChannelBackground("e+jets", "other",  5.0,  15.0);

	// set channel priors
 	model->SetChannelBackgroundPriorGauss("e+jets", "W+jets", 55.0,  4.8, 5.0); 
 	model->SetChannelBackgroundPriorGauss("e+jets", "Z+jets",  5.0,  1.5);
 	model->SetChannelBackgroundPriorGauss("e+jets", "QCD",    30.0,  3.7, 3.5);
 	model->SetChannelBackgroundPriorGauss("e+jets", "other",  10.0,  1.1, 1.0); 



	// add channel
	model->AddChannel("mu+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("mu+jets", 200.0); 
	model->SetChannelEfficiencyPriorGauss("mu+jets", 0.01, 0.001);
	model->SetChannelLuminosityPriorGauss("mu+jets", 1000.0, 50.0);
	model->SetChannelBR("mu+jets", 0.9);

	// add backgrounds
	model->AddChannelBackground("mu+jets", "W+jets", 0.0, 100.0);
	model->AddChannelBackground("mu+jets", "Z+jets", 5.0,  15.0); 
	model->AddChannelBackground("mu+jets", "QCD",    0.0,  50.0); 
	model->AddChannelBackground("mu+jets", "other",  0.0,  20.0); 

	// set channel priors
	model->SetChannelBackgroundPriorGauss("mu+jets", "W+jets", 60.0, 4.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "Z+jets", 11.0, 1.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "QCD",    13.0, 5.4, 5.2); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "other",  10.0, 1.5, 1.4); 

	//
	// add systematics
	//

	// add systematic
	model->AddSystError("JES");
	
	// define uncertainty for each channel and background source
	model->SetSystErrorChannelBackground("JES", "e+jets", "W+jets", 10.0, 10.0);
	model->SetSystErrorChannelBackground("JES", "e+jets", "Z+jets", 0.5, 0.7);
	model->SetSystErrorChannelBackground("JES", "e+jets", "QCD",    5.0, 4.9);
	model->SetSystErrorChannelBackground("JES", "e+jets", "other",  1.0, 1.1);

	model->SetSystErrorChannelBackground("JES", "mu+jets", "W+jets", 7, 6);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "Z+jets", 0.5, 0.7);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "QCD",    3.5, 3.6);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "other",  1.0, 1.1);

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	model->PerformAnalysis();

	// print summary plots
	BCSummaryTool* summary = new BCSummaryTool(model);
	summary->PrintParameterPlot("summary_parameters.ps");
	summary->PrintCorrelationPlot("summary_correlation.ps");
	summary->PrintKnowlegdeUpdatePlot("summary_update.ps"); 

	// print results
	model->PrintAllMarginalized("model_plots.ps");
	model->PrintResults("model_results.txt");
	model->PrintChannelOverview("channels.ps");

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

