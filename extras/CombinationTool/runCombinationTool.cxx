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
	BCLog::SetLogLevel(BCLog::detail);

	// create new CombinationModel object
	// and define the parameter region
	CombinationXSec * model = new CombinationXSec("#sigma [pb]", 0.0, 25.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNIterationsRun(1000000);
	model->MCMCSetNChains(5);
	model->SetNbins("#sigma [pb]", 200);

	// ----------------------------------------------------------
	// define cross-section contributions, background sources and
	// systematics here 
	// ----------------------------------------------------------

	//
	// set fitting options
	//
	model->SetFlagSystErrors(false);

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
	model->AddChannelBackground("ee", "fakes",   0.1); 

	// set channel priors
	//	model->SetChannelBackgroundPriorGauss("ee", "Z->ll", 8.5, 3.4); 
	//	model->SetChannelBackgroundPriorGauss("ee", "Diboson", 2.1, 0.8);
	//	model->SetChannelBackgroundPriorGauss("ee", "fakes", 0.1, 0.1, 0.2);

	/*
	// add channel
	model->AddChannel("mu+jets");

	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("mu+jets", 200.0); 
	model->SetChannelEfficiency("mu+jets", 0.01);
	model->SetChannelLuminosity("mu+jets", 1000.);
	model->SetChannelBR("mu+jets", 1.0);

	// add backgrounds
	model->AddChannelBackground("mu+jets", "W+jets", 0.0, 100.0);
	model->AddChannelBackground("mu+jets", "Z+jets", 5.0,  15.0); 
	model->AddChannelBackground("mu+jets", "QCD",    0.0,  50.0); 
	model->AddChannelBackground("mu+jets", "other",  0.0,  20.0); 

	// set channel priors
	model->SetChannelBackgroundPriorGauss("mu+jets", "W+jets", 60.0, 2.0); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "Z+jets", 10.0, 0.5); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "QCD",    15.0, 2.6, 2.5); 
	model->SetChannelBackgroundPriorGauss("mu+jets", "other",  15.0, 1.5, 1.4); 
	*/

	//
	// add systematics
	//

	/*
	// add systematic
	model->AddSystError("JES");

	// define uncertainty for each channel and background source
	model->SetSystErrorChannelSignal("JES", "e+jets", 0.5, 0.5);
	model->SetSystErrorChannelBackground("JES", "e+jets", "W+jets", 6.0, 6.0);	
	model->SetSystErrorChannelBackground("JES", "e+jets", "Z+jets", 0.5, 0.7);
	model->SetSystErrorChannelBackground("JES", "e+jets", "QCD",    5.0, 4.9);
	model->SetSystErrorChannelBackground("JES", "e+jets", "other",  1.0, 1.1);

	model->SetSystErrorChannelSignal("JES", "mu+jets", 0.5, 0.5);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "W+jets", 7, 6);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "Z+jets", 0.5, 0.7);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "QCD",    3.5, 3.6);
	model->SetSystErrorChannelBackground("JES", "mu+jets", "other",  1.0, 1.1);
	*/

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	//	model->PerformFullAnalysis();
	model->PerformAnalysis();

	// print summary plots
	BCSummaryTool* summary = new BCSummaryTool(model);
	summary->PrintParameterPlot("summary_parameters.ps");
	summary->PrintCorrelationPlot("summary_correlation.ps");
	//	summary->PrintKnowlegdeUpdatePlot("summary_update.ps"); 

	// print results
	model->PrintAllMarginalized("model_plots.ps");

	model->PrintResults("model_results.txt");
	//	model->PrintChannelOverview("channels.ps");
	//	model->PrintChannelSummary("summary.txt");

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

