// BAT 
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "CombinationModel.h"

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
	CombinationModel * model = new CombinationModel("#sigma [pb]", 5.0, 10.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNIterationsRun(100000);
	model->MCMCSetNChains(5);
	model->SetNbins("#sigma [pb]", 100);

	// ----------------------------------------------------------
	// define quantites here
	// ----------------------------------------------------------

	//
	// set fitting options
	//
	model->SetFlagSystErrors(false);

	//
	// add channels 
	// 

	// add channel
	model->AddChannel("e+jets");
	model->AddChannel("mu+jets");

	model->SetChannelSignalPriorGauss("e+jets",  6.77, 0.42, 0.42);
	model->SetChannelSignalPriorGauss("mu+jets", 8.03, 0.54, 0.55);

	//
	// add systematics
	// 

	// add systematic: higher order effects, 100% correlation among signal channels
	//	model->AddSystError("JES");
	//	model->SetSystErrorChannelSignal("JES", "e+jets", 0.5, 0.4);
	//	model->SetSystErrorChannelSignal("JES", "mu+jets", 0.5, 0.5);

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	model->PerformFullAnalysis();
	model->PerformAnalysis();

	// print results
	model->PrintAllMarginalized("model_plots.ps");

	model->PrintResults("model_results.txt");
	model->PrintChannelOverview("channels.ps");
	model->PrintChannelSummary("summary.txt");

	// ----------------------------------------------------------
	// clean-up and return
	// ----------------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// clean up memory
	delete model;

	return 0;

}

