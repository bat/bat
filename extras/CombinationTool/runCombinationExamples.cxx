// BAT 
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <TROOT.h>

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
	CombinationModel * model = new CombinationModel("#sigma [pb]", 0.0, 12.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNChains(5);
	model->MCMCSetNIterationsRun(10000000); // high precision
	//	model->MCMCSetNIterationsRun(1000000); // low precision
	model->SetNbins("#sigma [pb]", 400); // high precision 
	//	model->SetNbins("#sigma [pb]", 100); // low precision

	// ----------------------------------------------------------
	// define quantites here
	// ----------------------------------------------------------

	//
	// set fitting options
	//
	model->SetFlagSystErrors(true);

	int runtype = 9;

	if (runtype == 0) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);
	}

	if (runtype == 1) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, sqrt(2.), sqrt(2.));
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);
	}

	if (runtype == 2) {
		// add channel
		model->AddChannel("channel1");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.00, 1.00, 0.00);
	}

	if (runtype == 3) {
		// add channel
		model->AddChannel("channel1");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 0.00001, 0.00001, 1.00);
	}

	if (runtype == 4) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 0.00001, 0.00001, 1.00);
	}

	if (runtype == 5) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.0, 1.0, 0.00);
	}

	if (runtype == 6) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.0, 1.0, 0.00);
		model->SetSystErrorChannelSignal("syst1", "channel2", 1.0, 1.0, 0.00);
	}

	if (runtype == 7) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.0, 1.0, 1.00);

		model->AddSystError("syst2");
		model->SetSystErrorChannelSignal("syst2", "channel2", 1.0, 1.0, -1.00);
	}

	if (runtype == 8) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.0, 1.0, 1.00);
		model->SetSystErrorChannelSignal("syst1", "channel2", 1.0, 1.0, -1.00);
	}

	if (runtype == 9) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetNbins("syst1", 400); // high precision 
		model->SetSystErrorChannelSignal("syst1", "channel1", 0.7, 0.7,  0.00);
		model->SetSystErrorChannelSignal("syst1", "channel2", 1.0, 1.0, -1.00);
	}


	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	model->PerformFullAnalysis();
	model->PerformAnalysis();

	// print results
	model->PrintAllMarginalized(Form("model_plots_%i.ps", runtype));
	model->PrintResults(Form("model_results_%i.txt", runtype));
	model->PrintChannelOverview(Form("model_channels_%i.ps", runtype));
	model->PrintChannelSummary(Form("model_summary_%i.txt", runtype));

	// ----------------------------------------------------------
	// clean-up and return
	// ----------------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// clean up memory
	delete model;

	return 0;

}

