// BAT 
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

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
	CombinationModel * model = new CombinationModel("#sigma [pb]", -5.0, 19.0);

	// set mcmc options
	//	model->MCMCSetNLag(10);
	//	model->MCMCSetNChains(5);
	//	model->MCMCSetNIterationsRun(10000000); // high precision
	//	model->MCMCSetNIterationsRun(100000); // low precision
	//	model->SetNbins("#sigma [pb]", 400); // high precision 
	//	model->SetNbins("#sigma [pb]", 100); // low precision

	model->MCMCSetPrecision(BCEngineMCMC::kHigh);
	model->SetNbins("#sigma [pb]", 400); // high precision 


	// ----------------------------------------------------------
	// define quantites here
	// ----------------------------------------------------------

	//
	// set fitting options
	//
	model->SetFlagSystErrors(true);

	int runtype = 12;

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

	if (runtype == 10) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 5.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 7.0, 3.0, 3.0);
	}

	if (runtype == 11) {
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
		model->SetSystErrorChannelSignal("syst1", "channel1", 1.0, 1.0, 0.00);
		model->SetSystErrorChannelSignal("syst1", "channel2", 0.0, 0.0, 0.00);

		model->AddSystError("syst2");
		model->SetNbins("syst2", 400); // high precision 
		model->SetSystErrorChannelSignal("syst2", "channel1", 1.0, 1.0, 0.00);
		model->SetSystErrorChannelSignal("syst2", "channel2", 1.0, 1.0, 0.00);
	}

	if (runtype == 12) {
		// add channel
		model->AddChannel("channel1");
		model->AddChannel("channel2");
		
		// add channels
		// parameters: channel name, mean value, -sigma, +sigma
		model->SetChannelSignalPriorGauss("channel1", 6.0, 1.0, 1.0);
		model->SetChannelSignalPriorGauss("channel2", 6.0, 1.0, 1.0);

		// add systematics
		// parameters: uncertainty, channel, -sigma, +sigma, mean
		model->AddSystError("syst1");
		model->SetNbins("syst1", 400); // high precision 
		model->SetSystErrorChannelSignal("syst1", "channel1", 2.0, 2.0, -2.00);
		model->SetSystErrorChannelSignal("syst1", "channel2", 2.0, 2.0, -2.00);

		model->AddSystError("syst2");
		model->SetNbins("syst2", 400); // high precision 
		model->SetSystErrorChannelSignal("syst2", "channel1", 1.0, 1.0,  1.00);
		model->SetSystErrorChannelSignal("syst2", "channel2", 1.0, 1.0,  1.00);
	}

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform full analysis
	model->PerformFullAnalysis();

	// print results
	//	model->PrintChannelOverview(Form("model_channels_%i.ps", runtype), Form("model_systematics_%i.ps", runtype));
	model->PrintChannelSummary(Form("model_summary_%i.txt", runtype));

	// ----------------------------------------------------------
	// Perform fit again and summarize
	// ----------------------------------------------------------

	// create summary tool
	BCSummaryTool* st = new BCSummaryTool(model); 

	// perform analysis
	model->MarginalizeAll(); 
	model->FindMode(model->GetBestFitParameters());

	// print default plots
	model->PrintAllMarginalized(Form("model_plots_%i.ps", runtype));
	model->PrintResults(Form("model_results_%i.txt", runtype));

	st->PrintParameterPlot(Form("model_parameters_%i.ps", runtype)); 
	st->PrintCorrelationPlot(Form("model_correlation_%i.ps", runtype)); 
	//	st->PrintKnowledgeUpdatePlots(Form("model_update_%i.ps", runtype)); 
	st->PrintCorrelationMatrix(Form("model_matrix_%i.ps", runtype));

	// ----------------------------------------------------------
	// clean-up and return
	// ----------------------------------------------------------

	// close log file
	BCLog::CloseLog();

	// clean up memory
	delete model;

	// delete summary tool
	//	delete st; 

return 0;

}

