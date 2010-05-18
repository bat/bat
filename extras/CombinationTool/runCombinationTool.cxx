// BAT 
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

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
	CombinationXSec * model = new CombinationXSec("#sigma [pb]", 0.0, 20.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNIterationsRun(100000);
	model->MCMCSetNChains(5);
	model->SetNbins("#sigma [pb]", 200);

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
	model->SetChannelEfficiency("ee", 0.0110102);
	model->SetChannelLuminosity("ee", 4280.);
	model->SetChannelBR("ee", 0.1049760);

	// add backgrounds for this channel 
	model->AddChannelBackground("ee", "Z->ll",   8.46994);
	model->AddChannelBackground("ee", "Diboson", 2.12045); 
	
	// add channel
	model->AddChannel("emu");
	
	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("emu", 204 );                             
	model->SetChannelEfficiency("emu", 0.042805);                            
	model->SetChannelLuminosity("emu", 4280.);
	model->SetChannelBR("emu", 0.1049760);
	
	// add backgrounds for this channel
	model->AddChannelBackground("emu", "Z->ll", 11.9332);
	model->AddChannelBackground("emu", "Diboson", 6.54512);
	model->AddChannelBackground("emu", "fake e", 8.09789);
	model->AddChannelBackground("emu", "fake mu", 2.64443);

	// add channel
	model->AddChannel("mumu");
	
	// set channel observations, efficiency, luminosity and branching ratio
	model->SetChannelObservation("mumu", 72 );                             
	model->SetChannelEfficiency("mumu", 0.0134512);                            
	model->SetChannelLuminosity("mumu", 4280.);
	model->SetChannelBR("mumu", 0.1049760);
	
	// add backgrounds for this channel
	model->AddChannelBackground("mumu", "Z->ll", 21.7316);
	model->AddChannelBackground("mumu", "Diboson", 3.28345);
	model->AddChannelBackground("mumu", "fake mu", 3.21932);

	//
	// add systematics
	// note: the uncertainties are given in % of the efficiency for the signal, 
  //       and in % of the background contribution for the background.
	//

	// add systematic: higher order effects, 100% correlation among signal channels
	model->AddSystError("HO");
	model->SetSystErrorChannelSignal("HO", "ee", 0.02, 0.02);
	model->SetSystErrorChannelSignal("HO", "emu", 0.02, 0.02);
	model->SetSystErrorChannelSignal("HO", "mumu", 0.02, 0.02);

	// add systematic: branching ratio, 100% correlation among signal channels
	model->AddSystError("BR");
	model->SetSystErrorChannelSignal("BR", "ee", 0.017, 0.017);
	model->SetSystErrorChannelSignal("BR", "emu", 0.017, 0.017);
	model->SetSystErrorChannelSignal("BR", "mumu", 0.017, 0.017);

	// add systematic: jet energy scale, 100% correlation among all channels (signal and background)
	model->AddSystError("JES");
	model->SetSystErrorChannelSignal("JES", "ee",                0.0192, 0.0181);
	model->SetSystErrorChannelBackground("JES", "ee", "Z->ll",   0.1626, 0.1621);	
	model->SetSystErrorChannelBackground("JES", "ee", "Diboson", 0.0820, 0.0869);
	model->SetSystErrorChannelSignal("JES", "emu",                0.0169, 0.0152);
	model->SetSystErrorChannelBackground("JES", "emu", "Z->ll",   0.0323, 0.0431);	
	model->SetSystErrorChannelBackground("JES", "emu", "Diboson", 0.0504, 0.0654);
	model->SetSystErrorChannelSignal("JES", "mumu",                0.0173, 0.0196 );
	model->SetSystErrorChannelBackground("JES", "mumu", "Z->ll",   0.0128, 0.0033 );	
	model->SetSystErrorChannelBackground("JES", "mumu", "Diboson", 0.1341, 0.0449);

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	model->PerformFullAnalysis();
	//	model->PerformAnalysis();

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

