// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project CombinationTool
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h> 

#include "CombinationXSec.h"

int main()
{

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new CombinationModel object
	CombinationXSec * m = new CombinationXSec("sigma ttbar [pb]", 0.0, 20.0);

	// set mcmc options
	m->MCMCSetNLag(10);
	m->MCMCSetNIterationsRun(100000);

	// perform your analysis here

	// add channel
	m->AddChannel("e+jets");
	
	// set channel observations, efficiency and luminosity 
	m->SetChannelObservation("e+jets",  100.0); 
	m->SetChannelEfficiency("e+jets",     0.01);
	m->SetChannelLuminosity("e+jets",  1000.0);

	// add backgrounds for this channel 
	m->AddChannelBackground("e+jets", "W+jets", 0.0, 100.0); 
	m->AddChannelBackground("e+jets", "Z+jets", 0.0,  12.0);
	m->AddChannelBackground("e+jets", "QCD",    0.0,  50.0);
	m->AddChannelBackground("e+jets", "other",  0.0,  20.0);

	// set channel priors
	m->SetChannelSignalPriorGauss("e+jets",              100.0, 20.0);
 	m->SetChannelBackgroundPriorGauss("e+jets", "W+jets", 50.0,  5.0, 7.0); 
 	m->SetChannelBackgroundPriorGauss("e+jets", "Z+jets",  5.0,  1.5);
 	m->SetChannelBackgroundPriorGauss("e+jets", "QCD",    30.0,  4.0, 2.0);
 	m->SetChannelBackgroundPriorGauss("e+jets", "other",  10.0,  2.0, 3.0); 

	// add channel
	m->AddChannel("mu+jets");

	// set channel observations, efficiency and luminosity 
	m->SetChannelObservation("mu+jets", 100.0); 
	m->SetChannelEfficiency("mu+jets",    0.01);
	m->SetChannelLuminosity("mu+jets", 1000.0);

	// add backgrounds
	m->AddChannelBackground("mu+jets", "W+jets", 0.0, 100.0);
	m->AddChannelBackground("mu+jets", "Z+jets", 0.0,  12.0); 
	m->AddChannelBackground("mu+jets", "QCD",    0.0,  20.0); 
	m->AddChannelBackground("mu+jets", "other",  0.0,  20.0); 

	// set channel priors
	m->SetChannelSignalPriorGauss("mu+jets",              60.0, 5.0);
	m->SetChannelBackgroundPriorGauss("mu+jets", "W+jets", 60.0, 5.0); 
	m->SetChannelBackgroundPriorGauss("mu+jets", "Z+jets",  7.0, 1.0); 
	m->SetChannelBackgroundPriorGauss("mu+jets", "QCD",    10.0, 1.0, 2.0); 
	m->SetChannelBackgroundPriorGauss("mu+jets", "other",  10.0, 2.0, 1.0); 

	m -> MarginalizeAll();
	m -> FindMode( m -> GetBestFitParameters() );

	m -> PrintAllMarginalized("model_plots.ps");
	m -> PrintResults("model_results.txt");

	BCSummaryTool* summary = new BCSummaryTool(m);
	summary->PrintParameterPlot("summary_parameters.ps");
	summary->PrintCorrelationPlot("summary_correlation.ps");
	summary->PrintKnowlegdeUpdatePlot("summary_update.ps"); 

	// close log file
	BCLog::CloseLog();

	delete m;
	delete summary;

	return 0;

}

