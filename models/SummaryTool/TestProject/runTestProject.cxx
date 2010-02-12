#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "TestModel.h"
#include "SummaryTool.h"

int main()
{
	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create test model 
	TestModel * model = new TestModel(); 

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNIterationsRun(10000000);
	model->SetNbins((model->GetParameter(0)->GetName()).c_str(), 100);
	model->SetNbins((model->GetParameter(1)->GetName()).c_str(), 100);

	// run analysis 
	model->MarginalizeAll(); 
	model->FindMode(model->GetBestFitParameters()); 
	model->PrintAllMarginalized("model_plots.ps"); 
	model->PrintResults("model_results.txt"); 

	// create new SummaryTool object
	SummaryTool * summarytool = new SummaryTool();

	// set model
	summarytool->SetModel(model); 

	// print summaries 
 	summarytool->PrintParameterPlot("summary_parameters.eps"); 
 	summarytool->PrintCorrelationPlot("summary_correlation.eps"); 
	summarytool->PrintKnowlegdeUpdatePlot("summary_update.ps");

	// free memory
	delete summarytool;
	delete model; 

	return 0;

}

