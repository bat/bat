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
	CombinationModel * model = new CombinationModel("#sigma [pb]", 3.0, 15.0);

	// set mcmc options
	model->MCMCSetNLag(10);
	model->MCMCSetNChains(5);
	//	model->MCMCSetNIterationsRun(10000000); // high precision
	model->MCMCSetNIterationsRun(100000); // low precision
	//	model->SetNbins("#sigma [pb]", 400); // high precision 
	model->SetNbins("#sigma [pb]", 100); // low precision

	// ----------------------------------------------------------
	// define quantites here
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
	model->AddChannel("mu+jets");

	// parameters: channel name, mean value, -sigma, +sigma
	model->SetChannelSignalPriorGauss("e+jets",  6.49, 0.40, 0.41);
	model->SetChannelSignalPriorGauss("mu+jets", 7.94, 0.53, 0.53);

	//
	// add systematics
	// 

	// systematic: , not correlated among signal channels
	//	model->AddSystError("ID_e+jets");
	//	model->AddSystError("ID_mu+jets");
	// parameters: uncertinty, channel, -sigma, +sigma, mean
	//	model->SetSystErrorChannelSignal("ID_e+jets", "e+jets",               0.25, 0.28,  0.10 );

	// systematic: , not correlated among signal channels
	//	model->AddSystError("ID_mu+jets");
	// parameters: uncertinty, channel, -sigma, +sigma, mean
	//	model->SetSystErrorChannelSignal("ID_mu+jets", "mu+jets",             0.18, 0.21,  0.00 );

	/*
	// systematic: event preselection, not correlated among signal channels
	model->AddSystError("preselection_e+jets");
	model->SetSystErrorChannelSignal("preselection_e+jets", "e+jets",     0.11, 0.12,  0.02);
	
	// systematic: , not correlated among signal channels
	model->AddSystError("ID_e+jets");
	model->SetSystErrorChannelSignal("ID_e+jets", "e+jets",               0.25, 0.28,  0.10 );

	// systematic: , not correlated among signal channels
	model->AddSystError("Lumi_e+jets");
	model->SetSystErrorChannelSignal("Lumi_e+jets", "e+jets",             0.04, 0.00,  0.00 );

	// systematic: , not correlated among signal channels
	model->AddSystError("Z_pT_e+jets");
	model->SetSystErrorChannelSignal("Z_pT_e+jets", "e+jets",             0.00, 0.01,  0.01);

	// systematic: , not correlated among signal channels
	model->AddSystError("signal_modelling_e+jets");
	model->SetSystErrorChannelSignal("signal_modelling_e+jets", "e+jets", 0.16, 0.18, -0.14);

	// systematic: , not correlated among signal channels
	model->AddSystError("color_e+jets");
	model->SetSystErrorChannelSignal("color_e+jets", "e+jets",            0.04, 0.04,  0.00 );

	// systematic: , not correlated among signal channels
	model->AddSystError("ISR/FSR_e+jets");
	model->SetSystErrorChannelSignal("ISR/FSR_e+jets", "e+jets",          0.17, 0.19, -0.12 );

	// systematic: , not correlated among signal channels
	model->AddSystError("EM_trigger_e+jets");
	model->SetSystErrorChannelSignal("EM_trigger_e+jets", "e+jets",       0.06, 0.10,  0.01 );

	// systematic: , not correlated among signal channels
	model->AddSystError("JES_e+jets");
	model->SetSystErrorChannelSignal("JES_e+jets", "e+jets",              0.00, 0.00,  0.27 );

	// systematic: , not correlated among signal channels
	model->AddSystError("vertex_e+jets");
	model->SetSystErrorChannelSignal("vertex_e+jets", "e+jets",           0.04, 0.07, -0.22 );

	// systematic: , not correlated among signal channels
	model->AddSystError("b-JES_e+jets");
	model->SetSystErrorChannelSignal("b-JES_e+jets", "e+jets",            0.03, 0.02,  0.00);

	// systematic: , not correlated among signal channels
	model->AddSystError("JER_e+jets");
	model->SetSystErrorChannelSignal("JER_e+jets", "e+jets",              0.71, 0.01,  0.00);

	// systematic: , not correlated among signal channels
	model->AddSystError("jet_reco_e+jets");
	model->SetSystErrorChannelSignal("jet_reco_e+jets", "e+jets",         0.10, 0.11, -0.10 );

	// systematic: , not correlated among signal channels
	model->AddSystError("b-fragmenation_e+jets");
	model->SetSystErrorChannelSignal("b-fragementation_e+jets", "e+jets", 0.07, 0.07, -0.02 );

	// systematic: , not correlated among signal channels
	model->AddSystError("MM_e+jets");
	model->SetSystErrorChannelSignal("MM_e+jets", "e+jets",               0.01, 0.03, -0.11 );

	// systematic: , not correlated among signal channels
	model->AddSystError("MC_bkg_e+jets");
	model->SetSystErrorChannelSignal("MC_bkg_e+jets", "e+jets",           0.01, 0.00,  0.01 );

	// systematic: , not correlated among signal channels
	model->AddSystError("MC_br_e+jets");
	model->SetSystErrorChannelSignal("MC_br_e+jets", "e+jets",            0.05, 0.06,  0.00 );

	// systematic: , not correlated among signal channels
	model->AddSystError("MC_bkg_scale_e+jets");
	model->SetSystErrorChannelSignal("MC_bkg_scale_e+jets", "e+jets",     0.02, 0.01,  0.05 );

	// systematic: , not correlated among signal channels
	model->AddSystError("MC_stat_e+jets");
	model->SetSystErrorChannelSignal("MC_stat_e+jets", "e+jets",          0.03, 0.04,  0.00 );

	// systematic: , not correlated among signal channels
	model->AddSystError("W_e+jets");
	model->SetSystErrorChannelSignal("W_e+jets", "e+jets",                0.03, 0.03,  0.06 );

	// systematic: , not correlated among signal channels
	model->AddSystError("PDF_e+jets");
	model->SetSystErrorChannelSignal("PDF_e+jets", "e+jets",              0.09, 0.10,  0.00 );

	// systematic: , not correlated among signal channels
	model->AddSystError("Lumi_e+jets");
	model->SetSystErrorChannelSignal("Lumi_e+jets", "e+jets",             0.40, 0.47,  0.25 );
	*/

	// ----------------------------------------------------------
	// run analysis and plotting
	// ----------------------------------------------------------

	// perform analysis
	model->PerformFullAnalysis();
	model->PerformAnalysis();

	// print results
	model->PrintAllMarginalized("model_plots.ps");

	model->PrintResults("model_results.txt");
	model->PrintChannelOverview("channels.ps", "systematics.ps");
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

