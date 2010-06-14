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
   CombinationModel * model = new CombinationModel("#sigma [pb]", 3.0,
15.0);

   // set mcmc options
   model->MCMCSetNLag(10);
   model->MCMCSetNChains(5);
   //    model->MCMCSetNIterationsRun(10000000); // high precision
   model->MCMCSetNIterationsRun(100000); // low precision
   //    model->SetNbins("#sigma [pb]", 400); // high precision
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
//    model->SetChannelSignalPriorGauss("e+jets",  6.49, 0.40, 0.41);
//    model->SetChannelSignalPriorGauss("mu+jets", 7.94, 0.53, 0.53);
// stat error scaled to final cross section
   model->SetChannelSignalPriorGauss("e+jets",  6.53, 0.40, 0.41);
   model->SetChannelSignalPriorGauss("mu+jets", 8.37, 0.56, 0.56);

   //
   // add systematics
   //

   // parameters: uncertinty, channel, -sigma, +sigma, mean
//    model->SetSystErrorChannelSignal("ID_mu+jets", "mu+jets",             0.18, 0.21,  0.00 );

////////////////// all shifts set to 0        ////////////////
   // systematic: event preselection, correlated among signal channels
   model->AddSystError("preselection");
   model->SetSystErrorChannelSignal("preselection", "e+jets",     0.12, 0.12,  0.00);
   model->SetSystErrorChannelSignal("preselection", "mu+jets",    0.15, 0.15,  0.00);

   // ----------------------------------------------------------
   // run analysis and plotting
   // ----------------------------------------------------------

   // perform analysis
   model->PerformFullAnalysis();
   //   model->PerformAnalysis();

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

