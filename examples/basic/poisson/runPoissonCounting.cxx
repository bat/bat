#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "PoissonModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new PoissonModel object
  PoissonModel* m = new PoissonModel();

  // set precision
  m->MCMCSetPrecision(BCEngineMCMC::kHigh);

  // set number of observed events
  m->SetNObs(7);

  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters
  m->MarginalizeAll();

  // find mode starting from the best fit parameters
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized("PoissonModel_plots.eps");

	//	m->PrintSummary();

  // print results of the analysis into a text file
  m->PrintResults("PoissonModel_results.txt");

  // close log file
  BCLog::CloseLog();

  // free memory
  delete m;

  // no error
  return 0;

}

