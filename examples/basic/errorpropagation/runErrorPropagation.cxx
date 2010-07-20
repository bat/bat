#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "RatioModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new RatioModel object
  RatioModel * m = new RatioModel();

  // set precision
  m->MCMCSetPrecision(BCEngineMCMC::kHigh);

  // run the MCMC and marginalize w.r.t. to all parameters
  m->MarginalizeAll();

  // find mode using the best fit parameters as start values
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized("RatioModel_plots.ps");

  // print results of the analysis into a text file
  m->PrintResults("RatioModel_results.txt");

  // print ratio
  m->PrintHistogram();

  // close log file
  BCLog::CloseLog();

  // free memory
  delete m;

  // no error
  return 0;

}

