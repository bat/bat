#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "GaussModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new GaussModel object
  GaussModel * m = new GaussModel();

  // set marginalization method
  m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

  // set MCMC precision
  m->MCMCSetPrecision(BCEngineMCMC::kMedium);

  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters
  m->MarginalizeAll();

  // if MCMC was run before (MarginalizeAll()) it is
  // possible to use the mode found by MCMC as
  // starting point of Minuit minimization
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized("GaussModel_plots.pdf");

  // print results of the analysis into a text file
  m->PrintResults("GaussModel_results.txt");

  // close log file
  BCLog::CloseLog();

  delete m;

  return 0;

}

