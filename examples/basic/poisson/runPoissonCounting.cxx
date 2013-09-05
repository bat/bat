#include "PoissonModel.h"

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>

int main()
{
  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new PoissonModel object
  PoissonModel m;

  // set precision
  m.MCMCSetPrecision(BCEngineMCMC::kMedium);

  // set number of observed events
  m.SetNObs(7);

  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters
  m.MarginalizeAll();

  // find mode starting from the best fit parameters
  m.FindMode( m.GetBestFitParameters() );

  std::cout << " N "  << m.GetBestFitParameterErrors().size() << std::endl;

  // draw all marginalized distributions into a PDF file
  m.PrintAllMarginalized("PoissonModel_plots.pdf");

  // print summary on standard output
  m.PrintSummary();

  // print results of the analysis into a text file
  m.PrintResults("PoissonModel_results.txt");

  // close log file
  BCLog::CloseLog();

  // no error
  return 0;
}
