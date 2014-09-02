#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "GaussModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

  // create new GaussModel object
  GaussModel * m = new GaussModel("gausMod");

  // set marginalization method
  m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

  // set MCMC precision
  m->MCMCSetPrecision(BCEngineMCMC::kMedium);

  // switch writing of Markov Chains on
  m->WriteMarkovChain("output.root","RECREATE",true);

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

  // write marginalized distributions to output file
  m->WriteMarginalizedDistributions("output.root","UPDATE");

  // close log file
  BCLog::CloseLog();

  // free memory
  delete m;

  return 0;

}

