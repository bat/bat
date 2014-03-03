#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "BinomialModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new BinomialModel object
  BinomialModel * m = new BinomialModel();

  // set number of events
  m->SetNTotal(20);
  m->SetNSelected(3);

  // run the MCMC and calculate the marginalized distributions
  m->MarginalizeAll();

  // find mode starting from the best fit parameters
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a pdf file
  m->PrintAllMarginalized("BinomialModel_plots.pdf");

  // print results of the analysis into a text file
  m->PrintResults("BinomialModel_results.txt");

  // close log file
  BCLog::CloseLog();

  // free memory
  delete m;

  // no errors
  return 0;

}

