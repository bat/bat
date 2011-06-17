#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "CombinationModel.h"

int main()
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new CombinationModel object
  CombinationModel * m = new CombinationModel();

  // create a new summary tool object
  BCSummaryTool * summary = new BCSummaryTool(m);

  // set precision
  m->MCMCSetPrecision(BCEngineMCMC::kVeryHigh);

  // marginalize
  m->MarginalizeAll();

  // find mode
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized("CombinationModel_plots.ps");

  // print all summary plots
  summary->PrintParameterPlot("CombinationModel_parameters.ps");
  summary->PrintCorrelationPlot("CombinationModel_correlation.ps");
  summary->PrintKnowledgeUpdatePlots("CombinationModel_update.ps");

  // print results of the analysis into a text file
  m->PrintResults("CombinationModel_results.txt");

  // close log file
  BCLog::CloseLog();

  delete m;
  delete summary;

  return 0;

}

