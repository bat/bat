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

  // marginalize
  m->MarginalizeAll();

  // find mode
  m->FindMode( m->GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized("CombinationModel_plots.pdf");

  // print all summary plots
  summary->PrintParameterPlot("CombinationModel_parameters.pdf");
  summary->PrintCorrelationPlot("CombinationModel_correlation.pdf");
  summary->PrintKnowledgeUpdatePlots("CombinationModel_update.pdf");

  // print results of the analysis into a text file
  m->PrintResults("CombinationModel_results.txt");

  // close log file
  BCLog::CloseLog();

  delete m;
  delete summary;

  return 0;

}

