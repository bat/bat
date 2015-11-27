#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "CombinationModel.h"

int main()
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new CombinationModel object
    CombinationModel m("combMod",
                       35.7, 3.1, 152.3, 5.4, 0.7, // new measurements
                       39.4, 5.4, 150.3, 5.5);     // old measurements

    // marginalize
    m.MarginalizeAll();

    // find mode
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized("CombinationModel_plots.pdf");

    // print summary plots
    m.PrintParameterPlot("CombinationModel_parameters.pdf");
    m.PrintCorrelationPlot("CombinationModel_correlation.pdf");

    // m.SetKnowledgeUpdateDrawingStyle(BCModel::kKnowledgeUpdateDetailedPosterior);
    m.PrintKnowledgeUpdatePlots("CombinationModel_update.pdf");

    // print results of the analysis the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;

}

