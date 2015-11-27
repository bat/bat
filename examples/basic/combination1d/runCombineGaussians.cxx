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
    // new measurement = (35.7 +- 3.1) MeV
    // old measurement = (39.4 +- 5.4) MeV
    CombinationModel m("Combination Model", 35.7, 3.1, 39.4, 5.4);

    // marginalize
    m.MarginalizeAll();

    // find mode
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized("CombinationModel_plots.pdf");

    // print knowledge update plot, with detailed posterior
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    m.PrintKnowledgeUpdatePlots("CombinationModel_update.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;

}

