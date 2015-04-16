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
    CombinationModel* m = new CombinationModel("combMod");

    // marginalize
    m->MarginalizeAll();

    // find mode
    m->FindMode(m->GetGlobalMode());

    // draw all marginalized distributions into a PostScript file
    m->PrintAllMarginalized("CombinationModel_plots.pdf");

    // print all summary plots
    m->PrintParameterPlot("CombinationModel_parameters.pdf");
    m->PrintCorrelationPlot("CombinationModel_correlation.pdf");
    m->PrintKnowledgeUpdatePlots("CombinationModel_update.pdf");

    // print results of the analysis into a text file
    m->PrintResults("CombinationModel_results.txt");

    // close log file
    BCLog::CloseLog();

    delete m;

    return 0;

}

