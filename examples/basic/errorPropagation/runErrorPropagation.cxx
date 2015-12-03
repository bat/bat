#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "RatioModel.h"

int main()
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new RatioModel object
    RatioModel m("ratMod");

    // set Metropolis as marginalization method
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // run the MCMC and marginalize w.r.t. to all parameters
    m.MarginalizeAll();

    // find mode using the best fit parameters as start values
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized("RatioModel_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    // no error
    return 0;

}

