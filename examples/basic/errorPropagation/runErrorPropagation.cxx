#include <BAT/BCLog.h>

#include "RatioModel.h"

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new RatioModel object
    RatioModel m("Ratio Model");

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // run the MCMC and marginalize w.r.t. to all parameters
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // find mode using the best fit parameters as start values
    m.FindMode();

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    // no error
    return 0;

}

