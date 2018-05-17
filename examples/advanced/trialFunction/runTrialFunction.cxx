#include <BAT/BCLog.h>

#include "GaussModel.h"

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new GaussModel object
    GaussModel m("Gauss Model");

    // set MCMC precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // set scale factor upper limit
    m.SetScaleFactorUpperLimit(10);

    // run MCMC and marginalize posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // Optimize, starting from best fit point found by the marginalization
    m.FindMode();

    // draw all marginalized distributions into a file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;
}

