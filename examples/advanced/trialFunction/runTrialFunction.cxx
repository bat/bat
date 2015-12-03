#include <BAT/BCAux.h>
#include <BAT/BCLog.h>

#include "GaussModel.h"

int main()
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new GaussModel object
    GaussModel m("gausMod");

    // set marginalization method
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set MCMC precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // set scale factor upper limit
    m.SetScaleFactorUpperLimit(10);

    // run MCMC and marginalize posterior wrt. all parameters
    // and all combinations of two parameters
    m.MarginalizeAll();

    // if MCMC was run before (MarginalizeAll()) it is
    // possible to use the mode found by MCMC as
    // starting point of Minuit minimization
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized("GaussModel_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;

}

