#include <BAT/BCLog.h>

#include "GaussModel.h"

#include <string>

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new GaussModel object
    GaussModel m("Gauss Model");

    // set marginalization method
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set MCMC precision
    m.SetPrecision(BCEngineMCMC::kMedium);

    // switch writing of Markov Chains on
    m.WriteMarkovChain(m.GetSafeName() + ".root", "RECREATE");

    // run MCMC and marginalize posterior wrt. all parameters
    // and all combinations of two parameters
    m.MarginalizeAll();

    // run MCMC for another 100,000 iterations without rerunning the pre-run
    m.SetNIterationsRun(100000);
    m.MarginalizeAll();

    // find the mode, starting from best fit point found by the MCMC
    m.FindMode();

    // draw all marginalized distributions into a PostScript file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // write marginalized distributions to output file
    m.WriteMarginalizedDistributions(m.GetSafeName() + ".root", "UPDATE");

    // close log file
    BCLog::CloseLog();

    return 0;
}

