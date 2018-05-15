#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "GaussModel.h"

#include <string>

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new GaussModel object
    GaussModel m("Continue Gauss Model");

    // load in previous running from file
    m.PrepareToContinueMarginalization("GaussModel.root");

    // write out new samples to a new file
    m.WriteMarkovChain(m.GetSafeName() + ".root", "RECREATE");

    // continue marginalization
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // find the new mode
    m.FindMode(m.GetBestFitParameters());

    // print new marginalization plots
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;
}

