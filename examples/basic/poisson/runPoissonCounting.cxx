#include "PoissonModel.h"

#include <BAT/BCLog.h>

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new PoissonModel object
    PoissonModel m("Poison Model");

    // set number of observed events
    m.SetNObs(7);

    // run MCMC and marginalize posterior wrt. all parameters
    // and all combinations of two parameters
    m.MarginalizeAll();

    // find mode starting from the best fit parameters
    m.FindMode();

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    // no error
    return 0;
}
