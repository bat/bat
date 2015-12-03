#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "BinomialModel.h"

int main()
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new BinomialModel object
    BinomialModel m("Binomial Model", 20, 3);

    // run the MCMC and calculate the marginalized distributions
    m.MarginalizeAll();

    // find mode starting from the best fit parameters
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a pdf file
    m.PrintAllMarginalized("BinomialModel_plots.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // close log file
    BCLog::CloseLog();

    // no errors
    return 0;

}

