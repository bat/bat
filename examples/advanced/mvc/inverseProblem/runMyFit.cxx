#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <fstream>

#include "MyFit.h"

int main(int argc, char* argv[])
{

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // create new MyFit object
    MyFit* m = new MyFit();

    // set Metropolis as marginalization method
    m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    m->MCMCSetPrecision(BCIntegrate::kMedium);

    if (argc == 2) {
        int isopen = m->ReadInput(argv[1]);
        if (!isopen) {
            BCLog::OutError("Could not open file. Exit.");
            return 1;
        }
    } else {
        BCLog::OutError("No input file specified. Exit.");
        return 1;
    }

    // full 2D analysis
    m->MarginalizeAll();
    m->FindMode(m->GetGlobalMode());
    m->PrintAllMarginalized("MyFit_plots.pdf");
    m->PrintResults("MyFit_results.txt");

    // clean up
    delete m;

    // close log file
    BCLog::CloseLog();

    // no error
    return 0;
}

