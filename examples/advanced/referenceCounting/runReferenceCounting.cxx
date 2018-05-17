#include "ReferenceCounting.h"

#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCIntegrate.h>
#include <BAT/BCMath.h>
#include <BAT/BCModel.h>

#include <iostream>

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // cache factorials that may be used by Reference Prior
    BCMath::CacheFactorials(1000);

    // create new ReferenceCounting object for 20 measured events
    ReferenceCounting m("Reference Counting", 20);

    // set background expectation
    double bkg_exp = 10; // expectation value
    double bkg_std = 5;  // uncertainty on background

    double shape   = bkg_exp * bkg_exp / bkg_std / bkg_std;
    double rate    = bkg_exp / bkg_std / bkg_std;

    // set option of how to evaluate prior
    // uncomment one to choose method, kAnalytic is the slowest
    //    m.SetPrior(ReferenceCounting::kAnalytic, shape, rate);
    m.SetPrior(ReferenceCounting::kHistogram, shape, rate);
    //    m.SetPrior(ReferenceCounting::kApprox, shape, rate);

    // set precision level
    m.SetPrecision(BCEngineMCMC::kQuick);

    // perform sampling with MCMC
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // perform minimization with Minuit
    m.FindMode();

    // draw all marginalized distributions into a pdf file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");
    m.GetBCH1DdrawingOptions().SetLogy();
    m.PrintAllMarginalized(m.GetSafeName() + "_logy.pdf");

    // print knowledge update with detailed posteriors
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    m.GetBCH2DPosteriorDrawingOptions().SetNSmooth(0);
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis to the log
    m.PrintSummary();

    // print slice results to log
    m.GetMarginalized("s").PrintSummary();

    // close log file
    BCLog::CloseLog();

    return 0;
}

