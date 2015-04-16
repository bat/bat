#include "ReferenceCounting.h"

#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>

int main()
{
    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // create new ReferenceCounting object
    ReferenceCounting* m = new ReferenceCounting("refcountMod");

    // set option of how to evaluate prior
    // kHistogram : calculate prior first and fill into histogram
    // kAnalytic  : calculate analytic expression
    // kApprox    : calculate prior from a TF1 approximation fitted to a histgram
    m->SetPriorEvalOption(ReferenceCounting::kAnalytic);

    // set background expectation
    double bkg_exp = 10; // expectation value
    double bkg_std = 5;  // uncertainty on background

    double alpha   = bkg_exp * bkg_exp / bkg_std / bkg_std;
    double beta    = bkg_exp / bkg_std / bkg_std;

    m->SetAlphaBeta(alpha, beta);

    // set number of observed events
    m->SetNObs(20);

    // set parameter range
    m->GetParameter("s").SetLimits(0.0, 50); // signal
    m->GetParameter("b").SetLimits(0.0, 35); // background

    // perform sampling with MCMC
    m->MarginalizeAll();

    // perform minimization with Minuit
    m->FindMode(m->GetGlobalMode());

    // draw all marginalized distributions into a pdf file
    m->PrintAllMarginalized("ReferenceCounting_plots.pdf");

    // print histograms
    m->PrintAllMarginalized("ReferenceCounting.pdf");
    m->GetBCH1DdrawingOptions()->SetLogy();
    m->PrintAllMarginalized("ReferenceCounting_logy.pdf");

    // print priors
    m->PrintPriors("priors.pdf");

    // print summary plots
    m->PrintKnowledgeUpdatePlots("ReferenceCounting_update.pdf");

    // print results of the analysis into a text file
    m->PrintResults("ReferenceCounting_results.txt");

    // print slice results to screen
    m->GetMarginalized("s")->PrintToStream(std::cout, " ");

    // close log file
    BCLog::CloseLog();

    // free memory
    delete m;

    return 0;

}

