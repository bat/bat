//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x singleChannel.C
//
// or
//
//    root[1] .L singleChannel.C
//    root[2] singleChannel()
//
// or from the command line
//
//    $ root singleChannel.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x singleChannel.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro
// the #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#define COMPILER (!defined(__CINT__) && !defined(__CLING__))

#if defined(__MAKECINT__) || defined(__ROOTCLING__) || COMPILER

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCMTF.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCParameter.h>

#include <TFile.h>
#include <TH1D.h>

#include <string>

#endif

void singleChannel()
{
    // ---- set style and open log files ---- //

    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // ---- read histograms from a file ---- //

    // open file
    std::string fname = "templates.root";
    TFile* file = TFile::Open(fname.c_str(), "READ");

    // check if file is open
    if (!file || !file->IsOpen()) {
        BCLog::OutError(Form("Could not open file %s.", fname.c_str()));
        BCLog::OutError("Run macro CreateHistograms.C in Root to create the file.");
        return;
    }

    // read histograms
    TH1D* hist_signal     = (TH1D*)file->Get("hist_sgn");    // signal template
    TH1D* hist_background = (TH1D*)file->Get("hist_bkg");    // background template
    TH1D* hist_data       = (TH1D*)file->Get("hist_data");   // data

    if (!hist_signal || !hist_background || !hist_data) {
        BCLog::OutError("Could not find data histograms");
        return;
    }

    // ---- perform fitting ---- //

    // create new fitter object
    BCMTF* m = new BCMTF("SingleChannelMTF");

    // add channels
    m->AddChannel("channel1");

    // add processes
    m->AddProcess("background", 200., 400.);
    m->AddProcess("signal",       0., 200.);

    // set data
    m->SetData("channel1", *hist_data);

    // set template and histograms
    m->SetTemplate("channel1", "signal",     *hist_signal,     1.0);
    m->SetTemplate("channel1", "background", *hist_background, 1.0);

    // set priors
    m->GetParameter("background").SetPrior(new BCGaussianPrior(300., 10.));
    m->GetParameter("signal").SetPriorConstant();

    // set precision
    m->SetPrecision(BCEngineMCMC::kQuick);

    // marginalize
    m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // find global mode
    m->FindMode( m->GetBestFitParameters() );

    // print all marginalized distributions
    m->PrintAllMarginalized("marginalized.pdf");

    // print results of the analysis into a text file
    m->PrintSummary();

    // print templates and stacks
    BCMTFChannel* channel = m->GetChannel(0);
    channel->PrintTemplates(channel->GetSafeName() + "_templates.pdf");
    m->PrintStack(0, m->GetBestFitParameters(), channel->GetSafeName() + "_stack.pdf");

    // ---- clean up ---- //

    // free memory
    delete m;

}
