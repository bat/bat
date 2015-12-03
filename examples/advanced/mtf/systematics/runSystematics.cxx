#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

#include <BAT/BCMTFAnalysisFacility.h>
#include <BAT/BCMTF.h>
#include <BAT/BCMTFChannel.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

int main()
{
    // ---- set style and open log files ---- //

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // ---- read histograms from a file ---- //

    // open file
    std::string fname = "templates.root";
    TFile* file = TFile::Open(fname.c_str(), "READ");

    // check if file is open
    if (!file || file->IsZombie()) {
        BCLog::OutError(Form("Could not open file %s.", fname.c_str()));
        BCLog::OutError("Run macro CreateHistograms.C in Root to create the file.");
        return 1;
    }

    // read template histograms
    TH1D* hist_signal      = (TH1D*) file->Get("hist_sgn");    // signal template
    TH1D* hist_background  = (TH1D*) file->Get("hist_bkg");    // background template
    // read data histograms
    TH1D* hist_data1       = (TH1D*) file->Get("hist_data1");  // data for channel 1
    TH1D* hist_data2       = (TH1D*) file->Get("hist_data2");  // data for channel 2
    // read systematic histograms
    TH1D* hist_syst1_bkg_1 = (TH1D*) file->Get("hist_syst1_bkg_1");
    TH1D* hist_syst1_sgn_1 = (TH1D*) file->Get("hist_syst1_sgn_1");
    TH1D* hist_syst1_bkg_2 = (TH1D*) file->Get("hist_syst1_bkg_2");
    TH1D* hist_syst1_sgn_2 = (TH1D*) file->Get("hist_syst1_sgn_2");
    TH1D* hist_syst2_bkg_1 = (TH1D*) file->Get("hist_syst2_bkg_1");
    TH1D* hist_syst2_sgn_1 = (TH1D*) file->Get("hist_syst2_sgn_1");
    TH1D* hist_syst2_bkg_2 = (TH1D*) file->Get("hist_syst2_bkg_2");
    TH1D* hist_syst2_sgn_2 = (TH1D*) file->Get("hist_syst2_sgn_2");

    if (!hist_signal || !hist_background || !hist_data1 || !hist_data2
            || !hist_syst1_bkg_1 || !hist_syst1_sgn_1
            || !hist_syst1_bkg_2 || !hist_syst1_sgn_2
            || !hist_syst2_bkg_1 || !hist_syst2_sgn_1
            || !hist_syst2_bkg_2 || !hist_syst2_sgn_2) {
        BCLog::OutError("Could not open data histograms.");
        return 1;
    }

    // ---- perform fitting ---- //

    // create new fitter object
    BCMTF m;


    // set Metropolis as marginalization method
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set the required precision of the MCMC (kLow, kQuick, kMedium, kHigh)
    // the higher the precision the longer the MCMC run
    m.SetPrecision(BCEngineMCMC::kQuick);

    // add channels
    m.AddChannel("channel1");
    m.AddChannel("channel2");

    // add processes
    m.AddProcess("background", 150., 500.);
    m.AddProcess("signal",       0., 200.);

    // add systematics
    m.AddSystematic("systematic1", -5., 5.);
    m.AddSystematic("systematic2", -5., 5.);

    // set data
    m.SetData("channel1", *hist_data1);
    m.SetData("channel2", *hist_data2);

    // set template histograms
    m.SetTemplate("channel1", "background", *hist_background, 1.0);
    m.SetTemplate("channel1", "signal",     *hist_signal,     1.0);

    m.SetTemplate("channel2", "background", *hist_background, 1.0);
    m.SetTemplate("channel2", "signal",     *hist_signal,     0.5);

    // set systematic histograms
    m.SetSystematicVariation("channel1", "background", "systematic1", 0.1, 0.1);
    m.SetSystematicVariation("channel1", "signal", "systematic1", *hist_syst1_sgn_1, *hist_syst1_sgn_1);
    m.SetSystematicVariation("channel2", "background", "systematic1", 0.1, 0.1);
    m.SetSystematicVariation("channel2", "signal", "systematic1", *hist_syst1_sgn_2, *hist_syst1_sgn_2);

    m.SetSystematicVariation("channel1", "background", "systematic2", *hist_syst2_bkg_1, *hist_syst2_bkg_1);
    m.SetSystematicVariation("channel1", "signal", "systematic2", *hist_syst2_sgn_1, *hist_syst2_sgn_1);
    m.SetSystematicVariation("channel2", "background", "systematic2", *hist_syst2_bkg_2, *hist_syst2_bkg_2);
    m.SetSystematicVariation("channel2", "signal", "systematic2", *hist_syst2_sgn_2, *hist_syst2_sgn_2);

    // set priors
    m.GetParameter("signal").SetPriorConstant();
    m.GetParameter("background").SetPrior(new BCGaussianPrior(300., 30.));
    m.GetParameter("systematic1").SetPrior(new BCGaussianPrior(0., 1.));
    m.GetParameter("systematic2").SetPrior(new BCGaussianPrior(0., 1.));

    // run MCMC
    m.MarginalizeAll();

    // find global mode
    m.FindMode(m.GetBestFitParameters());

    // print all marginalized distributions
    m.PrintAllMarginalized("marginalized.pdf");

    // print results of the analysis into the log
    m.PrintSummary();

    // print summary results
    m.PrintParameterPlot("summary_parameters.pdf");
    m.PrintCorrelationPlot("summary_correlationplot.pdf");
    m.PrintCorrelationMatrix("summary_correlationmatrix.pdf");
    m.PrintKnowledgeUpdatePlots("summary_update.pdf");

    // print templates and stacks
    for (int i = 0; i < m.GetNChannels(); ++i) {
        BCMTFChannel* channel = m.GetChannel(i);
        channel->PrintTemplates(channel->GetSafeName() + "_templates.pdf");
        channel->PrintTemplate(0, Form("background_%i.pdf", i));
        channel->PrintTemplate(1, Form("signal_%i.pdf", i));

        m.PrintStack(i, m.GetBestFitParameters(), channel->GetSafeName() + "_stack.pdf");
    }

    // ---- perform single systematic analysis ---- //

    // create new analysis facility
    BCMTFAnalysisFacility facility(&m);

    // perform analysis
    facility.PerformSingleSystematicAnalyses("systematics");

    // ---- clean up ---- //

    // close log file
    BCLog::CloseLog();

    return 0;
}
