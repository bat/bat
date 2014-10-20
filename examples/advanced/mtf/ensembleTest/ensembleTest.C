//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x ensembleTest.C
//
// or
//
//    root[1] .L ensembleTest.C
//    root[2] ensembleTest()
//
// or from the command line
//
//    $ root ensembleTest.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x ensembleTest.C+
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
#include <BAT/BCLog.h>
#include <BAT/BCMTF.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCMTFAnalysisFacility.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <string>

#endif

void ensembleTest()
{
   // ---- set style and open log files ---- //

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // ---- read histograms from a file ---- //

   // open file
   std::string fname = "templates.root";
   TFile * file = TFile::Open(fname.c_str(), "READ");

   // check if file is open
   if (!file->IsOpen()) {
      BCLog::OutError(Form("Could not open file %s.",fname.c_str()));
      BCLog::OutError("Run macro CreateHistograms.C in Root to create the file.");
      return;
   }

   // read histograms
   TH1D hist_bkg1  = *(TH1D *)file->Get("hist_bkg1");   // background template for channel 1
   TH1D hist_bkg2  = *(TH1D *)file->Get("hist_bkg2");   // background template for channel 2
   TH1D hist_sgn1  = *(TH1D *)file->Get("hist_sgn1");   // signal template for channel 1
   TH1D hist_sgn2  = *(TH1D *)file->Get("hist_sgn2");   // signal template for channel 2
   TH1D hist_data1 = *(TH1D *)file->Get("hist_data1"); // data for channel 1
   TH1D hist_data2 = *(TH1D *)file->Get("hist_data2"); // data for channel 2

   // ---- perform fitting ---- //

   // create new fitter object
   BCMTF * m = new BCMTF("MyModel");

   // set Metropolis as marginalization method
   m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

   // set the required precision of the MCMC (kLow, kMedium, kHigh)
   // the higher the precision the longer the MCMC run
   m->MCMCSetPrecision(BCEngineMCMC::kMedium);

   // add channels
   m->AddChannel("channel1");
   m->AddChannel("channel2");

   // add processes
   m->AddProcess("background_channel1", 700., 900.);
   m->AddProcess("background_channel2", 300., 700.);
   m->AddProcess("signal",       0., 400.);

   // set data
   m->SetData("channel1", hist_data1);
   m->SetData("channel2", hist_data2);

   // set template and histograms
   // note: the process "background_channel2" is ignored in channel 1
   m->SetTemplate("channel1", "signal", hist_sgn1, 0.5);
   m->SetTemplate("channel1", "background_channel1", hist_bkg1, 1.0);
//   m->SetTemplate("channel1", "background_channel2", hist_bkg1, 0.0);

   // note: the process "background_channel1" is ignored in channel 2
   m->SetTemplate("channel2", "signal", hist_sgn2, 1.0);
   m->SetTemplate("channel2", "background_channel2", hist_bkg2, 1.0);
//   m->SetTemplate("channel2", "background_channel1", hist_bkg2, 0.0);

   // set priors
   m->SetPriorGauss("background_channel1", 800., 10.);
   m->SetPriorGauss("background_channel2", 500., 50.);
   m->SetPriorConstant("signal");

   // run MCMC
   m->MarginalizeAll();

   // find global mode
   m->FindMode( m->GetBestFitParameters() );

   // print all marginalized distributions
   m->PrintAllMarginalized("marginalized.pdf");

   // print results of the analysis into a text file
   m->PrintResults("results.txt");

   // print templates and stacks
   for (int i = 0; i < m->GetNChannels(); ++i) {
      BCMTFChannel * channel = m->GetChannel(i);
      channel->PrintTemplates(Form("%s_templates.pdf", channel->GetName().c_str()));
      m->PrintStack(i, m->GetBestFitParameters(), Form("%s_stack.pdf", channel->GetName().c_str()));
   }

   // ---- perform ensemble tests ---- //

   // create new analysis facility
   BCMTFAnalysisFacility * facility = new BCMTFAnalysisFacility(m);

   // settings
   facility->SetFlagMarginalize(false);
//   facility->SetFlagMarginalize(true);

   // open new file
   file = TFile::Open("ensembles.root", "RECREATE");
   file->cd();

   // create ensembles
   TTree * tree = facility->BuildEnsembles( m->GetBestFitParameters(), 2000 );

   // run ensemble test
   TTree * tree_out = facility->PerformEnsembleTest(tree, 2000);

   // write trees into file
   tree->Write();
   tree_out->Write();

   // close file
   file->Close();

   // free memory
   delete file;

   // ---- clean up ---- //

   // free memory
   delete m;

}
