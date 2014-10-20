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
#include <BAT/BCLog.h>
#include <BAT/BCMTF.h>
#include <BAT/BCMTFAnalysisFacility.h>
#include <BAT/BCMTFChannel.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <string>

#endif


void mcstat()
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
   TH1D hist_signal     = *(TH1D *)file->Get("hist_sgn");   // signal template
   TH1D hist_background = *(TH1D *)file->Get("hist_bkg");   // background template
   TH1D hist_data       = *(TH1D *)file->Get("hist_data");  // data

   // ---- perform fitting ---- //

   // create new fitter object
   BCMTF * m = new BCMTF("SingleChannelMTF");

   // add channels
   m->AddChannel("channel1");

   // add processes
   m->AddProcess("background", 200., 400.);
   m->AddProcess("signal",       0., 200.);

   // set data
   m->SetData("channel1", hist_data);

   // set template and histograms
   m->SetTemplate("channel1", "signal",     hist_signal,     1.0);
   m->SetTemplate("channel1", "background", hist_background, 1.0);

   // set priors
   m->SetPriorGauss("background", 300., 10.);
   m->SetPriorConstant("signal");

   // print templates
	 m->GetChannel(0)->PrintTemplates(Form("%s_templates.pdf", m->GetChannel(0)->GetName().c_str()));

   // ---- perform ensemble tests ---- //

   // create new analysis facility
   BCMTFAnalysisFacility * facility = new BCMTFAnalysisFacility(m);

   // settings
   facility->SetFlagMarginalize(true);

   // open new file
   file = TFile::Open("ensembles.root", "RECREATE");
   file->cd();

   // create ensembles; option "data" means that all ensembles equal the data set
   TTree * tree = facility->BuildEnsembles( std::vector<double>(0), 1000, "data");

   // run ensemble test; option "MCP" means that the templates are flucutated via a Poisson model
   TTree * tree_out = facility->PerformEnsembleTest(tree, 1000, 0, "MCP");

   // write trees into file
   tree->Write();
   tree_out->Write();

   // close file
   file->Close();

   // ---- clean up ---- //

   // free memory
   delete file;

   // free memory
   delete m;

}
