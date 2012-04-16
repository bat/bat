
void testSplusBmtf()
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
   TFile * file = new TFile(fname.c_str(), "READ");

   // check if file is open
   if (!file->IsOpen()) {
      BCLog::OutError(Form("Could not open file %s.",fname.c_str()));
      BCLog::OutError("Run macro CreateHistograms.C in Root to create the file.");
      return;
   }

   // read histograms
   TH1D hist_signal     = *(TH1D *)file->Get("hist_sgn");   // signal template
   TH1D hist_background = *(TH1D *)file->Get("hist_bkg");   // background template
   TH1D hist_data1      = *(TH1D *)file->Get("hist_data1"); // data for channel 1
   TH1D hist_data2      = *(TH1D *)file->Get("hist_data2"); // data for channel 2

   // ---- perform fitting ---- //

   // create new fitter object
   BCMTF * m = new BCMTF();

   // set the required precision of the MCMC (kLow, kMedium, kHigh)
   // the higher the precision the longer the MCMC run
   m->MCMCSetPrecision(BCEngineMCMC::kLow);

   // add channels
   m->AddChannel("channel1");
   m->AddChannel("channel2");

   // add processes
   m->AddProcess("background", 150., 500.);
   m->AddProcess("signal",       0., 200.);

   // set data
   m->SetData("channel1", hist_data1);
   m->SetData("channel2", hist_data2);

   // set template and histograms
   m->SetTemplate("channel1", "signal",     hist_signal,     1.0);
   m->SetTemplate("channel1", "background", hist_background, 1.0);

   m->SetTemplate("channel2", "signal",     hist_signal,     0.5);
   m->SetTemplate("channel2", "background", hist_background, 1.0);

   // set priors
   m->SetPriorGauss("background", 300., 30.);
   m->SetPriorConstant("signal");

   // run MCMC
   m->MarginalizeAll();

   // find global mode
   m->FindMode( m->GetBestFitParameters() );

   // print all marginalized distributions
   m->PrintAllMarginalized("marginalized.ps");

   // print results of the analysis into a text file
   m->PrintResults("results.txt");

   // print templates and stacks
   for (int i = 0; i < m->GetNChannels(); ++i) {
      BCMTFChannel * channel = m->GetChannel(i);
      channel->PrintTemplates(Form("%s_templates.ps", channel->GetName().c_str()));
      m->PrintStack(i, m->GetBestFitParameters(), Form("%s_stack.eps", channel->GetName().c_str()), "");
   }

   // ---- clean up ---- //

   // free memory
   delete m;

}
