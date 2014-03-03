#include "ReferenceCounting.h"

#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>
#include <BAT/BCSummaryTool.h>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>

int main()
{
   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new ReferenceCounting object
   ReferenceCounting* m = new ReferenceCounting();

   // set option of how to evaluate prior
   // kHistogram : calculate prior first and fill into histogram
   // kAnalytic  : calculate analytic expression
   // kApprox    : calculate prior from a TF1 approximation fitted to a histgram
   m->SetPriorEvalOption(ReferenceCounting::kAnalytic);

   // set background expectation
   double bkg_exp = 10; // expectation value
   double bkg_std = 5;  // uncertainty on background

   double alpha   = bkg_exp*bkg_exp/bkg_std/bkg_std;
   double beta    = bkg_exp/bkg_std/bkg_std;

   m->SetAlphaBeta(alpha, beta);

   // set number of observed events
   m->SetNObs(20);

   // set parameter range
   m->GetParameter("s")->SetLimits(0.0, 50); // signal
   m->GetParameter("b")->SetLimits(0.0, 35); // background

   // create a new summary tool object
   BCSummaryTool * summary = new BCSummaryTool(m);

   // perform sampling with MCMC
   m->MarginalizeAll();

   // perform minimization with Minuit
   m->FindMode( m->GetBestFitParameters() );

   // draw all marginalized distributions into a pdf file
   m->PrintAllMarginalized("ReferenceCounting_plots.pdf");

   // print individual histograms
   BCH1D* hist_s = m->GetMarginalized("s");
   hist_s->Print("ReferenceCounting_s.pdf", "BTulB3CS1D0pdf0L");
   hist_s->Print("ReferenceCounting_s_logy.pdf", "BTulB3CS1D0pdf0Llogy");

   BCH1D* hist_b = m->GetMarginalized("b");
   hist_b->Print("ReferenceCounting_b.pdf");
   hist_b->Print("ReferenceCounting_b_logy.pdf", "BTulB3CS1D0pdf0Llogy");

   BCH2D* hist_sb = m->GetMarginalized("s", "b");
   hist_sb->Print("ReferenceCounting_sb.pdf", "BTfB3CS1meangmodelmode");

   // print priors
   m->PrintPriors("priors.pdf");

   // print summary plots
   summary->PrintKnowledgeUpdatePlots("ReferenceCounting_update.pdf");

   // free memory
   delete summary;

   // print results of the analysis into a text file
   m->PrintResults("ReferenceCounting_results.txt");

   // print slice results to screen

   std::cout << " Mean   : " << hist_s->GetMean() << std::endl;
   std::cout << " Median : " << hist_s->GetMedian() << std::endl;
   std::cout << " Mode   : " << hist_s->GetMode() << std::endl;
   std::cout << " Std    : " << hist_s->GetSTD() << std::endl;
   std::cout << " Var    : " << hist_s->GetVariance() << std::endl;
   std::cout << " Q  5%  : " << hist_s->GetQuantile(0.05) << std::endl;
   std::cout << " Q 10%  : " << hist_s->GetQuantile(0.10) << std::endl;
   std::cout << " Q 16%  : " << hist_s->GetQuantile(0.16) << std::endl;
   std::cout << " Q 50%  : " << hist_s->GetQuantile(0.50) << std::endl;
   std::cout << " Q 84%  : " << hist_s->GetQuantile(0.84) << std::endl;
   std::cout << " Q 90%  : " << hist_s->GetQuantile(0.90) << std::endl;
   std::cout << " Q 95%  : " << hist_s->GetQuantile(0.95) << std::endl;

   // close log file
   BCLog::CloseLog();

   // free memory
   delete m;

   return 0;

}

