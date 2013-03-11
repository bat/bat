#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "ReferenceCounting.h"

int main()
{

   // run comparison between options of evaluation
   bool flag_comparison = false;

	 // run a superfine grid
	 bool flag_superfine = true;

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new ReferenceCounting object
   ReferenceCounting * m = new ReferenceCounting();

	 // BAT settings
	 m->SetNbins("s", 100);
	 m->SetNbins("b", 100);
	 m->MCMCSetPrecision(BCIntegrate::kMedium);

	 // set option of how to evaluate prior
	 m->SetPriorEvalOption(ReferenceCounting::kHistogram);
	 
	 // set background expectation
	 double bkg_exp = 10; // expectation value
	 double bkg_std = 5;  // uncertainty on background

	 double alpha   = bkg_exp*bkg_exp/bkg_std/bkg_std; 
	 double beta    = bkg_exp/bkg_std/bkg_std;

	 m->SetAlphaBeta(alpha, beta);

	 // set number of observed events
	 m->SetNObs(20);

	 // set parameter range
	 m->SetParameterRange(0, 0.0, 50); // signal 
	 m->SetParameterRange(1, 0.0, 35); // background

   // create a new summary tool object
   BCSummaryTool * summary = new BCSummaryTool(m);

	 // perform sampling with MCMC
	 m->MarginalizeAll();

	 // perform minimization with Minuit
	 m->FindMode( m->GetBestFitParameters() );

   // draw all marginalized distributions into a pdf file
	 m->PrintAllMarginalized("ReferenceCounting_plots.pdf");

   // print individual histograms
	 m->GetMarginalized("s")->Print("ReferenceCounting_s.pdf", "BTulB3CS1D0pdf0L");
	 m->GetMarginalized("s")->Print("ReferenceCounting_s_logy.pdf", "BTulB3CS1D0pdf0Llogy");
	 m->GetMarginalized("b")->Print("ReferenceCounting_b.pdf");
	 m->GetMarginalized("b")->Print("ReferenceCounting_b_logy.pdf", "BTulB3CS1D0pdf0Llogy");
	 m->GetMarginalized("s", "b")->Print("ReferenceCounting_sb.pdf", "BTfB3CS1meangmodelmode");

	 /*
	 BCH1D* hist_s = new BCH1D( m->GetSlice("s", "b", m->GetBestFitParameters(), 1000)->GetHistogram()->ProjectionX() );
	 hist_s->Print("s_proj.pdf", "BTulB3CS1D0pdf0L");
	 hist_s->Print("s_proj_logy.pdf", "BTulB3CS1D0pdf0Llogy");
	 BCH1D* hist_b = new BCH1D( m->GetSlice("s", "b", m->GetBestFitParameters(), 1000)->GetHistogram()->ProjectionY() );
	 hist_b->Print("b_proj.pdf", "BTulB3CS1D0pdf0L");
	 hist_b->Print("b_proj_loy.pdf", "BTulB3CS1D0pdf0Llogy");
	 */

   // print summary plots
	 summary->PrintKnowledgeUpdatePlots("ReferenceCounting_update.ps");

   // print results of the analysis into a text file
	 m->PrintResults("ReferenceCounting_results.txt");

	 delete summary;

	 if (flag_comparison) {
		 
		 TH1D* hist_slice = m->GetSlice("s", "b", m->GetBestFitParameters())->GetHistogram()->ProjectionX();
		 hist_slice->Scale(1.0/hist_slice->Integral());
		 hist_slice->SetLineColor(kBlack);
		 hist_slice->SetStats(kFALSE);

		 m->SetPriorEvalOption(ReferenceCounting::kAnalytic);
		 m->MarginalizeAll();
		 m->FindMode( m->GetBestFitParameters() );
		 TH1D* hist_analytic = new TH1D(*(m->GetMarginalized("s")->GetHistogram()));
		 hist_analytic->Scale(1.0/hist_analytic->Integral());
		 hist_analytic->SetLineColor(kGreen);

		 m->SetPriorEvalOption(ReferenceCounting::kHistogram);
		 m->MarginalizeAll();
		 m->FindMode( m->GetBestFitParameters() );
		 TH1D* hist_histogram = new TH1D(*(m->GetMarginalized("s")->GetHistogram()));
		 hist_histogram->Scale(1.0/hist_histogram->Integral());
		 hist_histogram->SetLineColor(kRed);

		 TCanvas* c1 = new TCanvas();

		 hist_slice->Draw("");
		 hist_analytic->Draw("SAME");
		 hist_histogram->Draw("SAME");
		 c1->Print("comparison.pdf");

		 c1->SetLogy();
		 c1->Print("comparison_log.pdf");
	 }

	 if (flag_superfine) {
		 BCH1D* hist_s = new BCH1D( m->GetSlice("s", "b", m->GetBestFitParameters(), 3000)->GetHistogram()->ProjectionX() );
		 hist_s->Print("superfine_s_proj.pdf", "BTulB3CS1D0pdf0L");
		 hist_s->Print("superfine_s_proj_logy.pdf", "BTulB3CS1D0pdf0Llogy");

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
	 }

   // close log file
   BCLog::CloseLog();

	 delete m;

   return 0;

}

