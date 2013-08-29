#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include <TMatrixT.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>

#include "MyFit.h"

int main(int argc, char *argv[])
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new MyFit object
   MyFit * m = new MyFit();

   // set precision
   m->MCMCSetPrecision(BCIntegrate::kMedium);
   
   if (argc==2) {
     int isopen = m->ReadInput(argv[1]);
     if (!isopen) {
       std::cout << "Could not open file. Exit." << std::endl;
       return 1;
     }
   }
   else {
     std::cout << "No input file specified. Exit." << std::endl;
     return 1;
   }

	 TCanvas* c1 = new TCanvas();
	 c1->cd();

	 // full 2D analysis

	 m->MarginalizeAll();
	 m->FindMode( m->GetBestFitParameters() );
	 m->PrintAllMarginalized("marginalized_all.pdf");
	 m->PrintResults("results_all.txt");
	 BCH2D* hist_all = m->GetSlice("#kappa_{1}", "#kappa_{2}", m->GetBestFitParameters(), 200);
	 hist_all->Draw("BTcB3CS1gmode");
	 c1->SetGrid(1,1);
	 c1->Print("2D_all.pdf");

	 // xs only

	 m->GetMeasurement(0)->SetFlagActive(true);
	 m->GetMeasurement(1)->SetFlagActive(true);
	 m->GetMeasurement(2)->SetFlagActive(false);
	 m->GetMeasurement(3)->SetFlagActive(false);
	 m->PrepareAnalysis();

	 m->MarginalizeAll();
	 m->FindMode( m->GetBestFitParameters() );
	 m->PrintAllMarginalized("marginalized_xs.pdf");
	 m->PrintResults("results_xs.txt");
	 BCH2D* hist_xs = m->GetSlice("#kappa_{1}", "#kappa_{2}", m->GetBestFitParameters(), 200);
	 hist_xs->Draw("BTcB3CS1gmode");
	 c1->SetGrid(1,1);
	 c1->Print("2D_xs.pdf");

	 // br only 

	 m->GetMeasurement(0)->SetFlagActive(false);
	 m->GetMeasurement(1)->SetFlagActive(false);
	 m->GetMeasurement(2)->SetFlagActive(true);
	 m->GetMeasurement(3)->SetFlagActive(true);
	 m->PrepareAnalysis();

	 m->MarginalizeAll();
	 m->FindMode( m->GetBestFitParameters() );
	 m->PrintAllMarginalized("marginalized_br.pdf");
	 m->PrintResults("results_br.txt");
	 BCH2D* hist_br = m->GetSlice("#kappa_{1}", "#kappa_{2}", m->GetBestFitParameters(), 200);
	 hist_br->Draw("BTcB3CS1gmode");
	 c1->SetGrid(1,1);
	 c1->Print("2D_br.pdf");

  // clean up
   delete m;
   
   // close log file
   BCLog::CloseLog();
   
   // no error
   return 0;
}

