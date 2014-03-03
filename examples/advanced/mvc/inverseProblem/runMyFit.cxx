#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include <TMatrixT.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h> 

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

  // set Metropolis as marginalization method
  m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

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

  // full 2D analysis
  m->MarginalizeAll();
  m->FindMode( m->GetBestFitParameters() );
  m->PrintAllMarginalized("MyFit_plots.pdf");
  m->PrintResults("MyFit_results.txt");

  BCH2D* hist_all = m->GetSlice("#kappa_{1}", "#kappa_{2}", m->GetBestFitParameters(), 200);

  double gmode[2];
  gmode[0] = m->GetBestFitParameter(0);
  gmode[1] = m->GetBestFitParameter(1);
  hist_all->SetGlobalMode(gmode);

  TCanvas* c1 = new TCanvas();
  c1->cd();
  hist_all->Draw("BTcB3CS1gmode");
  c1->SetGrid(1,1);
  c1->Print("kappa1_vs_kappa2.pdf");

  // clean up
  delete m;
  delete c1;

  // close log file
  BCLog::CloseLog();
   
  // no error
  return 0;
}

