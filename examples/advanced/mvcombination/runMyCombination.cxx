// ***************************************************************
// This file was created using the CreateProject.sh script
// for project MVCombination.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include <TMatrixT.h>

#include <iostream>
#include <fstream>

#include "MyCombination.h"
#include "ToyModel.h"

int main(int argc, char *argv[])
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new MyCombination object
   MyCombination * m = new MyCombination();

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

   // define histogram for FR
   TH1D* hist_fr = new TH1D("FR", ";FR;p", 100, 0., 0.4);
   hist_fr->SetStats(kFALSE);
   
   m->SetHistFR(hist_fr);
   
   // perform numerical analysis using MCMC
   m->MarginalizeAll();
   
   // find mode using Minuit
   m->FindMode( m->GetBestFitParameters() );
   
   m->PrintAllMarginalized("MyCombination_plots.pdf");
   
   BCH1D* histFR = new BCH1D(hist_fr);
   
   histFR->Print("FR.pdf", "BTulB3L");
   
   // print results of numerical analysis
   m->PrintResults("MyCombination_results.txt");
   
   // print summary to screen
   m->PrintSummary();

    // test goodness-of-fit

   ToyModel* toy = new ToyModel();
   
   TH1D* hist_chi2 = new TH1D("chi2", ";#chi^{2};p(#chi^{2})", 100, 0., 30.);
   hist_chi2->SetStats(kFALSE);
   
   toy->SetHistChi2(hist_chi2);
   toy->SetNMeasurements(m->GetNMeasurements(), 0., 1.3);
   toy->SetVectorMeasurements(m->GetVectorMeasurements());
   toy->SetVectorObservable(m->GetVectorObservable());
   toy->SetCovarianceMatrix(m->GetCovarianceMatrix());
   toy->SetParameters(m->GetBestFitParameters());
   
   toy->MarginalizeAll();
   toy->GetBestFitParameters();
   toy->PrintAllMarginalized("toys.pdf");
   toy->PrintResults("toys.txt");
   toy->PrintToys("chi2.pdf");

   // clean up
   delete m;
   delete toy; 

   // close log file
   BCLog::CloseLog();

   // no error
   return 0;

}

