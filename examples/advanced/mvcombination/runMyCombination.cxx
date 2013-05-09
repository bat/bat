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
#include <TH1D.h>
#include <TH2D.h>

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

   double gmode[2];
   gmode[0] = m->GetBestFitParameters().at(0);
   gmode[1] = m->GetBestFitParameters().at(1);

   //   BCH2D* hist_slice = m->GetSlice("F0", "FL", m->GetBestFitParameters(), 1000);
   BCH2D* hist_slice=m->GetMarginalized("F0", "FL");
   hist_slice->GetHistogram()->SetStats(kFALSE);
   hist_slice->SetGlobalMode(gmode);
   hist_slice->Print("F0_vs_FL.pdf", "BTfB3CS1gmode");
   //   BCH1D* hist_F0 = new BCH1D( hist_slice->GetHistogram()->ProjectionX() );
   BCH1D* hist_F0 = m->GetMarginalized("F0");
   hist_F0->GetHistogram()->SetStats(kFALSE);
   hist_F0->SetGlobalMode(gmode[0]);
   hist_F0->Print("F0.pdf", "BTciB3CS1D0pdf0Lmode");
   //   BCH1D* hist_FL = new BCH1D( hist_slice->GetHistogram()->ProjectionY() );
   BCH1D* hist_FL = m->GetMarginalized("FL");
   hist_FL->GetHistogram()->SetStats(kFALSE);
   hist_FL->SetGlobalMode(gmode[1]);
   hist_FL->Print("FL.pdf", "BTciB3CS1D0pdf0Lmode");

   BCH1D* hist_FR = new BCH1D(hist_fr);   
   hist_FR->Print("FR.pdf", "BTulB3L");

   // print results of numerical analysis to file
   m->PrintResults("MyCombination_results.txt");

	 // calculate BLUE
	 m->CalculateBLUE();
   
   // print BLUE results to file
   m->PrintBLUEResults("MyCombination_BLUE.txt");

   // test goodness-of-fit
   ToyModel* toy = new ToyModel(m);
   
   TH1D* hist_chi2 = new TH1D("chi2", ";#chi^{2};p(#chi^{2})", 100, 0., 30.);
   hist_chi2->SetStats(kFALSE);
   
   toy->SetHistChi2(hist_chi2);
   toy->SetMeasurementRanges(-0.5, 2.0); 
   std::vector<double> SM(2);
   SM[0]=0.687;
   SM[1]=0.311;
   toy->SetParameters(SM);
   
   toy->MarginalizeAll();
   toy->GetBestFitParameters();
   toy->PrintAllMarginalized("toys.pdf");
   toy->PrintResults("toys.txt");
   toy->PrintToys("chi2.pdf");
   toy->PrintSummary();

   // clean up
   delete m;
   delete toy; 

   // close log file
   BCLog::CloseLog();

   // no error
   return 0;

}

