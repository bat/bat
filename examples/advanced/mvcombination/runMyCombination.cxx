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

int main()
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

	 // open input file
	 ifstream infile;
	 infile.open("input.txt", std::ifstream::in);

	 int nobservables;
	 int nmeasurements;
	 int nuncertainties;
	 
	 infile >> nobservables >> nmeasurements >> nuncertainties;

	 std::vector<std::string> observable_names;
	 
	 for (int i = 0; i < nobservables; ++i) {
		 std::string name;
		 double min;
		 double max;
		 infile >> name >> min >> max;
		 
		 // add observable
		 m->AddObservable(name.c_str(), min, max);
	 }

	 for (int i = 0; i < nuncertainties; ++i) {
		 std::string name;
		 infile >> name;
		 
		 // add uncertainty
		 m->AddUncertainty(name);
	 }

	 for (int i = 0; i < nmeasurements; ++i) {
		 std::string name;
		 std::string observable;
		 double central;
		 std::vector<double> uncertainties(0);
		 
		 infile >> name;
		 infile >> observable;
		 infile >> central;

		 for (int j = 0; j < nuncertainties; ++j) {
			 double uncertainty;
			 infile >> uncertainty;
			 uncertainties.push_back(uncertainty);
		 }

		 // add measurement
		 m->AddMeasurement(name, observable, central, uncertainties);
	 }

	 for (int i = 0; i < nuncertainties; ++i) {
		 TMatrixD mat(nmeasurements, nmeasurements);

		 for (int j = 0; j < nmeasurements; ++j) 
			 for (int k = 0; k < nmeasurements; ++k) {
				 double corr;
				 infile >> corr;
				 mat[j][k] = corr;
			 }

		 // set correlation matrix
		 m->GetUncertainty(i)->SetCorrelationMatrix(mat);
	 }

	 // close input file
	 infile.close();

	 // prepare analysis
	 for (int i = 0; i < nuncertainties; ++i)
		 m->CalculateCorrelationMatrix(i);

	 m->CalculateCovarianceMatrix();

	 // calculate BLUE
	 m->CalculateBLUE();

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

	 // clean up
   delete m;

   // close log file
   BCLog::CloseLog();

   return 0;

}

