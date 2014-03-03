#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <TMatrixT.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>
#include <fstream>

#include "MyCombination.h"
#include <BAT/BCMVCombination.h>
#include <BAT/BCMVCDataModel.h>
#include <BAT/BCMVCMeasurement.h>
#include <BAT/BCMVCUncertainty.h>

int main(int argc, char *argv[])
{

  // flags 
  bool flag_phys = true;   // physical constrained on (true) or off (false)
  bool flag_full = true;   // run full analysis
  bool flag_meas = true;   // repeat analysis and remove one measurement at a time
  bool flag_unc  = true;   // repeat analysis and remove one uncertainty at a time
  bool flag_gof  = true;   // perform goodness-of-fit test

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // helper variables
  double rho_all;
  double area_all;

  std::vector<double> rho_single;
  std::vector<double> area_single;
  std::vector<double> rho_unc;
  std::vector<double> area_unc;

  // create new MyCombination object
  MyCombination * m = new MyCombination();

  // set Metropolis as marginalization method
  m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

  // set flag for physical constraints
  m->SetFlagPhysicalConstraints(flag_phys);

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
   
  // get number of measurements
  int nmeas = m->GetNMeasurements();

  // get number of uncertainties
  int nunc  = m->GetNUncertainties();

  // ----- Full 2D analysis --------------------

  if (flag_full) {
    // perform numerical analysis using MCMC
    m->MarginalizeAll();

    // reset histogram pointer
    m->SetHistFR(0);
   
    // find mode using Minuit
    m->FindMode( m->GetBestFitParameters() );
   
    m->PrintAllMarginalized("MyCombination_full_plots.pdf");

    double gmode[2];
    gmode[0] = m->GetBestFitParameters().at(0);
    gmode[1] = m->GetBestFitParameters().at(1);
    
    BCH2D* hist_slice = m->GetSlice("F0", "FL", m->GetBestFitParameters(), 300);
    hist_slice->GetHistogram()->SetStats(kFALSE);
    hist_slice->SetGlobalMode(gmode);
    hist_slice->Print("F0_vs_FL.pdf", "BTfB3CS1gmode");

    BCH1D* hist_F0 = new BCH1D( hist_slice->GetHistogram()->ProjectionX() );
    hist_F0->GetHistogram()->SetStats(kFALSE);
    hist_F0->SetGlobalMode(gmode[0]);
    hist_F0->Print("F0.pdf", "BTciB3CS1D0pdf0Lmode");

    BCH1D* hist_FL = new BCH1D( hist_slice->GetHistogram()->ProjectionY() );
    hist_FL->GetHistogram()->SetStats(kFALSE);
    hist_FL->SetGlobalMode(gmode[1]);
    hist_FL->Print("FL.pdf", "BTciB3CS1D0pdf0Lmode");
    
    BCH1D* hist_FR = new BCH1D(hist_fr);   
    hist_FR->Print("FR.pdf", "BTulB3L");

    // calculate correlation between F0 and FL
    rho_all = hist_slice->GetHistogram()->GetCorrelationFactor();

    // calculate size of uncertainty
    area_all = hist_slice->GetArea(0.39);

    // print results of numerical analysis to file
    m->PrintResults("MyCombination_full_results.txt");

    // calculate BLUE
    m->CalculateBLUE();
   
    // print BLUE results to file
    m->PrintBLUEResults("MyCombination_full_BLUE.txt");
  }

   // ----- Remove one measurement at a time --------------------

  // test single measurements
  if (flag_meas) {

    // switch off ith measurement
    for (int i = 0; i < nmeas; ++i) {
      for (int j = 0; j < nmeas; ++j) {
	if (i == j)
	  m->GetMeasurement(j)->SetFlagActive(false);
	else 
	  m->GetMeasurement(j)->SetFlagActive(true);
      }
      m->PrepareAnalysis();
      m->MarginalizeAll();
      m->FindMode( m->GetBestFitParameters() );
      m->PrintAllMarginalized(Form("MyCombination_measurement_%i_plots.pdf", i));
      m->PrintResults(Form("MyCombination_measurement_%i_results.txt", i));
      rho_single.push_back( m->GetMarginalized("F0", "FL")->GetHistogram()->GetCorrelationFactor());
      area_single.push_back(m->GetMarginalized("F0", "FL")->GetArea(0.39));
    }
    
    // switch all measurements on
    for (int j = 0; j < nmeas; ++j) 
      m->GetMeasurement(j)->SetFlagActive(true);
    m->PrepareAnalysis();
  }
  
   // ----- Remove one uncertainty at a time --------------------

  // test single uncertainties
  if (flag_unc) {

    // switch off ith uncertainty
    for (int i = 0; i < nunc; ++i) {
      for (int j = 0; j < nunc; ++j) {
	if (i == j)
	  m->GetUncertainty(j)->SetFlagActive(false);
	else 
	  m->GetUncertainty(j)->SetFlagActive(true);
      }
      m->PrepareAnalysis();
      m->MarginalizeAll();
      m->FindMode( m->GetBestFitParameters() );
      m->PrintAllMarginalized(Form("MyCombination_uncertainty_%i_plots.pdf", i));
      m->PrintResults(Form("MyCombination_uncertainty_%i_results.txt", i));
      rho_unc.push_back( m->GetMarginalized("F0", "FL")->GetHistogram()->GetCorrelationFactor());
      area_unc.push_back(m->GetMarginalized("F0", "FL")->GetArea(0.39));
    }

    // switch all uncertainties on
    for (int j = 0; j < nunc; ++j) 
      m->GetUncertainty(j)->SetFlagActive(true);
    m->PrepareAnalysis();
  }
  
  // ----- Goodness-of-fit --------------------
  
  if (flag_gof) {
    // test goodness-of-fit
    BCMVCDataModel* dm = new BCMVCDataModel(m);
   
    TH1D* hist_chi2 = new TH1D("chi2", ";#chi^{2};p(#chi^{2})", 100, 0., 30.);
    hist_chi2->SetStats(kFALSE);
   
    dm->SetHistChi2(hist_chi2);
    dm->SetMeasurementRanges(-0.5, 2.0); 
    std::vector<double> SM(2);
    SM[0]=0.687;
    SM[1]=0.311;
    dm->SetParameters(SM);
    
    dm->MarginalizeAll();
    dm->GetBestFitParameters();
    dm->PrintAllMarginalized("gof_plots.pdf");
    dm->PrintResults("gof_results.txt");
    dm->PrintToys("gof_chi2.pdf");
    dm->PrintSummary();

    // free memory
    delete dm; 
  }

  // ----- Print results to screen --------------------

  // print to screen
  if (flag_full) {
    std::cout << " rho (all)   : " << rho_all << std::endl;
    std::cout << " Area  (all) : " << area_all << std::endl;
    std::cout << std::endl;
  }
  
  if (flag_meas) {
    std::cout << " Measurements : " << std::endl;
    for (int i = 0; i < nmeas; ++i) {
      std::cout << " Area / rho  (" << i << ")   : " << (area_single[i] - area_all) / area_all * 100 << "%";
      std::cout << " / " << rho_single[i] << std::endl;
    }
    std::cout << std::endl;
  }
  
  if (flag_unc) {
    std::cout << " Uncertainties : " << std::endl;
    for (int i = 0; i < nunc; ++i) {
      std::cout << " Area / rho  (" << i << ")   : " << (area_unc[i] - area_all) / area_all * 100 << "%";
      std::cout << " / " << rho_unc[i] << std::endl;
    }
    std::cout << std::endl;
  }
  
  // free memory
  delete m;

  // close log file
  BCLog::CloseLog();

  // no error
  return 0;

}

