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
  bool flag_meas = true;  // repeat analysis and remove one measurement at a time
  bool flag_unc  = true;  // repeat analysis and remove one uncertainty at a time

  int nbins = 500;
  double rho_min   = 0.0; 
  double rho_max   = 4.5;
  double alpha_min = 1.5;
  double alpha_max = 3.0;
  double eta_min   = 0.0;
  double eta_max   = 4.0;

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // helper variables
  double rho_mean_all;
  double rho_std_all;

  std::vector<double> rho_mean_single;
  std::vector<double> rho_std_single;
  std::vector<double> rho_mean_unc;
  std::vector<double> rho_std_unc;

  // create new MyCombination object
  MyCombination * m = new MyCombination();

  // set flag for physical constraints
  m->SetFlagPhysicalConstraints(flag_phys);

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

  // get number of measurements
  int nmeas = m->GetNMeasurements();

  // get number of uncertainties
  int nunc  = m->GetNUncertainties();

  // ----- Full 2D analysis --------------------

  if (flag_full) {
    // define histogram for rho
    TH1D* hist_rho = new TH1D("hist_rho", ";#rho;p(#rho)", nbins, rho_min, rho_max);
    hist_rho->SetStats(kFALSE);
    TH2D* hist_rhoalpha = new TH2D("hist_rhoalpha", ";#rho;#alpha", nbins, rho_min, rho_max, nbins, alpha_min, alpha_max);
    hist_rhoalpha->SetStats(kFALSE);
    TH2D* hist_rhoeta = new TH2D("hist_rhoeta", ";#rho;#eta", nbins, rho_min, rho_max, nbins, eta_min, eta_max);
    hist_rhoeta->SetStats(kFALSE);
    
    m->SetHistRho(hist_rho);
    m->SetHistRhoAlpha(hist_rhoalpha);
    m->SetHistRhoEta(hist_rhoeta);
    
    // perform numerical analysis using MCMC
    m->MarginalizeAll();
    
    // find mode using Minuit
    m->FindMode( m->GetBestFitParameters() );
   
    m->PrintAllMarginalized("MyCombination_full_plots.pdf");

    double gmode[2];
    gmode[0] = m->GetBestFitParameters().at(0);
    gmode[1] = m->GetBestFitParameters().at(1);
    
    BCH2D* hist_slice = m->GetSlice("alpha", "eta", m->GetBestFitParameters(), 300);
    hist_slice->GetHistogram()->SetStats(kFALSE);
    hist_slice->SetGlobalMode(gmode);
    hist_slice->Print("alpha_vs_eta.pdf", "BTfB3CS1gmode");

    BCH1D* bchist_rho = new BCH1D(hist_rho);   
    bchist_rho->Print("rho.pdf", "BTfB3Lmodemeanrms");

    BCH2D* bchist_rhoalpha = new BCH2D(hist_rhoalpha);   
    bchist_rhoalpha->Print("rho_vs_alpha.pdf", "BTfB3Lmodemeanrms");

    BCH2D* bchist_rhoeta = new BCH2D(hist_rhoeta);   
    bchist_rhoeta->Print("rho_vs_eta.pdf", "BTfB3Lmodemeanrms");

    // calculate correlation between alpha and eta
    rho_mean_all = bchist_rho->GetMean();
    rho_std_all  = bchist_rho->GetSTD(); 

       // print results to screen
    std::cout << " E[rho]   = " << bchist_rho->GetMean() << std::endl;
    std::cout << " Std[rho] = " << bchist_rho->GetSTD() << std::endl;
    std::cout << " V[rho]   = " << bchist_rho->GetVariance() << std::endl;

    // print results of numerical analysis to file
    m->PrintResults("MyCombination_full_results.txt");

    // calculate BLUE
    m->CalculateBLUE();
   
    // print BLUE results to file
    m->PrintBLUEResults("MyCombination_full_BLUE.txt");

    // reset histogram pointer
    m->SetHistRho(0);
    m->SetHistRhoAlpha(0);
    m->SetHistRhoEta(0);
  }

   // ----- Remove one measurement at a time --------------------

  // test single measurements
  if (flag_meas) {

    // define histogram for rho
    TH1D* hist_rho = new TH1D("hist_rho", ";#rho;p(#rho)", nbins, rho_min, rho_max);
    hist_rho->SetStats(kFALSE);
    TH2D* hist_rhoalpha = new TH2D("hist_rhoalpha", ";#rho;#alpha", nbins, rho_min, rho_max, nbins, alpha_min, alpha_max);
    hist_rhoalpha->SetStats(kFALSE);
    TH2D* hist_rhoeta = new TH2D("hist_rhoeta", ";#rho;#eta", nbins, rho_min, rho_max, nbins, eta_min, eta_max);
    hist_rhoeta->SetStats(kFALSE);
    
    m->SetHistRho(hist_rho);
    m->SetHistRhoAlpha(hist_rhoalpha);
    m->SetHistRhoEta(hist_rhoeta);

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

      BCH1D* bchist_rho = new BCH1D(hist_rho);   
      bchist_rho->Print(Form("rho_measurement_%i.pdf", i), "BTfB3Lmodemeanrms");

      BCH2D* bchist_rhoalpha = new BCH2D(hist_rhoalpha);   
      bchist_rhoalpha->Print(Form("rho_vs_alpha_measurement_%i.pdf", i), "BTfB3Lmodemeanrms");

      BCH2D* bchist_rhoeta = new BCH2D(hist_rhoeta);   
      bchist_rhoeta->Print(Form("rho_vs_eta_measurement_%i.pdf", i), "BTfB3Lmodemeanrms");

      rho_mean_single.push_back(bchist_rho->GetMean());
      rho_std_single.push_back(bchist_rho->GetSTD()); 
    }
    
    // switch all measurements on
    for (int j = 0; j < nmeas; ++j) 
      m->GetMeasurement(j)->SetFlagActive(true);
    m->PrepareAnalysis();

    // reset histogram pointer
    m->SetHistRho(0);
    m->SetHistRhoAlpha(0);
    m->SetHistRhoEta(0);
  }
  
   // ----- Remove one uncertainty at a time --------------------

  // test single uncertainties
  if (flag_unc) {

    // define histogram for rho
    TH1D* hist_rho = new TH1D("hist_rho", ";#rho;p(#rho)", nbins, rho_min, rho_max);
    hist_rho->SetStats(kFALSE);
    TH2D* hist_rhoalpha = new TH2D("hist_rhoalpha", ";#rho;#alpha", nbins, rho_min, rho_max, nbins, alpha_min, alpha_max);
    hist_rhoalpha->SetStats(kFALSE);
    TH2D* hist_rhoeta = new TH2D("hist_rhoeta", ";#rho;#eta", nbins, rho_min, rho_max, nbins, eta_min, eta_max);
    hist_rhoeta->SetStats(kFALSE);

    m->SetHistRho(hist_rho);
    m->SetHistRhoAlpha(hist_rhoalpha);
    m->SetHistRhoEta(hist_rhoeta);

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

      BCH1D* bchist_rho = new BCH1D(hist_rho);   
      bchist_rho->Print(Form("rho_uncertainty_%i.pdf", i), "BTfB3Lmodemeanrms");

      BCH2D* bchist_rhoalpha = new BCH2D(hist_rhoalpha);   
      bchist_rhoalpha->Print(Form("rho_vs_alpha_uncertainty_%i.pdf", i), "BTfB3Lmodemeanrms");

      BCH2D* bchist_rhoeta = new BCH2D(hist_rhoeta);   
      bchist_rhoeta->Print(Form("rho_vs_eta_uncertainty_%i.pdf", i), "BTfB3Lmodemeanrms");

     rho_mean_unc.push_back(bchist_rho->GetMean());
      rho_std_unc.push_back(bchist_rho->GetSTD()); 
    }

    // switch all uncertainties on
    for (int j = 0; j < nunc; ++j) 
      m->GetUncertainty(j)->SetFlagActive(true);
    m->PrepareAnalysis();

    // reset histogram pointer
    m->SetHistRho(0);
    m->SetHistRhoAlpha(0);
    m->SetHistRhoEta(0);
  }
  
  // ----- Print results to screen --------------------

  // print to screen
  if (flag_full) {
    std::cout << " rho mean (all)   : " << rho_mean_all << std::endl;
    std::cout << " rho std (all)    : " << rho_std_all << std::endl;
    std::cout << std::endl;
  }
  
  if (flag_meas) {
    std::cout << " Measurements : " << std::endl;
    for (int i = 0; i < nmeas; ++i) {
      std::cout << " change in std  (" << i << ")   : " << (rho_std_single[i] - rho_std_all) / rho_std_all * 100 << "%" << std::endl;
    }
    std::cout << std::endl;
  }
  
  if (flag_unc) {
    std::cout << " Uncertainties : " << std::endl;
    for (int i = 0; i < nunc; ++i) {
      std::cout << " change in std  (" << i << ")   : " << (rho_std_unc[i] - rho_std_all) / rho_std_all * 100 << "%" << std::endl;
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

