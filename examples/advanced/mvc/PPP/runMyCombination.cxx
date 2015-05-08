#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <TMatrixT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>

#include "MyCombination.h"
#include <BAT/BCMVCombination.h>
#include <BAT/BCMVCDataModel.h>
#include <BAT/BCMVCMeasurement.h>
#include <BAT/BCMVCUncertainty.h>

int main(int argc, char* argv[])
{

    // flags
    bool flag_phys = true;   // physical constrained on (true) or off (false)
    bool flag_full = true;   // run full analysis
    bool flag_meas = true;  // repeat analysis and remove one measurement at a time
    bool flag_unc  = true;  // repeat analysis and remove one uncertainty at a time

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // helper variables
    double rho_mean_all = 0;
    double rho_std_all = 0;

    std::vector<double> rho_mean_single;
    std::vector<double> rho_std_single;
    std::vector<double> rho_mean_unc;
    std::vector<double> rho_std_unc;

    // create new MyCombination object
    MyCombination* m = new MyCombination();

    // set flag for physical constraints
    m->SetFlagPhysicalConstraints(flag_phys);

    // set Metropolis as marginalization method
    m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

    // set precision
    m->MCMCSetPrecision(BCIntegrate::kMedium);

    if (argc == 2) {
        int isopen = m->ReadInput(argv[1]);
        if (!isopen) {
            std::cout << "Could not open file. Exit." << std::endl;
            return 1;
        }
    } else {
        std::cout << "No input file specified. Exit." << std::endl;
        return 1;
    }

    // get number of measurements
    int nmeas = m->GetNMeasurements();

    // get number of uncertainties
    int nunc  = m->GetNUncertainties();

    // ----- Full 2D analysis --------------------

    if (flag_full) {

        // perform numerical analysis using MCMC
        m->MarginalizeAll();

        // find mode using Minuit
        m->FindMode(m->GetGlobalMode());

        m->PrintAllMarginalized("MyCombination_full_plots.pdf");

        // print results of numerical analysis to the log
        m->PrintSummary();

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
            m->FindMode(m->GetGlobalMode());
            m->PrintAllMarginalized(Form("MyCombination_measurement_%i_plots.pdf", i));
            BCLog::OutSummary(Form("Switching off %d-th measurement",i));
            m->PrintSummary();

            rho_mean_single.push_back(m->MCMCGetStatistics().mean[2]);
            rho_std_single.push_back(sqrt(m->MCMCGetStatistics().variance[2]));
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
            m->FindMode(m->GetGlobalMode());
            m->PrintAllMarginalized(Form("MyCombination_uncertainty_%i_plots.pdf", i));
            BCLog::OutSummary(Form("Switching off %d-th uncertainty",i));
            m->PrintSummary();

            rho_mean_unc.push_back(m->MCMCGetStatistics().mean[2]);
            rho_std_unc.push_back(sqrt(m->MCMCGetStatistics().variance[2]));
        }

        // switch all uncertainties on
        for (int j = 0; j < nunc; ++j)
            m->GetUncertainty(j)->SetFlagActive(true);
        m->PrepareAnalysis();

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

