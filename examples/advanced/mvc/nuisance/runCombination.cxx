#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH2D.h>

#include <BAT/BCMVCombination.h>

#include <iostream>
#include <fstream>

int main(int argc, char *argv[])
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new BCMVCombination object
  BCMVCombination * m = new BCMVCombination();

  // set Metropolis as marginalization method
  m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

  // set precision
  m->MCMCSetPrecision(BCEngineMCMC::kMedium);
   
	// read input from file
	int isopen = m->ReadInput(argv[1]);
	if (!isopen) {
		std::cout << "Could not open file. Exit." << std::endl;
		return 1;
	}

  // perform numerical analysis using MCMC
  m->MarginalizeAll();
   
  // find mode using Minuit
  m->FindMode( m->GetBestFitParameters() );

	// print all marginalized distributions   
  m->PrintAllMarginalized("BCMVCombination_plots.pdf");
   
  // print results of numerical analysis
  m->PrintResults("BCMVCombination_results.txt");

	// print a single 2D plot
	m->GetSlice("obs", "rho_1", m->GetBestFitParameters(), 200)->Print("rho1_vs_obs.pdf", "BTcB3CS1gmodeprofiley");

  // clean up
  delete m;
   
  // close log file
  BCLog::CloseLog();
   
  // no error
  return 0;
}

