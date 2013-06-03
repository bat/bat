//
// This ROOT macro is part of BAT and can only be run if BAT
// was installed correctly.
//
// The macro can be run from within ROOT via commands
//
//    root[1] .x singleChannel.C
//
// or
//
//    root[1] .L singleChannel.C
//    root[2] singleChannel()
//
// or from the command line
//
//    $ root singleChannel.C
//
// To improve the performance the macro can be run in a compiled
// mode. The commands are the same as above but with a '+' sign
// added to the name of the file, e.g.,
//
//    root[1] .x singleChannel.C+
//
// See ROOT documentation for details.
//
//
// Below are the includes needed for compilation of the macro
// the #if ... #endif directives around the includes allow to
// run the macro in both normal and compiled mode.
#if !defined(__CINT__) || defined(__MAKECINT__)

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <BAT/BCMVCombination.h>

#include <iostream>
#include <fstream>

#endif

int main(int argc, char *argv[])
{

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new BCMVCombination object
  BCMVCombination * m = new BCMVCombination();

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

  // calculate BLUE
  m->CalculateBLUE();
   
  // print BLUE results to file
  m->PrintBLUEResults("BCMVCombination_BLUE.txt");

  // clean up
  delete m;
   
  // close log file
  BCLog::CloseLog();
   
  // no error
  return 0;
}

