#include "MyCombination.h"

#include <iomanip>

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MyCombination::MyCombination() : MVCombination()
															 , fHistFR(0)
{
}

// ---------------------------------------------------------
MyCombination::~MyCombination()
{
}

// ---------------------------------------------------------
double MyCombination::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob =   MVCombination::LogLikelihood(parameters);

  double F0 = parameters[0];
  double FL = parameters[1];
  
  if (F0+FL > 1)
    return -1e55;
  
  return logprob;
}

// ---------------------------------------------------------
void MyCombination::MCMCIterationInterface()
{
  // get number of chains
  int nchains = MCMCGetNChains();

  // get number of parameters
  int npar = GetNParameters();

  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    // get the current values of the parameters x and y. These are
    // stored in fMCMCx.
    double F0 = fMCMCx.at(i * npar + 0);
    double FL = fMCMCx.at(i * npar + 1);
		double FR = 1. - F0 - FL;

    // fill the ratio histogram
    fHistFR->Fill(FR);
  }
}

// ---------------------------------------------------------
