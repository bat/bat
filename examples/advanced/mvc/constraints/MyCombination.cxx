#include <iomanip>

#include <TMath.h>
#include <TH1D.h>

#include "MyCombination.h"
#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MyCombination::MyCombination() 
  : BCMVCombination()
  , fFlagPhysicalConstraints(true)
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
  double logprob =   BCMVCombination::LogLikelihood(parameters);

  double F0 = parameters[0];
  double FL = parameters[1];
  
  // check physical constraints
  if (fFlagPhysicalConstraints && (F0+FL > 1))
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
    double F0 = fMCMCx.at(i * npar + 0);
    double FL = fMCMCx.at(i * npar + 1);
    double FR = 1. - F0 - FL;

    // fill the ratio histogram
    if (fHistFR)
      fHistFR->Fill(FR);
  }
}

// ---------------------------------------------------------
