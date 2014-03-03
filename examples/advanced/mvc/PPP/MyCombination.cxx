#include <iomanip>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>

#include "MyCombination.h"
#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MyCombination::MyCombination() 
  : BCMVCombination()
  , fFlagPhysicalConstraints(true)
  , fHistRho(0)
  , fHistRhoAlpha(0)
  , fHistRhoEta(0)
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

  //  double alpha = parameters[0];
  double eta   = parameters[1];
  
  // check physical constraints
  if (fFlagPhysicalConstraints && (eta < 0))
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
    double alpha = fMCMCx.at(i * npar + 0);
    double eta   = fMCMCx.at(i * npar + 1);
    double rho = alpha * eta;

    // fill the ratio histogram
    if (fHistRho)
      fHistRho->Fill(rho);
    if (fHistRhoAlpha)
      fHistRhoAlpha->Fill(rho, alpha);
    if (fHistRhoEta)
      fHistRhoEta->Fill(rho,eta);
  }
}

// ---------------------------------------------------------
