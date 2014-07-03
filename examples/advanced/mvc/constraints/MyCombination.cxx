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
double MyCombination::FR(const std::vector<double> & pars)
{
	return 1. - pars[0] - pars[1];
}

// ---------------------------------------------------------
