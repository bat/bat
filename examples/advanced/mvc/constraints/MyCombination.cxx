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

	// Add observable "FR"
	AddObservable("FR",0.,0.4);

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
void MyCombination::CalculateObservables(const std::vector<double> & pars)
{
	// calculate FR
	GetObservable(0) -> Value( 1. - pars[0] - pars[1] );
}

// ---------------------------------------------------------
