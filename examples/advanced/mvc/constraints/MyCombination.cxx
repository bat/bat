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
    AddObservable("FR", 0., 0.4, "F_{R}");
}

// ---------------------------------------------------------
MyCombination::~MyCombination()
{
}

// ---------------------------------------------------------
double MyCombination::LogLikelihood(const std::vector<double>& parameters)
{
    // check physical constraints
    if (fFlagPhysicalConstraints && (parameters[0] + parameters[1] > 1))
        return -1e55;

    return BCMVCombination::LogLikelihood(parameters);
}

// ---------------------------------------------------------
void MyCombination::CalculateObservables(const std::vector<double>& pars)
{
    // calculate FR
    GetObservable(0).Value( 1. - pars[0] - pars[1] );
}

// ---------------------------------------------------------
