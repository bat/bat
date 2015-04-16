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
{
    AddObservable("rho", 0, 4.5, "#rho");
    GetObservable("rho").SetNbins(500);
}

// ---------------------------------------------------------
MyCombination::~MyCombination()
{
}

// ---------------------------------------------------------
double MyCombination::LogLikelihood(const std::vector<double>& parameters)
{
    // check physical constraints
    if (fFlagPhysicalConstraints && (parameters[1] < 0))
        return -1e55;

    return BCMVCombination::LogLikelihood(parameters);
}

// ---------------------------------------------------------
void MyCombination::CalculateObservables(const std::vector<double>& pars)
{
    fObservables[0].Value(pars[0]*pars[1]);
}
