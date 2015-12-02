#include "PoissonModel.h"

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
PoissonModel::PoissonModel(const std::string& name)
    :	BCModel(name),
      fNObs(0)
{
    // add a parameter for the number of expected events. The range will
    // be adjusted later according to the number of observed events.
    AddParameter("lambda", 0., 7., "#lambda"); // index 0
    GetParameter("lambda").SetPriorConstant();
}


// ---------------------------------------------------------
void PoissonModel::SetNObs(unsigned nobs)
{
    // set number of observed events
    fNObs = nobs;

    // take [0,2*nobs] if nobs small
    if (fNObs < 30)
        GetParameter(0).SetLimits(0, 2 * fNObs);
    // take 5*sigma window under Gaussian equivalent if nobs large
    else
        GetParameter(0).SetLimits(fNObs - 5 * sqrt(fNObs), nobs + 5 * sqrt(fNObs));

}

// ---------------------------------------------------------
double PoissonModel::LogLikelihood(const std::vector<double>& parameters)
{
    // Likelihood = Poisson(fNobs | expectation = lambda)
    return BCMath::LogPoisson(fNObs, parameters[0]);
}
