#include "GaussModel.h"

#include <BAT/BCMath.h>

#include <TRandom3.h>

// ---------------------------------------------------------
GaussModel::GaussModel(const std::string& name) : BCModel(name)
{
    // add parameters x and y
    AddParameter("x", -10.0, 50.0); // index 0
    AddParameter("y",  -5.0,  5.0); // index 1

    SetPriorConstantAll();
}

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double>& parameters)
{
    double logprob = 0.;

    // Likelihood = Gaus(x | mean = 0, std.dev. = 2) * Gaus(y | mean = 0, std.dev. = 1)
    logprob += BCMath::LogGaus(parameters[0], 0.0, 2.0);
    logprob += BCMath::LogGaus(parameters[1], 0.0, 1.0);

    return logprob;
}

// --------------------------------------------------------
double GaussModel::MCMCTrialFunctionSingle(unsigned ichain, unsigned ipar)
{
    // no check of range for performance reasons

    // get scale factor from an array of scale factors. the size of the
    // array is number of chains times number of parameters.
    // double scale = fMCMCTrialFunctionScaleFactor[ichain * GetNParameters() + ipar];

    // choose trial function by uncommenting choices below

    // Gaussian with fixed width
    //  return fRandom.Gaus(0.0, 1.0);

    // Gaussian with adjustable width
    //  return fRandom.Gaus(0.0, scale);

    // Breit-Wigner with adjustable width
    //  return fRandom.BreitWigner(0.0, scale);

    // Flat function with adjustable width
    //  return scale * 2. * (0.5 - fRandom.Uniform());

    // Flat function with fixed width
    return 0.02 * (0.5 - fRandom.Uniform());
}


