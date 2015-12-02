#include "GaussModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
GaussModel::GaussModel(const std::string& name) : BCModel(name)
{
    AddParameter("x", -10.0, 50.0); // index 0
    AddParameter("y",  -5.0,  5.0); // index 1
    AddParameter("z",  -5.0,  5.0); // index 2

    GetParameters().SetPriorConstantAll();
}

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double>& parameters)
{
    double logprob = 0.;

    // Likelihood = Gaus(x | mean = 0, std. dev. = 2)
    //            * Gaus(y | mean = 0, std. dev. = 1)
    //            * Gaus(z | mean = 0, std. dev. = 1)
    logprob += BCMath::LogGaus(parameters[0], 0.0, 2.0); // x
    logprob += BCMath::LogGaus(parameters[1], 0.0, 1.0); // y
    logprob += BCMath::LogGaus(parameters[2], 0.0, 1.0); // z

    return logprob;
}
