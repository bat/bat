#include "GaussModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
GaussModel::GaussModel(const char* name) : BCModel(name)
{
    AddParameter("x", -10.0, 50.0); // index 0
    AddParameter("y",  -5.0,  5.0); // index 1
    AddParameter("z",  -5.0,  5.0); // index 2

    GetParameters().SetPriorConstantAll();
};

// ---------------------------------------------------------
GaussModel::~GaussModel()
{
};

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double>& parameters)
{
    double logprob = 0.;

    // Gaussian Likelihoods
    logprob += BCMath::LogGaus(parameters[0], 0.0, 2.0); // x
    logprob += BCMath::LogGaus(parameters[1], 0.0, 1.0); // y
    logprob += BCMath::LogGaus(parameters[2], 0.0, 1.0); // z

    return logprob;
}

