// ***************************************************************
// This file was created using the ((PROGRAM)) script.
// ((PROGRAM)) is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "((MODEL)).h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
((MODEL))::((MODEL))(const std::string& name)
    : BCModel(name)
{
    // Define parameters here in the constructor. For example:
    // AddParameter("mu",-2,1,"#mu");
    // And set priors, if using built-in priors. For example:
    // GetParamater("mu").SetPrior(new BCPriorGaus(-1, 0.25));
}

// ---------------------------------------------------------
((MODEL))::~((MODEL))()
{
    // destructor
}

// ---------------------------------------------------------
double ((MODEL))::LogLikelihood(const std::vector<double>& parameters)
{
    // This returns the log of the conditional probability p(data|pars)
    // This is where you define your model.
    // BCMath contains many functions you will find helpful

    return -1;
}

// ---------------------------------------------------------
// double ((MODEL))::LogAPrioriProbability(const std::vector<double> & parameters) {
// This returns the log of the prior probability for the parameters
// If you use built-in 1D priors, don't uncomment this function.
// }
