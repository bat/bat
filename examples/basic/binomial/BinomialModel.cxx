#include "BinomialModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
BinomialModel::BinomialModel(const std::string& name, unsigned ntotal, unsigned nselected)
    : BCModel(name)
    , fNTotal(ntotal)
    , fNSelected(nselected)
{
    // add a parameter which is the efficiency of observing a certain
    // amount of events starting from a larger set of events.
    // and set it's prior flat
    AddParameter("epsilon", 0., 1., "#varepsilon");
    SetPriorConstant("epsilon");
}

// ---------------------------------------------------------
BinomialModel::~BinomialModel()
{
}

// ---------------------------------------------------------
double BinomialModel::LogLikelihood(const std::vector<double>& parameters)
{
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.

    // calculate the binomial probability for observing a certain amount
    // of events given a larger number of events and an efficiency
    return BCMath::LogApproxBinomial(fNTotal, fNSelected, parameters[0]);
}
