#include "BinomialModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
BinomialModel::BinomialModel(const char * name)
   : BCModel(name)
   , fNTotal(10)
   , fNSelected(5)
{
   // add a parameter which is the efficiency of observing a certain
   // amount of events starting from a larger set of events.
	 // and set it's prior flat
   AddParameter("epsilon", 0., 1., "#varepsilon"); // index 0
	 SetPriorConstant("epsilon");
}

// ---------------------------------------------------------
BinomialModel::~BinomialModel()
{
}

// ---------------------------------------------------------
double BinomialModel::LogLikelihood(const std::vector<double> &parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   double logprob = 0.;

   double eps = parameters.at(0);

   // calculate the binomial probability for observing a certain amount
   // of events given a larger number of events and an efficiency
   logprob += BCMath::LogApproxBinomial(fNTotal, fNSelected, eps);

   // return log likelihood
   return logprob;
}
