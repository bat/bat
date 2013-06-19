#include "BinomialModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
BinomialModel::BinomialModel()
   : BCModel()
   , fNTotal(10)
   , fNSelected(5)
{
  DefineParameters();
}

// ---------------------------------------------------------
BinomialModel::BinomialModel(const char * name)
   : BCModel(name)
{
  DefineParameters();
}

// ---------------------------------------------------------
BinomialModel::~BinomialModel()
{
}

// ---------------------------------------------------------
void BinomialModel::DefineParameters()
{
   // add a parameter which is the efficiency of observing a certain
   // amount of events starting from a larger set of events.
   AddParameter("epsilon", 0., 1., "#varepsilon"); // index 0
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

// ---------------------------------------------------------
double BinomialModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;

   // get width of the parameter range
   double deps = GetParameter(0)->GetRangeWidth();

   // add a flat prior probability
   logprob += log(1./deps); // flat prior

   // return log prior
   return logprob;
}
// ---------------------------------------------------------

