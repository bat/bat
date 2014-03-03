#include "PoissonModel.h"

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
PoissonModel::PoissonModel(const char * name) :
BCModel(name),
fNObs(0)
{
   DefineParameters();
};

// ---------------------------------------------------------
PoissonModel::~PoissonModel()
{
};

// ---------------------------------------------------------
void PoissonModel::DefineParameters()
{
   // add a parameter for the number of expected events. The range will
   // be adjusted later according to the number of observed events.
   AddParameter("#lambda", 0., 7.); // index 0
}

// ---------------------------------------------------------
void PoissonModel::SetNObs(unsigned nobs)
{
   // set number of observed events
   fNObs = nobs;

   // adjust parameter ranges
   double lambdamin = 0;
   double lambdamax = 10;

   // the adjustment depends on the number of observed events
   if (nobs >= 5 && nobs < 10) {
      lambdamin = 0;
      lambdamax = 20;
   }
   else if (nobs >= 10 && nobs < 20) {
      lambdamin = 0;
      lambdamax = 40;
   }
   else if (nobs >= 20 && nobs < 30) {
      lambdamin = 5;
      lambdamax = 55;
   }
   else if (nobs >= 30) {
      lambdamin = double(nobs) - 5 * sqrt(double(nobs));
      lambdamax = double(nobs) + 5 * sqrt(double(nobs));
   }

   // re-set the parameter range
   GetParameter(0)->SetLimits(lambdamin, lambdamax);
}

// ---------------------------------------------------------
double PoissonModel::LogLikelihood(const std::vector<double> & parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   double logprob = 0.;

   double lambda = parameters.at(0);

   // log Poisson term
   logprob += BCMath::LogPoisson(fNObs, lambda);

   // return log likelihood
   return logprob;
}

// ---------------------------------------------------------
double PoissonModel::LogAPrioriProbability(const std::vector<double> & parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;

   // add a flat prior probability
   logprob += log(1. / GetParameter(0)->GetRangeWidth()); // flat prior

   // return log prior
   return logprob;
}
