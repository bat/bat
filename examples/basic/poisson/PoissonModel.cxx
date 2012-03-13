#include "PoissonModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
PoissonModel::PoissonModel()
 : BCModel()
 , fNObs(0)
{
  DefineParameters();
};

// ---------------------------------------------------------
PoissonModel::PoissonModel(const char * name) : BCModel(name)
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
  AddParameter("lambda", 0., 7.); // index 0
}

// ---------------------------------------------------------
int PoissonModel::SetNObs(int nobs)
{
  // check that number is positive and greater than 0
  if (nobs < 0)
    return 0;

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
    lambdamin = double(nobs) - 5*sqrt(double(nobs));
    lambdamax = double(nobs) + 5*sqrt(double(nobs));
  }

  // re-set the parameter range
  SetParameterRange(0, lambdamin, lambdamax);

  // no error
  return 1;
}

// ---------------------------------------------------------
double PoissonModel::LogLikelihood(const std::vector<double> &parameters)
{
  // This methods returns the logarithm of the conditional probability
  // p(data|parameters). This is where you have to define your model.

  double logprob = 0.;

  double lambda = parameters.at(0);

  // log Poisson term
  logprob += BCMath::LogPoisson( double(fNObs), lambda );

  // return log likelihood
  return logprob;
}

// ---------------------------------------------------------
double PoissonModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  // This method returns the logarithm of the prior probability for the
  // parameters p(parameters).

  double logprob = 0.;

  // get width of the parameter range
  double dlambda = GetParameter(0)->GetRangeWidth();

  // add a flat prior probability
  logprob += log(1./dlambda); // flat prior

  // return log prior
  return logprob;
}
// ---------------------------------------------------------

