#include "PoissonModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
PoissonModel::PoissonModel() : BCModel()
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
  AddParameter("lambda", 0., 7.); // index 0
}

// ---------------------------------------------------------
int PoissonModel::SetNObs(int nobs)
{
  // check number 
  if (nobs < 0)
    {
      return 0; 
    }

  // set number of observed events
  fNObs = nobs; 

  // adjust ranges
  double lambdamin = 0;
  double lambdamax = 10;

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

  // set parameter range
  SetParameterRange(0, lambdamin, lambdamax);

  // no error
  return 1;
}

// ---------------------------------------------------------
double PoissonModel::LogLikelihood(std::vector <double> parameters)
{
  // This methods returns the logarithm of the conditional probability
  // p(data|parameters). This is where you have to define your model.

  double logprob = 0.;

  double lambda = parameters.at(0);

  // log Poisson term
  logprob += BCMath::LogPoisson(double(fNObs), lambda); 

  return logprob;
}

// ---------------------------------------------------------
double PoissonModel::LogAPrioriProbability(std::vector <double> parameters)
{
  // This method returns the logarithm of the prior probability for the
  // parameters p(parameters).

  double logprob = 0.;

  double dlambda = GetParameter(0)->GetRangeWidth(); 

  logprob += log(1./dlambda); // flat prior

  return logprob;
}
// ---------------------------------------------------------

