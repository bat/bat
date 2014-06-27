#include "GaussModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
GaussModel::GaussModel(const char * name) : BCModel(name)
{
  AddParameter("x", -10.0, 50.0); // index 0
  AddParameter("y",  -5.0,  5.0); // index 1
  AddParameter("z",  -5.0,  5.0); // index 2

	SetPriorConstantAll();
};

// ---------------------------------------------------------
GaussModel::~GaussModel()
// default destructor
{
};

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double> &parameters)
{
  // assume a simple Gaussian Likelihood with two independent
  // variables
  double logprob = 0.;

  double x = parameters.at(0);
  double y = parameters.at(1);
  double z = parameters.at(2);

  // Gaussian Likelihood
  logprob += BCMath::LogGaus(x, 0.0, 2.0);
  logprob += BCMath::LogGaus(y, 0.0, 1.0);
  logprob += BCMath::LogGaus(z, 0.0, 1.0);

  return logprob;
}

