#include "RatioModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
RatioModel::RatioModel(const char * name)
 : BCModel(name)
{
  // define the parameters x and y
  AddParameter("x", 0., 8.); // index 0
  AddParameter("y", 0., 16.); // index 1

	SetPriorConstantAll();

	AddObservable("r", 0, 2, "#frac{x}{y}");
};

// ---------------------------------------------------------
RatioModel::~RatioModel()
{
};

// ---------------------------------------------------------
double RatioModel::LogLikelihood(const std::vector<double> &parameters)
{
  // This methods returns the logarithm of the conditional probability
  // p(data|parameters). This is where you have to define your model.

  double logprob = 0.;

  double x = parameters[0];
  double y = parameters[1];

  // calculate the Gaussian probability densities
  logprob += BCMath::LogGaus(x, 4., 1.);
  logprob += BCMath::LogGaus(y, 8., 2.0);

  // return the log likelihood
  return logprob;
}

// ---------------------------------------------------------
void RatioModel::CalculateObservables(const std::vector<double> &parameters)
{
	// store ratio
	GetObservable(0) -> Value( (parameters[1]!=0) ? parameters[0] / parameters[1] : 0 );
}

