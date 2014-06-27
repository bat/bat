#include "CombinationModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char * name) : BCModel(name)
{
  AddParameter("mass", 15.0, 65.0); // mass of a particle
	SetPriorGauss("mass",39.4, 5.4);	// Gaussian prior for the mass
};

// ---------------------------------------------------------
CombinationModel::~CombinationModel()
{
};

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob = 0.;

  double mass = parameters.at(0);

  logprob += BCMath::LogGaus(mass, 35.7, 3.1);

  return logprob;
}
