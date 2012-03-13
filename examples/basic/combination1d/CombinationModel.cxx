#include "CombinationModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel() : BCModel()
{
  DefineParameters();
};

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char * name) : BCModel(name)
{
  DefineParameters();
};

// ---------------------------------------------------------
CombinationModel::~CombinationModel()
{
};

// ---------------------------------------------------------
void CombinationModel::DefineParameters()
{
  AddParameter("mass", 15.0, 65.0); // mass of a particle
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob = 0.;

  double mass = parameters.at(0);

  logprob += BCMath::LogGaus(mass, 35.7, 3.1);

  return logprob;
}

// ---------------------------------------------------------
double CombinationModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  double logprob = 0.;

  double mass = parameters.at(0);

  logprob += BCMath::LogGaus(mass, 39.4, 5.4); // Gaussian prior for the mass

  return logprob;
}
// ---------------------------------------------------------

