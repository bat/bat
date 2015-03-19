#include "CombinationModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char* name)
    : BCModel(name)
{
    AddParameter("mass", 15.0, 65.0); // mass of a particle
    GetParameter("mass")->SetPrior(new BCGaussianPrior(39.4, 5.4)); // Gaussian prior for the mass
}

// ---------------------------------------------------------
CombinationModel::~CombinationModel()
{
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double>& parameters)
{
    return BCMath::LogGaus(parameters[0], 35.7, 3.1);
}
