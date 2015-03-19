#include "CombinationModel.h"

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>

// ---------------------------------------------------------
CombinationModel::CombinationModel(const char* name)
    : BCModel(name)
{
    AddParameter("mass", 15., 65., "M"); // mass of a particle
    GetParameter("mass")->SetPrior(new BCGaussianPrior(39.4, 5.4)); // Gaussian prior for the mass

    AddParameter("cross section", 120., 180., "#sigma"); // cross section for a certain reaction
    GetParameter("cross section")->SetPrior(new BCGaussianPrior(150.3, 5.5)); // Gaussian prior for the cross section
};

// ---------------------------------------------------------
CombinationModel::~CombinationModel()
{
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double>& parameters)
{
    double m = parameters[0];
    double c = parameters[1];

    double mu_m  = 35.7;
    double sig_m = 3.1;

    double mu_c  = 152.3;
    double sig_c =   5.4;

    double rho = 0.7;

    double pre = 1. / 2. / TMath::Pi() / sig_m / sig_c / sqrt(1 - rho * rho);
    double exp1 = -1. / 2 / (1 - rho * rho);
    double exp2 = (m - mu_m) * (m - mu_m) / sig_m / sig_m
                  + (c - mu_c) * (c - mu_c) / sig_c / sig_c
                  - 2 * rho * (m - mu_m) * (c - mu_c) / sig_m / sig_c;

    return exp1 * exp2 + log(pre);
}
