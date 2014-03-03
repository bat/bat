#include "CombinationModel.h"

#include <TMath.h>
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
  AddParameter("mass", 15., 65.); // mass of a particle
  AddParameter("cross section", 120., 180.); // cross section for a certain reaction
}

// ---------------------------------------------------------
double CombinationModel::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob = 0.;

  double m = parameters.at(0);
  double c = parameters.at(1);

  double mu_m  = 35.7;
  double sig_m = 3.1;

  double mu_c  = 152.3;
  double sig_c =   5.4;

  double rho = 0.7;

  double pre = 1./2./TMath::Pi()/sig_m/sig_c/sqrt(1-rho*rho);
  double exp1 = -1./2/(1-rho*rho);
  double exp2 = (m-mu_m)*(m-mu_m)/sig_m/sig_m
    + (c-mu_c)*(c-mu_c)/sig_c/sig_c
    - 2*rho*(m-mu_m)*(c-mu_c)/sig_m/sig_c;

  double gaus = pre*exp(exp1*exp2);

  logprob += log(gaus);

  return logprob;
}

// ---------------------------------------------------------
double CombinationModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  double logprob = 0.;

  double mass = parameters.at(0);
  double crosssection = parameters.at(1);

  logprob += BCMath::LogGaus(mass, 39.4, 5.4); // Gaussian prior for the mass
  logprob += BCMath::LogGaus(crosssection, 150.3, 5.5); // Gaussian prior for the mass

  return logprob;
}
// ---------------------------------------------------------

