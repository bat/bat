#include "GaussModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <TRandom3.h>

#include <cmath>
#include <iostream>

// ---------------------------------------------------------
GaussModel::GaussModel() : BCModel()
{
  // default constructor
  DefineParameters();
};

// ---------------------------------------------------------
GaussModel::GaussModel(const char * name) : BCModel(name)
{
  // constructor
  DefineParameters();
};

// ---------------------------------------------------------
GaussModel::~GaussModel()
// default destructor
{
};

// ---------------------------------------------------------
void GaussModel::DefineParameters()
{
  // add parameters x and y
  AddParameter("x", -10.0, 50.0); // index 0
  AddParameter("y",  -5.0,  5.0); // index 1
}

// ---------------------------------------------------------
double GaussModel::LogLikelihood(const std::vector<double> &parameters)
{
  // assume a simple Gaussian Likelihood with two independent
  // variables
  double logprob = 0.;

  double x = parameters.at(0);
  double y = parameters.at(1);

  // Gaussian Likelihood
  logprob += BCMath::LogGaus(x, 0.0, 2.0);
  logprob += BCMath::LogGaus(y, 0.0, 1.0);

  return logprob;
}

// ---------------------------------------------------------
double GaussModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  // assume flat prior in both variables
  double logprob = 0.;

  double dx = GetParameter(0)->GetRangeWidth();
  double dy = GetParameter(1)->GetRangeWidth();

  logprob += log(1./dx); // flat prior for x
  logprob += log(1./dy); // flat prior for y

  return logprob;
}

// --------------------------------------------------------
double GaussModel::MCMCTrialFunctionSingle(unsigned int ichain, unsigned int ipar)
{
  // no check of range for performance reasons

  // get scale factor from an array of scale factors. the size of the
  // array is number of chains times number of parameters.
  // double scale = fMCMCTrialFunctionScaleFactor[ichain * GetNParameters() + ipar];

  // choose trial function by uncommenting any of the lines below

  // Gaussian with fixed width
  //  return fRandom->Gaus(0.0, 1.0);

  // Gaussian with adjustable width
  //  return fRandom->Gaus(0.0, scale);

  // Breit-Wigner with adjustable width
  //  return fRandom->BreitWigner(0.0, scale);

  // Flat function with adjustable width
  //  return scale * 2. * (0.5 - fRandom->Uniform());

  // Flat function with fixed width
  return 0.02 * (0.5 - fRandom->Uniform());

}
// ---------------------------------------------------------

