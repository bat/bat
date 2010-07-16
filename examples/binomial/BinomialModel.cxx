#include "BinomialModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
BinomialModel::BinomialModel() : BCModel()
			       , fNTotal(10)
			       , fNSelected(5)
{  
  DefineParameters();
};

// ---------------------------------------------------------
BinomialModel::BinomialModel(const char * name) : BCModel(name)
{ 
  DefineParameters();
};

// ---------------------------------------------------------
BinomialModel::~BinomialModel()
{
};

// ---------------------------------------------------------
void BinomialModel::DefineParameters()
{
  AddParameter("epsilon", 0., 1.); // the efficiency
}

// ---------------------------------------------------------
double BinomialModel::LogLikelihood(std::vector <double> parameters)
{
  double logprob = 0.;

  double eps = parameters.at(0);

  logprob += BCMath::LogApproxBinomial(fNTotal, fNSelected, eps);

  return logprob;
}

// ---------------------------------------------------------
double BinomialModel::LogAPrioriProbability(std::vector <double> parameters)
{
  double logprob = 0.;

  return logprob;
}
// ---------------------------------------------------------

