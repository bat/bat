#include "PoissonModel.h"

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
PoissonModel::PoissonModel(const char * name)
	:	BCModel(name)
	, fNObs(0)
{
	// add a parameter for the number of expected events. The range will
	// be adjusted later according to the number of observed events.
	AddParameter("lambda", 0., 7., "#lambda"); // index 0
	GetParameter("lambda") -> SetPriorConstant();
}

// ---------------------------------------------------------
PoissonModel::~PoissonModel() {
}

// ---------------------------------------------------------
void PoissonModel::SetNObs(unsigned nobs) {
   // set number of observed events
   fNObs = nobs;

   // adjust parameter ranges
   double lambdamin = 0;
   double lambdamax = 10;

   // the adjustment depends on the number of observed events
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
      lambdamin = double(nobs) - 5 * sqrt(double(nobs));
      lambdamax = double(nobs) + 5 * sqrt(double(nobs));
   }

   // re-set the parameter range
   GetParameter(0) -> SetLimits(lambdamin, lambdamax);
}

// ---------------------------------------------------------
double PoissonModel::LogLikelihood(const std::vector<double> & parameters) {
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   // log Poisson term
	 return BCMath::LogPoisson(fNObs, parameters[0]);
}
