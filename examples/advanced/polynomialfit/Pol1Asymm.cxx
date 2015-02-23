#include "Pol1Asymm.h"

#include <BAT/BCDataPoint.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <TString.h>

#include <cmath>

// ---------------------------------------------------------
Pol1Asymm::Pol1Asymm(const char * name)
	: BCFitter(name)
{
	DefineParameters();
}

// ---------------------------------------------------------
Pol1Asymm::~Pol1Asymm() {
}

// ---------------------------------------------------------
void Pol1Asymm::DefineParameters() {
	AddParameter("p0", -0.2   ,  1.2);    // index 0
	AddParameter("p1",  0.015 ,  0.045);   // index 1
	
	SetPriorConstantAll();
	
	PrintSummary();
}

// ---------------------------------------------------------
double Pol1Asymm::FitFunction(const std::vector<double> & x, const std::vector<double> & par) {
	// first-order polynomial
	double r = par[0];

	// n-th order polynomial
	double X = 1;
	for (unsigned i=1; i<par.size(); ++i) {
		X *= x[0];
		r += par[i]*X;
	}

	return r;
}

// ---------------------------------------------------------
double Pol1Asymm::LogLikelihood(const std::vector<double> & par) {
   double logl = 0.;

   // loop over the data points
   for(unsigned i=0 ; i < GetNDataPoints(); i++) {

		 // get data point
		 std::vector<double> x = GetDataSet() -> GetDataPoint(i) -> GetValues();
		 double y    = x[1];
		 double eylo = x[2];
		 double eyhi = x[3];
		 
		 // calculate the expectation value of the function at this point
		 double yexp = FitFunction(x,par);
		 
		 // We define asymmetric errors as two half gaussians with different
		 // widths around the expectation value. In this definition the central
		 // value remains central (it is equal to mean) and the uncertainties
		 // represent central 68% interval.
		 // Other definitions are possible.
		 
		 // if measured value is below the expectation value we use 'eylo' as
		 // the width of the half-gaussian, if it is above we use 'eyhi'
		 double yerr = (y > yexp) ? eylo : eyhi;
		 
		 logl += BCMath::LogGaus(y, yexp, yerr, true);
   }

   return logl;
}
