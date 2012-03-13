
#include <TString.h>

#include <BAT/BCLog.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>
#include <BAT/BCMath.h>

#include "Pol1Asymm.h"

// ---------------------------------------------------------
Pol1Asymm::Pol1Asymm() : BCModel()
{ // default constructor
}

// ---------------------------------------------------------
Pol1Asymm::Pol1Asymm(const char * name) : BCModel(name)
{ // constructor
}

// ---------------------------------------------------------
Pol1Asymm::~Pol1Asymm()
{ // destructor
}

// ---------------------------------------------------------
void Pol1Asymm::DefineParameters()
{
	// Adding all the parameters of the model and their allowed ranges
	// Keep track of index of the parameter as you need it to be able
	// to use the parameter later on. Usage by using the parameter
	// name is also possible but due to the array searching it is very
	// CPU time expensive and therefore not recommended in routines
	// called for every iteration, like LogLikelihood() or
	// LogAPrioriProbability()
	this -> AddParameter("p0", -2.  ,  2.);    // index 0
	this -> AddParameter("p1", -0.05,  0.05);  // index 1

	// Print parameter summary
	BCLog::OutSummary(
			Form("Model \'%s\' has %d parameters:",this->GetName().data(),this -> GetNParameters()));
	for(unsigned int i=0; i< this -> GetNParameters(); i++)
		BCLog::OutSummary(Form("   %d. %s    range: %g - %g",
				i,
				this -> GetParameter(i) -> GetName().data(),
				this -> GetParameter(i) -> GetLowerLimit(),
				this -> GetParameter(i) -> GetUpperLimit() ) );
}

// ---------------------------------------------------------
double Pol1Asymm::FitFunction(const std::vector<double> & x, const std::vector<double> & par)
{
	// function to fit with
	// get the parameters of the function
	double p0 = par[0];
	double p1 = par[1];

	// calculate the value of the function at x[0] for a given set of parameters
	return p0 + p1*x[0];
}

// ---------------------------------------------------------
double Pol1Asymm::LogLikelihood(const std::vector<double> & par)
{
	double logl = 0.;

	// loop over the data points
	for(int i=0 ; i < this -> GetNDataPoints(); i++)
	{
		// get data point
		std::vector<double> x = this -> GetDataPoint(i) -> GetValues();
		double y    = x[1];
		double eylo = x[2];
		double eyhi = x[3];

		// calculate the expectation value of the function at this point
		double yexp = this -> FitFunction(x,par);

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

// ---------------------------------------------------------
double Pol1Asymm::LogAPrioriProbability(const std::vector<double> & parameters)
{
	// Definition of the Prior of the model.
	// For flat prior it's very easy.
	double logprob = 0.;
	for(unsigned int i=0; i < this -> GetNParameters(); i++)
		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}

// ---------------------------------------------------------
