#include "BCModelPol1.h" 

#include <TMath.h> 

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1() : BCModel()
{
	// define parameters
	this -> DefineParameters();
}

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1(const char* name) : BCModel(name)
{
	// define parameters
	this -> DefineParameters();
}

// --------------------------------------------------------- 

void BCModelPol1::DefineParameters()
{
	this -> AddParameter("constant", 1.0, 3.0);  // index 0
	this -> AddParameter("slope",    0.0, 0.03); // index 1
}

// --------------------------------------------------------- 

double BCModelPol1::LogAPrioriProbability(std::vector <double> parameters)
{
	// get parameter ranges
	double offset_lower = this -> GetParameter(0) -> GetLowerLimit();
	double offset_upper = this -> GetParameter(0) -> GetUpperLimit();

	double slope_lower = this -> GetParameter(1) -> GetLowerLimit();
	double slope_upper = this -> GetParameter(1) -> GetUpperLimit();

/*
	double probability_offset = 1.0 / (offset_upper - offset_lower);
	double probability_slope = 1.0 / (slope_upper - slope_lower);
	double probability = probability_offset * probability_slope;
*/

	// calculate natural logarithm of the flat probability
	// this is equivalent to logarithm of the probability calculated
	// useing the three lines of code above

	double logprob = 0.;

	logprob -= TMath::Log(offset_upper - offset_lower);
	logprob -= TMath::Log(slope_upper - slope_lower);

	return logprob;
}

// --------------------------------------------------------- 

double BCModelPol1::LogLikelihood(std::vector <double> parameters)
{
	// define log of probability
	double logprob = 0.0;

	// get parameter values
	double offset = parameters.at(0);
	double slope  = parameters.at(1);

	// loop over all data points
	int npoints = this -> GetNDataPoints();

	for (int ipoint = 0; ipoint < npoints; ipoint++)
	{
		// define data values
		double x       = this -> GetDataPoint(ipoint) -> GetValue(0);
		double y       = this -> GetDataPoint(ipoint) -> GetValue(1);
		double sigma_y = this -> GetDataPoint(ipoint) -> GetValue(2);

		// calculate probability assuming a Gaussian distribution for each point
		logprob += BCMath::LogGaus(y, (offset + x * slope), sigma_y, true);
	}

	return logprob; 
}

// --------------------------------------------------------- 
