#include "BCModelPol1.h" 

#include "BCMath.h"

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1() : BCModel()
{

}

// --------------------------------------------------------- 

BCModelPol1::BCModelPol1(const char* name) : BCModel(name)
{

}

// --------------------------------------------------------- 

double BCModelPol1::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.; 

	// in this case we define flat a priopi probability across the whole parameter space

	double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
	double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

	logprob -= log(offset_upper - offset_lower);
	logprob -= log(slope_upper - slope_lower);

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelPol1::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

	// define data values 

	//	double x       = datapoint -> GetValue(0); 
	double y       = datapoint -> GetValue(1); 
	double sigma_y = datapoint -> GetValue(2); 

	// calculate expectation value 

	double yex = this -> FitFunction(datapoint -> GetValues(), parameters); 

	// calculate probability for a single measurement 

	return BCMath::LogGaus(y, yex, sigma_y, true); 

}

// --------------------------------------------------------- 

double BCModelPol1::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

double BCModelPol1::FitFunction(std::vector <double> x, std::vector <double> parameters)
{
	
	// get parameter values

	double offset = parameters.at(0);
	double slope  = parameters.at(1);

	return offset + x.at(0) * slope; 

}

// --------------------------------------------------------- 

