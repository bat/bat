#include "BCModelPol2.h" 

#include "BCMath.h"

// --------------------------------------------------------- 

BCModelPol2::BCModelPol2() : BCModel()
{

}

// --------------------------------------------------------- 

BCModelPol2::BCModelPol2(const char* name) : BCModel(name)
{

}

// --------------------------------------------------------- 

double BCModelPol2::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.; 

	double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
	double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

	double quad_lower =  this -> GetParameter(2) -> GetLowerLimit(); 
	double quad_upper =  this -> GetParameter(2) -> GetUpperLimit(); 

	// in this case we define flat a priopi probability across the whole parameter space

	logprob -= log(offset_upper - offset_lower);
	logprob -= log(slope_upper - slope_lower);
	logprob -= log(quad_upper - quad_lower);

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelPol2::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
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

double BCModelPol2::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

double BCModelPol2::FitFunction(std::vector <double> x, std::vector <double> parameters)
{
	
	// get parameter values

	double offset = parameters.at(0); 
	double slope  = parameters.at(1); 
	double quad   = parameters.at(2); 

	return offset + x.at(0) * slope + x.at(0) * x.at(0) * quad; 

}

// --------------------------------------------------------- 

