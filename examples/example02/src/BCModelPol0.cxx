#include "BCModelPol0.h" 

#include <TMath.h> 

// --------------------------------------------------------- 

BCModelPol0::BCModelPol0() : BCModel()
{
	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 

BCModelPol0::BCModelPol0(const char* name) : BCModel(name)
{
	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 

void BCModelPol0::DefineParameters()
{
	this -> AddParameter("constant", 0.0, 4.0); // index 0 
}

// --------------------------------------------------------- 

double BCModelPol0::LogAPrioriProbability(std::vector <double> parameters)
{
	double logprob = 0.; 

	// in this case we define flat a priopi probability across the whole parameter space

	double constant_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double constant_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	logprob -= TMath::Log(constant_upper - constant_lower);

	return logprob;
}

// --------------------------------------------------------- 

double BCModelPol0::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{
	// define data values 
	double y       = datapoint -> GetValue(1); 
	double sigma_y = datapoint -> GetValue(2); 

	// define parameters and limits 
	double constant = parameters.at(0); 

	// calculate probability for a single measurement 
	return BCMath::LogGaus(y * 100.0, constant * 100.0, sigma_y * 100.0, true); 
}

// --------------------------------------------------------- 

double BCModelPol0::LogPoissonProbability(int nentries, std::vector <double> parameters)
{
	// Poisson term is 1. => Log of it is 0.
	return 0.;
}

// --------------------------------------------------------- 

