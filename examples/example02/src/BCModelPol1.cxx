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
	this -> AddParameter("constant", 0.0, 4.0);   // index 0 
	this -> AddParameter("slope",   0.02,  0.04); // index 1 
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

	logprob -= TMath::Log(offset_upper - offset_lower);
	logprob -= TMath::Log(slope_upper - slope_lower);

	return logprob; 
}

// --------------------------------------------------------- 

double BCModelPol1::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{
	// define data values 
	double x       = datapoint -> GetValue(0); 
	double y       = datapoint -> GetValue(1); 
	double sigma_y = datapoint -> GetValue(2); 

	// define parameters and limits 
	double offset = parameters.at(0); 
	double slope  = parameters.at(1); 

	// calculate probability for a single measurement 
	return BCMath::LogGaus(y * 100.0, (offset + x * slope) * 100.0 , sigma_y * 100.0, true); 
}

// --------------------------------------------------------- 

double BCModelPol1::LogPoissonProbability(int nentries, std::vector <double> parameters)
{
	// Poisson term is 1. => Log of it is 0.
	return 0.;
}

// --------------------------------------------------------- 

