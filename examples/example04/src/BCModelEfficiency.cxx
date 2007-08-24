#include "BCModelEfficiency.h" 

#include <BCMath.h> 

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency() : BCModel()
{
	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency(const char* name) : BCModel(name)
{
	// define parameters 
	this -> DefineParameters(); 
}

// --------------------------------------------------------- 

void BCModelEfficiency::DefineParameters()
{
	this -> AddParameter("efficiency", 0.0,  1.0); // index 0 
}

// --------------------------------------------------------- 

double BCModelEfficiency::LogAPrioriProbability(std::vector <double> parameters)
{
	// get parameter ranges 
	double eff_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double eff_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	// calculate probabilities 
	double logprob = 0.;

	logprob -= log(eff_upper - eff_lower);

	return logprob;
}

// --------------------------------------------------------- 

double BCModelEfficiency::LogLikelihood(std::vector <double> parameters)
{
	// get parameters
	double eff = parameters.at(0);

	// get data values
	int k = BCMath::Nint(this -> GetDataPoint(0) -> GetValue(0));
	int n = BCMath::Nint(this -> GetDataPoint(0) -> GetValue(1));

	return BCMath::LogApproxBinomial(n, k, eff);
}

// --------------------------------------------------------- 
