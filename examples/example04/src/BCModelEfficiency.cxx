#include "BCModelEfficiency.h" 

#include <TMath.h> 

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
//	this -> AddParameter("lambda",  0.0, 20.0); // index 1 
}

// --------------------------------------------------------- 

double BCModelEfficiency::LogAPrioriProbability(std::vector <double> parameters)
{
	// get parameter ranges 
	double eff_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double eff_upper = this -> GetParameter(0) -> GetUpperLimit(); 

//	double lambda_lower = this -> GetParameter(1) -> GetLowerLimit(); 
//	double lambda_upper = this -> GetParameter(1) -> GetUpperLimit(); 

	// calculate probabilities 
	double logprob = 0.;

	logprob -= TMath::Log(eff_upper - eff_lower);

	return logprob;
}

// --------------------------------------------------------- 

double BCModelEfficiency::LogLikelihood(std::vector <double> parameters)
{
	// get parameters
	double eff = parameters.at(0);
//	double lambda  = parameters.at(1);

	// get data values
	int k = TMath::Nint(this -> GetDataPoint(0) -> GetValue(0));
	int n = TMath::Nint(this -> GetDataPoint(0) -> GetValue(1));

//	double poisson  = TMath::Poisson(double(y), lambda);
//	double binomial = 1.0;

//	if (y < 20)
//		binomial = TMath::Binomial(y, x) * TMath::Power(epsilon, x) * TMath::Power(1.0 - epsilon, y - x);

//	return TMath::Log(poisson) + TMath::Log(binomial);

	return BCMath::LogApproxBinomial(n, k, eff);
}

// --------------------------------------------------------- 
