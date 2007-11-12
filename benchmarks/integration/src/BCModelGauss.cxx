#include "BCModelGauss.h" 

#include <math.h>

// --------------------------------------------------------- 

BCModelGauss::BCModelGauss() : BCModel()
{

	// define parameters

	this -> DefineParameters();

}

// --------------------------------------------------------- 

BCModelGauss::BCModelGauss(const char* name) : BCModel(name)
{

	// define parameters

	this -> DefineParameters();

}

// --------------------------------------------------------- 

void BCModelGauss::DefineParameters()
{

	this -> AddParameter("xrew", -10.0, 10.0);  // index 0
	this -> AddParameter("y", -10.0, 10.0);  // index 0
	this -> AddParameter("z", -10.0, 10.0);  // index 0

}

// --------------------------------------------------------- 

double BCModelGauss::LogAPrioriProbability(std::vector <double> parameters)
{

	return 0;

}

// --------------------------------------------------------- 

double BCModelGauss::LogLikelihood(std::vector <double> parameters)
{

	double logprob = 0; 

	logprob += BCMath::LogGaus(parameters.at(0), 0.0, 0.5, true);
	logprob += BCMath::LogGaus(parameters.at(1), 0.0, 0.5, true);
	logprob += BCMath::LogGaus(parameters.at(2), 0.0, 0.5, true);
	logprob += log(5.0); 

	return logprob; 

}

// --------------------------------------------------------- 
