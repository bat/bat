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

	this -> AddParameter("x1", 0.0, 100.0);  // index 0
	this -> AddParameter("x2", 0.0, 100.0);  // index 0
	this -> AddParameter("x3", 0.0, 100.0);  // index 0
	this -> AddParameter("x4", 0.0, 100.0);  // index 0
	this -> AddParameter("x5", 0.0, 100.0);  // index 0

	// debug
	std::vector<double> position; 
	
	position.push_back(1.0); 
	position.push_back(2.0); 
	position.push_back(3.0); 
	position.push_back(4.0); 
	position.push_back(5.0); 

	this -> SetMarkovChainInitialPosition(position); 

	this -> SetMarkovChainStepSize(0.01); 

}

// --------------------------------------------------------- 

double BCModelGauss::LogAPrioriProbability(std::vector<double> parameters)
{

	return 0;

}

// --------------------------------------------------------- 

double BCModelGauss::LogLikelihood(std::vector<double> parameters)
{

	double logprob = 0; 

	logprob += BCMath::LogGaus(parameters.at(0), 1.0, 0.5 * sqrt(parameters.at(0)), true); 
	logprob += BCMath::LogGaus(parameters.at(1), 2.0, 0.5 * sqrt(parameters.at(1)), true); 
	logprob += BCMath::LogGaus(parameters.at(2), 3.0, 0.5 * sqrt(parameters.at(2)), true); 
	logprob += BCMath::LogGaus(parameters.at(3), 4.0, 0.5 * sqrt(parameters.at(3)), true); 
	logprob += BCMath::LogGaus(parameters.at(4), 5.0, 0.1 * sqrt(parameters.at(4)), true); 

	return logprob; 

}

// --------------------------------------------------------- 
