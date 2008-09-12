#include "BCModelParticleDecay.h" 

// --------------------------------------------------------- 

BCModelParticleDecay::BCModelParticleDecay() : BCModel()
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelParticleDecay::BCModelParticleDecay(const char* name) : BCModel(name)
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelParticleDecay::DefineParameters()
{

	// adds the parameters which define the model and their limits. 

	this -> AddParameter("energy", 0.0, 100.0); // index 0
	this -> AddParameter("mass",   0.0, 100.0); // index 1 
	this -> AddParameter("width",  0.0, 5.0);   // index 2

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogAPrioriProbability(std::vector <double> parameters)
{

	return 0; 

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

	double logprob = 0.0; 

	// add resolution 

	logprob += this -> logEnergyResolution(datapoint -> GetValue(0), parameters.at(0)); 

	// Breit-Wigner 

	//	logprob += BCMath::LogBreitWignerNonRel(parameters.at(1), parameters.at(0), 2.0, true); 

	logprob += this -> logMassDistribution(parameters.at(0), parameters.at(1), parameters.at(2)); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

double BCModelParticleDecay::logEnergyResolution(double Emeasured, double Etrue)
{

	double logprob; 

	//	double sigma = 0.5 * sqrt(Etrue); 

	double sigma = 5.0; 

	logprob = BCMath::LogGaus(Emeasured, Etrue, sigma, true); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelParticleDecay::logMassDistribution(double mass, double polemass, double width)
{

	double logprob; 

	logprob = BCMath::LogGaus(mass, polemass, width, true); 

	return logprob; 

}

// --------------------------------------------------------- 
