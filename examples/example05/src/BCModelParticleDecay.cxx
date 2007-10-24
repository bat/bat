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

	this -> AddParameter("mass",   0.0, 400.0);    // index 0 
	this -> AddParameter("energy", 0.0, 200.0); // index 1

	// debug
	// only one energy since they are correlated due to the mass

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.; 

	// in this case we define flat a priopi probability across the whole parameter space

	double mass_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double mass_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	logprob -= log(mass_upper - mass_lower);

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

	double logprob = 0.0; 

	// define data values 

	double E1 = datapoint -> GetValue(0); 
	double E2 = datapoint -> GetValue(1); 

	// get parameters 

	// add resolution 

	logprob += this -> logW(E1, parameters.at(1)); 
	logprob += this -> logW(E2, parameters.at(1)); 

	// Breit-Wigner 

	double Esum = 2.0 * parameters.at(1); 

	logprob += BCMath::LogBreitWignerNonRel(Esum, parameters.at(0), 2.0, true); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelParticleDecay::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

double BCModelParticleDecay::logW(double Emeasured, double Etrue)
{

	double logprob; 

	double sigma = 0.5 * sqrt(Etrue); 

	logprob = BCMath::LogGaus(Emeasured, Etrue, sigma, true); 

	return logprob; 

}

// --------------------------------------------------------- 
