#include "BCModelPol2.h" 

// --------------------------------------------------------- 

BCModelPol2::BCModelPol2() : BCModel()
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelPol2::BCModelPol2(const char* name) : BCModel(name)
{

	// define parameters 

	this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelPol2::DefineParameters()
{

	// define the parameters which define the model and their limits. 

	this -> AddParameter("constant", 0.0, 2.0); // index 0 
	this -> AddParameter("slope", -0.05, 0.05); // index 1 
	this -> AddParameter("quad", 0.0, 0.0008);  // index 2 
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

	double x       = datapoint -> GetValue(0); 
	double y       = datapoint -> GetValue(1); 
	double sigma_y = datapoint -> GetValue(2); 

	// define parameters and limits 

	double offset = parameters.at(0); 
	double slope  = parameters.at(1); 
	double quad   = parameters.at(2); 

	// calculate probability for a single measurement 

	return BCMath::LogGaus(y * 100.0, (offset + x * slope + x * x * quad) * 100.0, sigma_y * 100.0, true);

}

// --------------------------------------------------------- 

double BCModelPol2::LogPoissonProbability(int nentries, std::vector <double> parameters)
{

	// Poisson term is 1. => Log of it is 0.

	return 0.;

}

// --------------------------------------------------------- 

