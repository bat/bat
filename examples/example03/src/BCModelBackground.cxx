#include "BCModelBackground.h" 

#include "BCMath.h"

// --------------------------------------------------------- 

BCModelBackground::BCModelBackground() : BCModel() 
{
  
	// define parameters

	this -> DefineParameters();

}

// ---------------------------------------------------------

BCModelBackground::BCModelBackground(const char* name) : BCModel(name)
{

	// define parameters

	this -> DefineParameters();

}

// ---------------------------------------------------------

void BCModelBackground::DefineParameters()
{

	// add the one parameter which defines the model 

	this -> AddParameter("background", 0.0, 300.0);  // index 0

}

// --------------------------------------------------------- 

double BCModelBackground::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.0; 

	//  double constant = parameters.at(0); 

	double background_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double background_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	// calculate probability 

	logprob -= log(background_upper - background_lower); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelBackground::LogLikelihood(std::vector <double> parameters)
{

	double logprob = 0.0; 

	// get parameter values 

	double background = parameters.at(0); 

	// loop over all data points 

	int nbins = this -> GetNDataPoints() - 1; 

	// get histogram parameters 

	double Emin = this -> GetDataPoint(0) -> GetValue(0); 
	double Emax = this -> GetDataPoint(0) -> GetValue(1); 
	double dE   = (Emax - Emin) / double(nbins); 

	for (int ipoint = 1; ipoint < nbins + 1; ipoint++)
		{      
			// define data values 

			double energy = this -> GetDataPoint(ipoint) -> GetValue(0); 
			double events = this -> GetDataPoint(ipoint) -> GetValue(1); 

			// calculate probability 

			double expected = background * this -> integral_f_B(energy, dE, nbins) * dE; 

			logprob += events * log (expected) - BCMath::ApproxLogFact(events) - expected; 
		}

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelBackground::integral_f_B(double E, double dE, int nbins)
{

	return 1.0 / double(nbins) / dE; 

}

// --------------------------------------------------------- 

double BCModelBackground::FitFunction(std::vector <double> x, std::vector <double> parameters)
{

	int nbins = int(this -> GetNDataPoints() - 1); 

	return parameters.at(0)/double(nbins); 
	
}

// --------------------------------------------------------- 


