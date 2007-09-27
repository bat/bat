#include <TMath.h>

#include "BCModelSignal.h" 

// --------------------------------------------------------- 

BCModelSignal::BCModelSignal() : BCModel() 
{
  
	// define parameters

	this -> DefineParameters();

}

// ---------------------------------------------------------

BCModelSignal::BCModelSignal(const char* name) : BCModel(name)
{

	// define parameters

	this -> DefineParameters();

}

// ---------------------------------------------------------

void BCModelSignal::DefineParameters()
{

	// define the two parameters which define the model 

	this -> AddParameter("background", 0.0, 200.0);  // index 0
	this -> AddParameter("signal",     0.0, 100.0);  // index 1

}

// --------------------------------------------------------- 

double BCModelSignal::LogAPrioriProbability(std::vector <double> parameters)
{

	double logprob = 0.0; 

	double background_lower = this -> GetParameter(0) -> GetLowerLimit(); 
	double background_upper = this -> GetParameter(0) -> GetUpperLimit(); 

	double signal_lower = this -> GetParameter(1) -> GetLowerLimit(); 
	double signal_upper = this -> GetParameter(1) -> GetUpperLimit(); 

	// calculate probability 

	logprob -= log(background_upper - background_lower); 
	logprob -= log(signal_upper - signal_lower); 

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelSignal::LogLikelihood(std::vector <double> parameters)
{

	double logprob = 0.0; 

	// get parameter values 

	double background = parameters.at(0); 
	double signal     = parameters.at(1); 

	// loop over all data points 

	int nbins = this -> GetNDataPoints(); 

	// get histogram parameters 

	double Emin = this -> GetDataPoint(0) -> GetValue(0); 
	double Emax = this -> GetDataPoint(0) -> GetValue(1); 
	double dE   = (Emax - Emin) / double(nbins); 

	for (int ipoint = 1; ipoint < nbins; ipoint++)
		{      
			// define data values 

			double energy = this -> GetDataPoint(ipoint) -> GetValue(0); 
			double events = this -> GetDataPoint(ipoint) -> GetValue(1); 

			// calculate probability 

			double expected = (background * this -> integral_f_B(energy, dE, nbins) * dE 
												 + signal * this -> integral_f_S(energy, dE, nbins)); 
			
			logprob += events * log (expected) - BCMath::ApproxLogFact(events) - expected; 
		}

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelSignal::integral_f_B(double E, double dE, int nbins)
{

	return 1.0 / double(nbins) / dE; 

}

// --------------------------------------------------------- 

double BCModelSignal::integral_f_S(double E, double dE, int nbins)
{

	double mean  = 2039.0; 
	double sigma =    2.5; 
	double t1     = (E - mean) / (sqrt(2.0) * sigma); 
	double t2     = (E + dE - mean) / (sqrt(2.0) * sigma); 

	return 0.5 * (TMath::Erf(t2) - TMath::Erf(t1)); 

}

// --------------------------------------------------------- 

