#include "BCModelPol1.h" 

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

	// add two parameters which define a line

	this -> AddParameter("constant", 1.0, 3.0);  // index 0
	this -> AddParameter("slope",    0.0, 0.03); // index 1

}

// --------------------------------------------------------- 

double BCModelPol1::LogAPrioriProbability(std::vector <double> parameters)
{

	// get parameter ranges

	double offset_lower = this -> GetParameter(0) -> GetLowerLimit();
	double offset_upper = this -> GetParameter(0) -> GetUpperLimit();

	double slope_lower = this -> GetParameter(1) -> GetLowerLimit();
	double slope_upper = this -> GetParameter(1) -> GetUpperLimit();

	// calculate natural logarithm of the flat probability

	double logprob = 0.;

	logprob -= log(offset_upper - offset_lower);
	logprob -= log(slope_upper - slope_lower);

	return logprob;

}

// --------------------------------------------------------- 

double BCModelPol1::LogLikelihood(std::vector <double> parameters)
{

	// define log of probability

	double logprob = 0.0;

	// loop over all data points

	int npoints = this -> GetNDataPoints();

	for (int ipoint = 0; ipoint < npoints; ipoint++)
		{
			// define data values

			//			double x       = this -> GetDataPoint(ipoint) -> GetValue(0);
			double y       = this -> GetDataPoint(ipoint) -> GetValue(1);
			double sigma_y = this -> GetDataPoint(ipoint) -> GetValue(2);

			// calculate expectation value 

			double yex = this -> FitFunction(this -> GetDataPoint(ipoint) -> GetValues(), parameters); 

			// calculate probability assuming a Gaussian distribution for each point

			logprob += BCMath::LogGaus(y, yex, sigma_y, true);
		}

	return logprob; 

}

// --------------------------------------------------------- 

double BCModelPol1::FitFunction(std::vector <double> x, std::vector <double> parameters)
{
	
	// get parameter values

	double offset = parameters.at(0);
	double slope  = parameters.at(1);

	return offset + x.at(0) * slope; 

}

// --------------------------------------------------------- 
