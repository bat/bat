#include "BCModelPol1.h" 

#include <TMath.h> 

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
  this -> AddParameter("constant", 1.0, 3.0);  // index 0 
  this -> AddParameter("slope",    0.0, 0.03); // index 1 

}

// --------------------------------------------------------- 

double BCModelPol1::APrioriProbability(std::vector <double> parameters)
{

  // get parameter ranges 

  double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
  double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

  // calculate probabilities 

  double probability_offset = 1.0 / (offset_upper - offset_lower); 
  double probability_slope = 1.0 / (slope_upper - slope_lower);   

  double probability = probability_offset * probability_slope; 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol1::Likelihood(std::vector <double> parameters)
{

  // define log of probability 

  double logprobability = 0.0; 

  // get parameter values 

  double offset = parameters.at(0); 
  double slope  = parameters.at(1); 

  // loop over all data points 

  int npoints = this -> GetNDataPoints(); 

  for (int ipoint = 0; ipoint < npoints; ipoint++)
    {

      // define data values 

      double x       = this -> GetDataPoint(ipoint) -> GetValue(0); 
      double y       = this -> GetDataPoint(ipoint) -> GetValue(1); 
      double sigma_y = this -> GetDataPoint(ipoint) -> GetValue(2); 
      
      // calculate probability assuming a Gaussian distribution for each point 

      logprobability += TMath::Log10(TMath::Gaus(y, (offset + x * slope), sigma_y, true)); 
    }

  // calculate probability 

  double probability = TMath::Power(10.0, logprobability); 

  return probability; 

}

// --------------------------------------------------------- 
