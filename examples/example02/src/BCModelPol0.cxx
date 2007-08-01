#include "BCModelPol0.h" 

#include <TMath.h> 

// --------------------------------------------------------- 

BCModelPol0::BCModelPol0() : BCModel()
{

  // define parameters 

  this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelPol0::BCModelPol0(const char* name) : BCModel(name)
{

  // define parameters 

  this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelPol0::DefineParameters()
{

  this -> AddParameter("constant", 0.0, 4.0); // index 0 

}

// --------------------------------------------------------- 

double BCModelPol0::APrioriProbability(std::vector <double> parameters)
{

  double probability = 1.0; 

  double constant_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double constant_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  double probability_constant = 1.0 / (constant_upper - constant_lower); 

  probability = probability_constant; 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol0::ConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

  double probability = 1.0; 

  // define data values 

  double y       = datapoint -> GetValue(1); 
  double sigma_y = datapoint -> GetValue(2); 

  // define parameters and limits 

  double constant = parameters.at(0); 

  // calculate probability for a single measurement 
  
  probability = TMath::Gaus(y * 100.0, constant * 100.0, sigma_y * 100.0, true); 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol0::PoissonProbability(int nentries, std::vector <double> parameters)
{

  double probability = 1.0; 

  return probability; 

}

// --------------------------------------------------------- 

