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

  this -> AddParameter("constant", 0.0, 4.0);   // index 0 
  this -> AddParameter("slope",   0.02,  0.04); // index 1 

}

// --------------------------------------------------------- 

double BCModelPol1::APrioriProbability(std::vector <double> parameters)
{

  double probability = 1.0; 

  //  double offset = parameters.at(0); 
  //  double slope  = parameters.at(0); 

  double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
  double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

  double probability_offset = 1.0 / (offset_upper - offset_lower); 
  double probability_slope = 1.0 / (slope_upper - slope_lower); 

  probability = probability_offset * probability_slope; 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol1::ConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

  double probability = 1.0; 

  // define data values 

  double x       = datapoint -> GetValue(0); 
  double y       = datapoint -> GetValue(1); 
  double sigma_y = datapoint -> GetValue(2); 

  // define parameters and limits 

  double offset = parameters.at(0); 
  double slope  = parameters.at(1); 

  // calculate probability for a single measurement 
  
  probability = TMath::Gaus(y * 100.0, (offset + x * slope) * 100.0 , sigma_y * 100.0, true); 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol1::PoissonProbability(int nentries, std::vector <double> parameters)
{

  double probability = 1.0; 

  return probability; 

}

// --------------------------------------------------------- 

