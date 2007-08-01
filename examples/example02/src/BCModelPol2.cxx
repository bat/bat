#include "BCModelPol2.h" 

#include <TMath.h> 

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

  this -> AddParameter("constant", 0.0, 2.0); // index 0 
  this -> AddParameter("slope", -0.05, 0.05); // index 1 
  this -> AddParameter("quad", 0.0, 0.0008);  // index 2 

}

// --------------------------------------------------------- 

double BCModelPol2::APrioriProbability(std::vector <double> parameters)
{

  double probability = 1.0; 

  //  double offset = parameters.at(0); 
  //  double slope  = parameters.at(0); 

  double offset_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double offset_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  double slope_lower = this -> GetParameter(1) -> GetLowerLimit(); 
  double slope_upper = this -> GetParameter(1) -> GetUpperLimit(); 

  double quad_lower =  this -> GetParameter(2) -> GetLowerLimit(); 
  double quad_upper =  this -> GetParameter(2) -> GetUpperLimit(); 

  double probability_offset = 1.0 / (offset_upper - offset_lower); 
  double probability_slope = 1.0 / (slope_upper - slope_lower); 
  double probability_quad = 1.0 / (quad_upper - quad_lower); 

  probability = probability_offset * probability_slope * probability_quad; 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol2::ConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters)
{

  double probability = 1.0; 

  // define data values 

  double x       = datapoint -> GetValue(0); 
  double y       = datapoint -> GetValue(1); 
  double sigma_y = datapoint -> GetValue(2); 

  // define parameters and limits 

  double offset = parameters.at(0); 
  double slope  = parameters.at(1); 
  double quad   = parameters.at(2); 

  // calculate probability for a single measurement 
  
  probability = TMath::Gaus(y * 100.0, (offset + x * slope + x * x * quad) * 100.0, sigma_y * 100.0, true); 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelPol2::PoissonProbability(int nentries, std::vector <double> parameters)
{

  double probability = 1.0; 

  return probability; 

}

// --------------------------------------------------------- 

