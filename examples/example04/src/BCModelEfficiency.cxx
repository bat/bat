#include "BCModelEfficiency.h" 

#include <TMath.h> 

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency() : BCModel()
{

  // define parameters 

  this -> DefineParameters(); 

}

// --------------------------------------------------------- 

BCModelEfficiency::BCModelEfficiency(const char* name) : BCModel(name)
{

  // define parameters 

  this -> DefineParameters(); 

}

// --------------------------------------------------------- 

void BCModelEfficiency::DefineParameters()
{
  this -> AddParameter("epsilon", 0.0,  1.0); // index 0 
  this -> AddParameter("lambda",  0.0, 20.0); // index 1 

}

// --------------------------------------------------------- 

double BCModelEfficiency::APrioriProbability(std::vector <double> parameters)
{

  // get parameter ranges 

  double epsilon_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double epsilon_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  double lambda_lower = this -> GetParameter(1) -> GetLowerLimit(); 
  double lambda_upper = this -> GetParameter(1) -> GetUpperLimit(); 

  // calculate probabilities 

  double probability_epsilon = 1.0 / (epsilon_upper - epsilon_lower); 
  double probability_lambda = 1.0 / (lambda_upper - lambda_lower);   

  double probability = probability_epsilon * probability_lambda; 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelEfficiency::Likelihood(std::vector <double> parameters)
{

  // get parameters 

  double epsilon = parameters.at(0); 
  double lambda  = parameters.at(1); 

  // get data values 

  int x = TMath::Nint(this -> GetDataPoint(0) -> GetValue(0)); 
  int y = TMath::Nint(this -> GetDataPoint(0) -> GetValue(1)); 
      
  // calculate probability assuming a Gaussian distribution for each point 
  
  double poisson  = TMath::Poisson(double(y), lambda); 
  double binomial = 1.0; 

  if (y < 20) 
    binomial = TMath::Binomial(y, x) * TMath::Power(epsilon, x) * TMath::Power(1.0 - epsilon, y - x);       

  double probability = poisson * binomial; 

  return probability; 

}

// --------------------------------------------------------- 
