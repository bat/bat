#include "BCModelBackground.h" 

#include <TMath.h> 

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

  this -> AddParameter("background", 0.0, 40.0);  // index 0

}

// --------------------------------------------------------- 

double BCModelBackground::APrioriProbability(std::vector <double> parameters)
{

  double probability = 1.0; 

  //  double constant = parameters.at(0); 

  double background_lower = this -> GetParameter(0) -> GetLowerLimit(); 
  double background_upper = this -> GetParameter(0) -> GetUpperLimit(); 

  // calculate probability 

  probability = 1.0 / (background_upper - background_lower); 

  return probability; 

}

// --------------------------------------------------------- 

double BCModelBackground::Likelihood(std::vector <double> parameters)
{

  double logprobability = 0.0; 

  // get parameter values 

  double background = parameters.at(0); 

  // loop over all data points 

  int nbins = this -> GetNDataPoints() - 1; 

  // get histogram parameters 

  double Emin = this -> GetDataPoint(0) -> GetValue(0); 
  double Emax = this -> GetDataPoint(0) -> GetValue(1); 
  double dE   = (Emax - Emin) / double(nbins); 

  // debug
  //  cout << endl; 
  //  cout << background << " " << nbins << " " << dE << endl; 

  for (int ipoint = 1; ipoint < nbins + 1; ipoint++)
    {      
      // define data values 

      double energy = this -> GetDataPoint(ipoint) -> GetValue(0); 
      double events = this -> GetDataPoint(ipoint) -> GetValue(1); 

      // calculate probability 

      double expected = background * this -> integral_f_B(energy, dE, nbins) * dE; 

      logprobability += TMath::Log10(TMath::PoissonI(events, expected));       

      // debug
      //      cout << energy << " " << events << " " << expected << " " << TMath::PoissonI(events, expected) << endl;       
    }

  // debug 
  //  cout << " prob : " << logprobability << " " << TMath::Power(10.0, logprobability) << endl; 


  return TMath::Power(10.0, logprobability); 

}

// --------------------------------------------------------- 

double BCModelBackground::integral_f_B(double E, double dE, int nbins)
{

  return 1.0 / double(nbins) / dE; 

}

// --------------------------------------------------------- 

double BCModelBackground::SamplingFunction(std::vector <double> parameters)
{

  return 1.0; 

}

// --------------------------------------------------------- 

