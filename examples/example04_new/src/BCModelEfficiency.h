#ifndef __BCMODELEFFICIENCY__H
#define __BCMODELEFFICIENCY__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelEfficiency : public BCModel 
{

 public: 

  // constructor 

  BCModelEfficiency(); 

  BCModelEfficiency(const char* name); 

  ~BCModelEfficiency()
    { ;};  

  // methods 

  double LogAPrioriProbability(std::vector <double> parameters); 
  
  double LogLikelihood(std::vector <double> parameters); 

  TH1D * hist_bestfit; 
  TH1D * hist_efficiency; 

}; 

// --------------------------------------------------------- 

#endif 
