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

  void DefineParameters(); 

  double APrioriProbability(std::vector <double> parameters); 
  
  double Likelihood(std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
