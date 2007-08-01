#ifndef __BCMODELBACKGROUND__H
#define __BCMODELBACKGROUND__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelBackground : public BCModel 
{

 public: 

  // constructor 

  BCModelBackground(); 

  BCModelBackground(const char* name); 

  // destructor 

  ~BCModelBackground()
    { ;}; 

  // methods 

  void DefineParameters(); 

  virtual double APrioriProbability(std::vector <double> parameters); 

  virtual double Likelihood(std::vector <double> parameters); 
  
  virtual double SamplingFunction(std::vector <double> parameters); 

 private: 

  double integral_f_B(double E, double dE, int nbins); 

}; 

// --------------------------------------------------------- 

#endif 
