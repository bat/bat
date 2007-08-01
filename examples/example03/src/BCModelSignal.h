#ifndef __BCMODELSIGNAL__H
#define __BCMODELSIGNAL__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelSignal : public BCModel 
{

 public: 

  // constructor 

  BCModelSignal(); 

  BCModelSignal(const char* name); 

  // destructor 

  ~BCModelSignal()
    { ;}; 

  // methods 

  void DefineParameters(); 

  virtual double APrioriProbability(std::vector <double> parameters); 

  virtual double Likelihood(std::vector <double> parameters); 
  
  virtual double SamplingFunction(std::vector <double> parameters); 

 private: 

  double integral_f_B(double E, double dE, int nbins); 

  double integral_f_S(double E, double dE, int nbins); 

}; 

// --------------------------------------------------------- 

#endif 
