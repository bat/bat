#ifndef __BCMODELPOL0__H
#define __BCMODELPOL0__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelPol0 : public BCModel 
{

 public: 

  // constructor 

  BCModelPol0(); 

  BCModelPol0(const char* name); 

  ~BCModelPol0() 
    { ;}; 

  // methods 

  virtual void DefineParameters(); 

  virtual double APrioriProbability(std::vector <double> parameters); 
  
  virtual double ConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

  virtual double PoissonProbability(int nentries, std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
