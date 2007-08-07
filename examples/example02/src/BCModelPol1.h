#ifndef __BCMODELPOL1__H
#define __BCMODELPOL1__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelPol1 : public BCModel 
{

 public: 

  // constructor 

  BCModelPol1(); 

  BCModelPol1(const char* name); 

  ~BCModelPol1()
    { ;};  

  // methods 

  virtual void DefineParameters(); 

  virtual double LogAPrioriProbability(std::vector <double> parameters); 
  
  virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

  virtual double LogPoissonProbability(int nentries, std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
