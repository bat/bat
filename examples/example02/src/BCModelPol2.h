#ifndef __BCMODELPOL2__H
#define __BCMODELPOL2__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelPol2 : public BCModel 
{

 public: 

  // constructor 

  BCModelPol2(); 

  BCModelPol2(const char* name); 

  ~BCModelPol2() 
    { ;}; 

  // methods 

  virtual void DefineParameters(); 

  virtual double LogAPrioriProbability(std::vector <double> parameters); 

  virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

  virtual double LogPoissonProbability(int nentries, std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
