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

  void DefineParameters(); 

  double LogAPrioriProbability(std::vector <double> parameters); 
  
  double LogLikelihood(std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
