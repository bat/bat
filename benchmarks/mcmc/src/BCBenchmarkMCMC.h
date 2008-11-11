#ifndef __BCBENCHMARKMCMC__H
#define __BCBENCHMARKMCMC__H

#include "BCModel.h" 

#include "TF1.h" 

// --------------------------------------------------------- 

class BCBenchmarkMCMC : public BCModel 
{

 public: 

  // constructor 

  BCBenchmarkMCMC(const char* name); 

  ~BCBenchmarkMCMC()
    { ;};  

  // methods 

  double LogAPrioriProbability(std::vector <double> parameters)
	{ return 0; }; 
  
  double LogLikelihood(std::vector <double> parameters); 

	void PerformTest(std::vector<double> parameters, 
									 int index, 
									 BCH1D * hist, 
									 double * chi2, 
									 bool flag_print = true, 
									 const char * filename = "test.ps"); 

	void SetTestFunction(TF1 * testfunction)
		{ fTestFunction = testfunction; }; 

	TF1 * fTestFunction; 

}; 

// --------------------------------------------------------- 

#endif 
