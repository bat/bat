// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a constant
// correlation relation between measured points. One parameter is
// defined within the model, a constant. The data are points (x,y)
// with an uncertainty on y, s_y. The uncertainty is assumed to be
// Gaussian.
// 
// --------------------------------------------------------- 

#ifndef __BCMODELPOL0__H
#define __BCMODELPOL0__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelPol0 : public BCModel 
{

 public: 

	// constructors 

	BCModelPol0(); 

	BCModelPol0(const char* name); 

	// destructor 

	~BCModelPol0() 
		{ ;}; 

	// methods 

	virtual double LogAPrioriProbability(std::vector <double> parameters); 

	virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

	virtual double LogPoissonProbability(int nentries, std::vector <double> parameters); 

	double FitFunction(std::vector <double> x, std::vector <double> parameters); 	

}; 

// --------------------------------------------------------- 

#endif 
