// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a quadratic
// correlation relation between measured points. Three parameters are
// defined within the model, an offset, a slope and a quadtratic
// term. The data are points (x,y) with an uncertainty on y, s_y. The
// uncertainty is assumed to be Gaussian.
// 
// --------------------------------------------------------- 

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

	// destructor 

	~BCModelPol2() 
		{ ;}; 

	// methods 

	virtual double LogAPrioriProbability(std::vector <double> parameters); 

	virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

	virtual double LogPoissonProbability(int nentries, std::vector <double> parameters); 

	double FitFunction(std::vector <double> x, std::vector <double> parameters); 	

}; 

// --------------------------------------------------------- 

#endif 
