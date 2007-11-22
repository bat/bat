// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a linear correlation
// relation between measured points. Two parameters are defined within
// the model, an offset, b, and a slope, m. The data are points (x,y)
// with an uncertainty on y, s_y. The uncertainty is assumed to be
// Gaussian. 
// 
// --------------------------------------------------------- 

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

	// destructor 

	~BCModelPol1()
		{ ;};  

	// methods 

	virtual void DefineParameters(); 

	virtual double LogAPrioriProbability(std::vector <double> parameters); 

	virtual double LogConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

	virtual double LogPoissonProbability(int nentries, std::vector <double> parameters); 

	virtual void GetRandomVectorMetro(std::vector <double> &x); 

	double FitFunction(std::vector <double> x, std::vector <double> parameters); 	

}; 

// --------------------------------------------------------- 

#endif 
