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

	// constructors 

	BCModelPol1(); 

	BCModelPol1(const char* name); 

	// destructor 

	~BCModelPol1()
		{ ;};  

	// methods 

	void DefineParameters(); 

	double LogAPrioriProbability(std::vector <double> parameters); 

	double LogLikelihood(std::vector <double> parameters); 

}; 

// --------------------------------------------------------- 

#endif 
