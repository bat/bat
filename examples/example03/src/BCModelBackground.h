// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a flat background
// spectrum with Poissonian fluctuations in each bin. Only one
// parameter, the number of background events, is defined in the
// model. 
// 
// --------------------------------------------------------- 

#ifndef __BCMODELBACKGROUND__H
#define __BCMODELBACKGROUND__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelBackground : public BCModel 
{

 public: 

	// constructor 

	BCModelBackground(); 

	BCModelBackground(const char* name); 

	// destructor 

	~BCModelBackground()
		{ ;}; 

	// methods 

	void DefineParameters(); 

	virtual double LogAPrioriProbability(std::vector <double> parameters); 

	virtual double LogLikelihood(std::vector <double> parameters); 

 private: 

	double integral_f_B(double E, double dE, int nbins); 

}; 

// --------------------------------------------------------- 

#endif 
