// --------------------------------------------------------- 
//
// This class derives from BCModel. It describes a spectrum with a
// flat (background) and Gaussian (signal) contribution. Each bin
// content fluctuates with a Poissonian distribution. Two parameters,
// the number of signal and background events, are defined in the
// model.
// 
// --------------------------------------------------------- 

#ifndef __BCMODELSIGNAL__H
#define __BCMODELSIGNAL__H

#include "BCModel.h" 

// --------------------------------------------------------- 

class BCModelSignal : public BCModel 
{

 public: 

	// constructor 

	BCModelSignal(); 

	BCModelSignal(const char* name); 

	// destructor 

	~BCModelSignal()
		{ ;}; 

	// methods 

	void DefineParameters(); 

	virtual double LogAPrioriProbability(std::vector <double> parameters); 

	virtual double LogLikelihood(std::vector <double> parameters); 

	double FitFunction(std::vector <double> x, std::vector <double> parameters); 

 private: 

	double integral_f_B(double E, double dE, int nbins); 

	double integral_f_S(double E, double dE, int nbins); 

}; 

// --------------------------------------------------------- 

#endif 
