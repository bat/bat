#include "MyFit.h"

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MyFit::MyFit() : BCMVCPhysicsModel()
{
	DefineParameters();
	DefineObservables();
}

// ---------------------------------------------------------
MyFit::~MyFit()
{
}

// ---------------------------------------------------------
void MyFit::DefineParameters()
{
	// add parameters here
	AddParameter("#kappa_{1}", -2.5, 2.5);
	AddParameter("#kappa_{2}", -2.5, 2.5);
	SetPriorConstant("#kappa_{1}");
	SetPriorConstant("#kappa_{2}");
}

// ---------------------------------------------------------
void MyFit::DefineObservables()
{
	AddObservable("br_e",  0.,   1.);
	AddObservable("br_mu", 0.,   1.);
	AddObservable("xs",    0., 100.);
}

// ---------------------------------------------------------
double MyFit::CalculateObservable(int index, const std::vector<double> &parameters)
{
	double kappa1 = parameters.at(0);
	double kappa2 = parameters.at(1);

	if (index == 0) { // br e
		double br = 0.60 + 0.1*kappa1; 

		return br;
	}
	else if (index == 1) { // br mu
		double br = 0.40 - 0.1*kappa1; 

		return br;
	}
	else if (index == 2) { // cross-section

		double xs = 40 + 10*(kappa1-0.6)*(kappa1-0.6) + 20*kappa2*kappa2;

		return xs;
	}

	BCLog::OutWarning("CalculateObservable. Index of observable out of range. Return 0");

	// if unknown, return 0;
	return 0;
}

// ---------------------------------------------------------
