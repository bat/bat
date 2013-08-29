#include "MyFit.h"

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MyFit::MyFit() : MVFit()
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
	AddParameter("alpha", 0.0, 100.);
	SetPriorConstant("alpha");
}

// ---------------------------------------------------------
void MyFit::DefineObservables()
{
	AddObservable("xs", 0., 100.);
}

// ---------------------------------------------------------
double MyFit::CalculateObservable(int index, const std::vector<double> &parameters)
{
	double alpha = parameters.at(0);

	if (index == 0) { // cross-section

		double xs = 0. + 0.1*(alpha-50) *(alpha-50);

		return xs;
	}
	else if (index == 1) {
		; // ...
	}

	BCLog::OutWarning("CalculateObservable. Index of observable out of range. Return 0");

	// if unknown, return 0;
	return 0;
}

// ---------------------------------------------------------
