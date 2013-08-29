#ifndef __BAT__MYFIT__H
#define __BAT__MYFIT__H

#include "MVFit.h"

// ---------------------------------------------------------
class MyFit : public MVFit
{
 public:

  // Constructor
  MyFit();

  // Destructor
  ~MyFit();
	
	// define the physics parameters
	virtual void DefineParameters();

	// define the observables
	virtual void DefineObservables();

	// return a value for an observable
	// index: the index of the variable
	// parameters: the physics parameters
  double CalculateObservable(int index, const std::vector<double> &parameters);

 private:

};
// ---------------------------------------------------------

#endif

