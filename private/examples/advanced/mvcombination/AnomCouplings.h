#ifndef __BAT__ANOMCOUPLINGS__H
#define __BAT__ANOMCOUPLINGS__H

#include "MVFit.h"

// ---------------------------------------------------------
class AnomCouplings : public MVFit
{
 public:

  // Constructor
  AnomCouplings();

  // Destructor
  ~AnomCouplings();
	
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

