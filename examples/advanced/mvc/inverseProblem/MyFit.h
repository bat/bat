#ifndef __MYFIT__H
#define __MYFIT__H

#include <BAT/BCMVCPhysicsModel.h>

// ---------------------------------------------------------
class MyFit : public BCMVCPhysicsModel
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

