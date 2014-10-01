/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCPHYSICSMODEL__H
#define __BCMVCPHYSICSMODEL__H

#include "BCMVCombination.h"

// ---------------------------------------------------------
class BCMVCPhysicsModel : public BCMVCombination
{
 public:

  // Constructor
  BCMVCPhysicsModel();

  // Destructor
  ~BCMVCPhysicsModel();

  // Add a parameter
  // name: the name of the parameter
  // min:  the minimum value of the parameter
  // max:  the maximum value of the parameter
  void AddObservable(std::string name, double min, double max);

  // return a value for an observable
  // index: the index of the variable
  // parameters: the physics parameters
  virtual double CalculateObservable(int index, const std::vector<double> &parameters)
  { (void) index; // suppress compiler warning about unused parameters
    (void) parameters; // suppress compiler warning about unused parameters
    return 0; };

  // BAT methods

  // the log of the likelihood
  double LogLikelihood(const std::vector<double> &parameters);

 private:

};
// ---------------------------------------------------------

#endif

