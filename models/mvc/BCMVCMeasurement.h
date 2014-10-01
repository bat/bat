/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCMEASUREMENT__H
#define __BCMVCMEASUREMENT__H

#include <string>
#include <vector>

// ---------------------------------------------------------
class BCMVCMeasurement
{
 public:

  // constructor
  // name: the name of the measurement
  BCMVCMeasurement(std::string name);

  // destructor
  ~BCMVCMeasurement();

  // getters

  // return the name of the observable
  std::string GetName()
    { return fName; };

  // return the index of the observable
  int GetObservable()
  { return fObservable; };

  // return the central value
  double GetCentralValue()
  { return fCentralValue; };

  // return the set of uncertainties
  std::vector<double> GetUncertainties()
    { return fUncertainties; };

  // return a single uncertainty
  double GetUncertainty(int index)
  { return fUncertainties.at(index); };

  // return the total uncertainty
  double GetTotalUncertainty();

  // return the flag if the measurement is active or not
  bool GetFlagActive()
  { return fFlagActive; };

  // setters

  // set the (index of the) observable this measurement corresponds to
  void SetObservable(int index)
  { fObservable = index; };

  // set the central value of the measurement
  void SetCentralValue(double value)
  { fCentralValue = value; };

  // set the uncertainties on the measurement
  void SetUncertainties(std::vector<double> uncertainties)
  { fUncertainties = uncertainties; };

  // set flag if measurement is active for the combination
  void SetFlagActive(bool flag)
  { fFlagActive = flag; };

 private:

  // the name of the measurement
  std::string fName;

  // the (index of the) observale this measurement corresponds to
  int fObservable;

  // the central value of the measurement
  double fCentralValue;

  // the uncertainties on the measurements
  std::vector<double> fUncertainties;

  // flag: active in combination (true) or not (false)
  bool fFlagActive;

};
// ---------------------------------------------------------

#endif

