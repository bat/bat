/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCUNCERTAINTY__H
#define __BCMVCUNCERTAINTY__H

#include <string>
#include <vector>
#include <iostream>

#include <TMatrixT.h>
#include <TMatrixD.h>

// ---------------------------------------------------------
class BCMVCUncertainty
{
 public:

  // constructor
  // name: the name of the uncertainty
  BCMVCUncertainty(std::string name);

  // destructor
  ~BCMVCUncertainty();

  // getters

  // return the name of the uncertainty
  std::string GetName()
    { return fName; };

  // return the correlation matrix
  TMatrixD GetCorrelationMatrix()
  { return fCorrelationMatrix; };

  // return the covariance matrix
  TMatrixD GetCovarianceMatrix()
  { return fCovarianceMatrix; };

  // return the inverse of the covariance matrix
  TMatrixD GetInvCovarianceMatrix()
  { return fInvCovarianceMatrix; };

  // return the flag if the uncertainty is active or not
  bool GetFlagActive()
  { return fFlagActive; };

  // setters
  void SetCorrelationMatrix(const TMatrixD &matrix);

  // setters
  void SetCovarianceMatrix(const TMatrixT<double> &matrix);

  // set flag if uncertainty is active for the combination
  void SetFlagActive(bool flag)
  { fFlagActive = flag; };

 private:

  // the name of the uncertainty
  std::string fName;

  // the symmetric correlation matrix; the number of columns and rows
  // is equal to the number of measurements
  TMatrixD fCorrelationMatrix;

  // the symmetric covariance matrix; the number of columns and rows
  // is equal to the number of measurements
  TMatrixD fCovarianceMatrix;

  // the inverse of the covariance matrix
  TMatrixD fInvCovarianceMatrix;

  // flag: active in combination (true) or not (false)
  bool fFlagActive;
};
// ---------------------------------------------------------

#endif
