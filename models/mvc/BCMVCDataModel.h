/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCDATAMODEL__H
#define __BCMVCDATAMODEL__H

#include "../../BAT/BCModel.h"

#include "BCMVCombination.h"

#include <TH1D.h>
#include <TMatrixT.h>
#include <TMatrixD.h>
#include <TVectorT.h>

// ---------------------------------------------------------
class BCMVCDataModel : public BCModel
{
 public:

  // Constructor
  BCMVCDataModel(BCMVCombination* mvc);

  // Destructor
  ~BCMVCDataModel();

  // setters

  // set the number of measurements
  void SetNMeasurements(int n, double min, double max);

  // set the vector of measurements
  void SetVectorMeasurements(TVectorD measurements)
  { fVectorMeasurements.Clear();
    fVectorMeasurements.ResizeTo(measurements);
    fVectorMeasurements=measurements; };

  // set parameters
  void SetParameters(std::vector<double> parameters);

  // set the ranges of the measurements for each measurement
  // individually
  void SetMeasurementRanges(const std::vector<double> & min, const std::vector<double> & max);

  // set the ranges of the measurements for all measurements
  void SetMeasurementRanges(double min, double max);

  // set vector of observables
  void SetVectorObservable(std::vector<int> vec)
  { fVectorObservable = vec; };

  // set the covariance matrix
  void SetCovarianceMatrix(TMatrixD matrix);

  // set the chi2 histogram
  void SetHistChi2(TH1D* hist)
  { fHistChi2 = hist; };

  // misc

  // print scatter plots of toy models and the chi2 indicating the
  // position of the observed data
  void PrintToys(std::string filename);

  // print a summary to the screen
  void PrintSummary();

  // calculate the chi2
  double Chi2(TVectorD observables, TVectorD measurements);

  // BAT methods

  // Methods to overload, see file BCMVCDataModel.cxx
  double LogAPrioriProbability(const std::vector<double> &parameters);

  double LogLikelihood(const std::vector<double> &parameters);

  void MCMCUserIterationInterface();

 private:

  // the vector of measurements
  TVectorD fVectorMeasurements;

  // the total covariance matrix
  TMatrixD fCovarianceMatrix;

  // the inverse of the covariance matrix
  TMatrixD fInvCovarianceMatrix;

  // the determinant of the covariance matrix
  double fDetCovariance;

  // the parameters
  std::vector<double> fPars;

  // the vector of observables
  TVectorD fVectorObservables;

  // the vector of the index of the observables being measured
  std::vector<int> fVectorObservable;

  // the chi2 histogram
  TH1D* fHistChi2;

};
// ---------------------------------------------------------

#endif

