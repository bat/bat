/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#ifndef __BCMVCOMBINATION__H
#define __BCMVCOMBINATION__H

#include "../../BAT/BCModel.h"

#include <TMath.h>
#include <TMatrixDEigen.h>
#include <TMatrixT.h>
#include <TVectorT.h>

class BCMVCMeasurement;
class BCMVCUncertainty;
class BCMVCObservable;

// ---------------------------------------------------------
class BCMVCombination : public BCModel
{

 public:

  // Constructor
  BCMVCombination();

  // Destructor
  ~BCMVCombination();

  // Add a parameter
  // name: the name of the parameter
  // min:  the minimum value of the parameter
  // max:  the maximum value of the parameter
  void AddObservable(std::string name, double min, double max);

  // Add a source of systematic uncertainty
  // name: the name of the uncertainty
  void AddUncertainty(std::string name);

  // Add a measurement
  // name:          the name of the measurement
  // observable:    the name of the observable
  // central:       the central value of the measurement
  // uncertainties  the uncertainties
  void AddMeasurement(std::string name, std::string observable, double value, std::vector<double> uncertainties);

  // getters

  // return the number of observables
  int GetNObservables()
  { return fNObservables; };

  // return the number of uncertainties
  int GetNUncertainties()
  { return int(fUncertainties.size()); };

  // return the number of measurements
  int GetNMeasurements()
  { return int(fMeasurements.size()); };

  // return the number of measurements
  int GetNActiveMeasurements();

  // return a specific uncertainty
  BCMVCUncertainty* GetUncertainty(int index)
    { return fUncertainties.at(index); }

  // return a specific measurement
  BCMVCMeasurement* GetMeasurement(int index)
    { return fMeasurements.at(index); }

  // return the total covariance matrix
  TMatrixD GetCovarianceMatrix()
  { return fCovarianceMatrix; };

  // return the BLUE weights
  TMatrixD GetBLUEWeights()
  { return fBLUEWeights; };

  // return the BLUE central values
  TVectorD GetBLUECentralValues()
  { return fBLUECentral; };

  // return the BLUE uncertainties
  TVectorD GetBLUEUncertainties()
  { return fBLUEUncertainties; };

  // return the BLUE uncertainties for a certain source of
  // uncertainty
  // index: the index of the uncertainty source
  TVectorD GetBLUEUncertainties(int index)
  { return fBLUEUncertaintiesPerSource.at(index); };

  // return the BLUE covariance matrix
  TMatrixD GetBLUECovarianceMatrix()
  { return fBLUECovarianceMatrix; };

  // return the BLUE covariance matrix for a certain source of
  // uncertainty
  // index: the index of the uncertainty source
  TMatrixD GetBLUECovarianceMatrix(int index)
  { return fBLUECovarianceMatrices.at(index); };

  // return the BLUE correlation matrix for a certain source of
  // uncertainty
  // index: the index of the uncertainty source
  TMatrixD GetBLUECorrelationMatrix(int index)
  { return fBLUECorrelationMatrices.at(index); };

  // return the BLUE correlation matrix
  TMatrixD GetBLUECorrelationMatrix()
  { return fBLUECorrelationMatrix; };

  // return the vector of observables
  std::vector<int> GetVectorObservable()
    { return fVectorObservable; };

  // return the vector of measurements
  TVectorD GetVectorMeasurements()
  { return fVectorMeasurements; };

  // return the index of the measurements
  int GetIndexMeasurement(std::string measurement, std::string observable);

  // return the index of the uncertainty
  int GetIndexUncertainty(std::string name);

  // return the index of the observable
  int GetIndexObservable(std::string name);

  // misc

  // read input file
  int ReadInput(std::string filename);

  // calculate the correlation matrix for a particular uncertainty
  void CalculateCorrelationMatrix(int index);

  // calculate the total covariance matrix
  void CalculateCovarianceMatrix(std::vector<double> nuisance = std::vector<double>(0));

  // calculate helper vectors
  void CalculateHelperVectors();

  // check for positive definiteness of the covariance matrix
  bool PositiveDefinite();

  // calculate BLUE
  void CalculateBLUE();

  // calculate all necessary matrices
  void PrepareAnalysis();

  // output

  // print summary to screen
  void PrintBLUEResults(std::string filename);

  void PrintMatrix(TMatrixT<double> &matrix, std::string name="matrix");

  void PrintVector(TVectorD &vector, std::string name="vector");

  // BAT methods

  double LogLikelihood(const std::vector<double> &parameters);

 protected:

  struct NuisanceParameter {
    int index_uncertainty;
    int index_measurement1;
    int index_measurement2;
    int index_rhoparameter;
    double pre;
  };

  // the names of the uncertainties
  std::vector<BCMVCUncertainty*> fUncertainties;

  // the measurements
  std::vector<BCMVCMeasurement*> fMeasurements;

  // the total covariance matrix
  TMatrixD fCovarianceMatrix;

  // the inverse of the covariance matrix
  TMatrixD fInvCovarianceMatrix;

  // the determinant of the covariance matrix
  double fDetCovariance;

  // helper: the vector of measurements
  TVectorD fVectorMeasurements;

  // helper: the vector of active measurements
  TVectorD fVectorActiveMeasurements;

  // the vector of the index of the observables being measured
  std::vector<int> fVectorObservable;

  // the vector of the index of the observables being measured if active
  std::vector<int> fVectorActiveObservable;

  // the BLUE matrix
  TMatrixD fBLUEWeights;

  // the BLUE central values
  TVectorD fBLUECentral;

  // the BLUE uncertainties
  TVectorD fBLUEUncertainties;

  // the BLUE uncertainties per source
  std::vector<TVectorD> fBLUEUncertaintiesPerSource;

  // the BLUE covariance matrix
  TMatrixD fBLUECovarianceMatrix;

  // the BLUE covariance matrices for all uncertainties
  std::vector<TMatrixD> fBLUECovarianceMatrices;

  // the BLUE covariance matrices for all uncertainties
  std::vector<TMatrixD> fBLUECorrelationMatrices;

  // the BLUE correlation matrix
  TMatrixD fBLUECorrelationMatrix;

  // number of observables
  int fNObservables;

  // number of nuisance parameters for correlations
  int fNNuisanceCorrelation;

  // nuisance parameters
  std::vector<NuisanceParameter> fNuisanceCorrelation;

  // the observables
  std::vector<BCMVCObservable*> fObservables;

};
// ---------------------------------------------------------

#endif

