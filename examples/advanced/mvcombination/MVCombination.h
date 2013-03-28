#ifndef __BAT__MVCOMBINATION__H
#define __BAT__MVCOMBINATION__H

#include <BAT/BCModel.h>

#include "MVMeasurement.h"
#include "MVUncertainty.h"

#include <TMatrixT.h>
#include <TVectorT.h>

// This is a MVCombination header file.
// Model source code is located in file MVCombination/MVCombination.cxx

// ---------------------------------------------------------
class MVCombination : public BCModel
{
 public:

  // Constructor
  MVCombination();

  // Destructor
  ~MVCombination();

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
  { return GetNParameters(); };

  // return the number of uncertainties
  int GetNUncertainties() 
  { return int(fUncertainties.size()); };

  // return the number of measurements
  int GetNMeasurements() 
  { return int(fMeasurements.size()); };

  // return a specific uncertainty
  MVUncertainty* GetUncertainty(int index) 
    { return fUncertainties.at(index); }

  // return a specific measurement
  MVMeasurement* GetMeasurement(int index) 
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

  // misc

  // read input file
  int ReadInput(std::string filename);
			
  // calculate the correlation matrix for a particular uncertainty
  void CalculateCorrelationMatrix(int index);

  // calculate the total covariance matrix
  void CalculateCovarianceMatrix();
			
  // calculate BLUE
  void CalculateBLUE();

  // output

  // print summary to screen
  void PrintBLUEResults(std::string filename);

  void PrintMatrix(TMatrixT<double> &matrix, std::string name="matrix");

  void PrintVector(TVectorD &vector, std::string name="vector");

  // BAT methods

  // Methods to overload, see file MVCombination.cxx
  double LogAPrioriProbability(const std::vector<double> &parameters);

  double LogLikelihood(const std::vector<double> &parameters);

 private:

  int GetIndexObservable(std::string name); 

  // the names of the uncertainties
  std::vector<MVUncertainty*> fUncertainties;

  // the measurements
  std::vector<MVMeasurement*> fMeasurements;

  // the total covariance matrix
  TMatrixD fCovarianceMatrix;

  // the inverse of the covariance matrix
  TMatrixD fInvCovarianceMatrix;

  // the determinant of the covariance matrix
  double fDetCovariance;

  // the vector of measurements
  TVectorD fVectorMeasurements;

  // the vector of the index of the observables being measured
  std::vector<int> fVectorObservable;

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

};
// ---------------------------------------------------------

#endif

