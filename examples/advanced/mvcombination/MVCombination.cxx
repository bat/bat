#include "MVCombination.h"

#include <iomanip>
#include <iostream>
#include <fstream>

#include <TMath.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

// ---------------------------------------------------------
MVCombination::MVCombination() : BCModel("MVCombination")
{
}

// ---------------------------------------------------------
MVCombination::~MVCombination()
{
}

// ---------------------------------------------------------
void MVCombination::AddObservable(std::string name, double min, double max)
{
  AddParameter(name.c_str(), min, max); 
}

// ---------------------------------------------------------
void MVCombination::AddUncertainty(std::string name)
{
  MVUncertainty* u = new MVUncertainty(name);

  fUncertainties.push_back(u);
}

// ---------------------------------------------------------
void MVCombination::AddMeasurement(std::string name, std::string observable, double value, std::vector<double> uncertainties)
{
  // get index of the corresponding observable
  int index = GetIndexObservable(observable);
	
  // check if observable exists
  if (index < 0) {
    BCLog::OutWarning("Observable does not exist. Measurement was not added.");
    return;
  }

  MVMeasurement* m = new MVMeasurement(name);
  m->SetObservable(index);
  m->SetCentralValue(value);
  m->SetUncertainties(uncertainties);

  fMeasurements.push_back(m);

  fVectorObservable.push_back(index);

  int n = GetNMeasurements();
  fVectorMeasurements.ResizeTo(n);
  fVectorMeasurements[n-1]=value;
}

// ---------------------------------------------------------
double MVCombination::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob = 0.;

  int nmeas = GetNMeasurements();

  // copy parameters into a vector
  TVectorD observables(nmeas);

  for (int i = 0; i < nmeas; ++i) {
    observables[i] = parameters[fVectorObservable[i]];
  }

  TVectorD prod1 = observables - fVectorMeasurements;
  TVectorD prod2 = fInvCovarianceMatrix * prod1;
  double prod = prod1 * prod2;

  logprob = -0.5 * prod - log(TMath::Power(2*TMath::Pi(), nmeas/2.) * sqrt(fDetCovariance));

  return logprob;
}

// ---------------------------------------------------------
double MVCombination::LogAPrioriProbability(const std::vector<double> &parameters)
{
  double logprob = 0.;

  return logprob;
}

// ---------------------------------------------------------
int MVCombination::ReadInput(std::string filename)
{
  // open input file
  ifstream infile;
  infile.open(filename.c_str(), std::ifstream::in);
  
  // check if file is open
  if (!infile.is_open()) {
    BCLog::OutWarning(Form("MVCombination::ReadInput. Could not open input file %s.", filename.c_str()));
    return 0;
  }

  int nobservables;
  int nmeasurements;
  int nuncertainties;
	 
  infile >> nobservables >> nmeasurements >> nuncertainties;

  std::vector<std::string> observable_names;
	 
  for (int i = 0; i < nobservables; ++i) {
    std::string name;
    double min;
    double max;
    infile >> name >> min >> max;
		 
    // add observable
    AddObservable(name.c_str(), min, max);
  }

  for (int i = 0; i < nuncertainties; ++i) {
    std::string name;
    infile >> name;
		 
    // add uncertainty
    AddUncertainty(name);
  }

  for (int i = 0; i < nmeasurements; ++i) {
    std::string name;
    std::string observable;
    double central;
    std::vector<double> uncertainties(0);
		 
    infile >> name;
    infile >> observable;
    infile >> central;

    for (int j = 0; j < nuncertainties; ++j) {
      double uncertainty;
      infile >> uncertainty;
      uncertainties.push_back(uncertainty);
    }

    // add measurement
    AddMeasurement(name, observable, central, uncertainties);
  }

  for (int i = 0; i < nuncertainties; ++i) {
    TMatrixD mat(nmeasurements, nmeasurements);

    for (int j = 0; j < nmeasurements; ++j) 
      for (int k = 0; k < nmeasurements; ++k) {
	double corr;
	infile >> corr;
	mat[j][k] = corr;
      }

    // set correlation matrix
    GetUncertainty(i)->SetCorrelationMatrix(mat);
  }

  // close input file
  infile.close();

  // prepare analysis
  for (int i = 0; i < nuncertainties; ++i)
    CalculateCorrelationMatrix(i);

  CalculateCovarianceMatrix();

  // calculate BLUE
  CalculateBLUE();

  // no error
  return 1;
}

// ---------------------------------------------------------
void MVCombination::CalculateCorrelationMatrix(int index)
{
  MVUncertainty* u = fUncertainties.at(index);

  int n = GetNMeasurements();

  TMatrixD cov(n, n);

  TMatrixD corr = u->GetCorrelationMatrix();

  for (int i = 0; i < n; ++i) 
    for (int j = 0; j < n; ++j) {
      double sigma_i = GetMeasurement(i)->GetUncertainty(index);
      double sigma_j = GetMeasurement(j)->GetUncertainty(index);
			
      double corr_ij = corr[i][j];
      cov[i][j] = corr_ij*sigma_i*sigma_j;
    }
  u->SetCovarianceMatrix(cov);
}

// ---------------------------------------------------------
void MVCombination::CalculateCovarianceMatrix()
{
  int n = GetNUncertainties();

  fCovarianceMatrix.Clear();
  fCovarianceMatrix.ResizeTo(GetNMeasurements(), GetNMeasurements());
  for (int i = 0; i < n; ++i) {
    MVUncertainty* u = fUncertainties.at(i);
    TMatrixD cov = u->GetCovarianceMatrix();
    fCovarianceMatrix += cov;
  }
  fInvCovarianceMatrix.Clear();
  fInvCovarianceMatrix.ResizeTo(fCovarianceMatrix);
  fInvCovarianceMatrix = fCovarianceMatrix;
  fInvCovarianceMatrix.Invert();

  fDetCovariance = fCovarianceMatrix.Determinant();
}

// ---------------------------------------------------------
void MVCombination::CalculateBLUE()
{
  int nmeas = GetNMeasurements();
  int nobs = GetNObservables();
  int nunc = GetNUncertainties();

  // calculate U matrix
  TMatrixD u(nmeas, nobs);

  for (int i = 0; i < nmeas; ++i)
    for (int j = 0; j < nobs; ++j) {
      MVMeasurement* m = GetMeasurement(i);
      if (m->GetObservable() == j)
	u[i][j] = 1;
      else 
	u[i][j] = 0;
    }

  // calculate weight matrix
  TMatrixD ut = u;
  ut.Transpose(ut); 

  TMatrixD m1 = ut * fInvCovarianceMatrix;
  TMatrixD m2 = m1 * u;
  m2.Invert();

  fBLUEWeights.Clear();
  fBLUEWeights.ResizeTo(nobs, nmeas);
	
  fBLUEWeights = m2*m1;

  // calculate central values
  fBLUECentral.Clear();
  fBLUECentral.ResizeTo(nobs);

  fBLUECentral = fBLUEWeights * fVectorMeasurements;

  // calculate covariance matrix
  fBLUECovarianceMatrix.Clear();
  fBLUECovarianceMatrix.ResizeTo(nobs, nobs);

  TMatrixD weightt = fBLUEWeights;
  weightt.Transpose(weightt); 

  fBLUECovarianceMatrix = fBLUEWeights * fCovarianceMatrix * weightt;

  // calculate uncertainties, covariance and correlation matrices for each uncertainty
  for (int i = 0; i < nunc; ++i) {

    // calculate covariance matrix
    MVUncertainty* u = GetUncertainty(i);
    TMatrixD cov = u->GetCovarianceMatrix();
    TMatrixD mat = fBLUEWeights * cov * weightt;
    fBLUECovarianceMatrices.push_back(mat);

    // calculate uncertainties
    TVectorD vec(nobs);
    for (int j = 0; j < nobs; ++j) {
      double sigma_j = sqrt(mat[j][j]);
      vec[j] = sigma_j;
    }
    fBLUEUncertaintiesPerSource.push_back(vec);

    TMatrixD mat2(nobs, nobs);
    for (int j = 0; j < nobs; ++j) 
      for (int k = 0; k < nobs; ++k) {
	mat2[j][k] = mat[j][k]/vec[j]/vec[k];
      }
    fBLUECorrelationMatrices.push_back(mat2);
  }

  // calculate uncertainties
  fBLUEUncertainties.Clear();
  fBLUEUncertainties.ResizeTo(nobs);

  for (int i = 0; i < nobs; ++i)
    fBLUEUncertainties[i] = sqrt(fBLUECovarianceMatrix[i][i]);
	
  // calculate correlation matrix
  fBLUECorrelationMatrix.Clear();
  fBLUECorrelationMatrix.ResizeTo(nobs, nobs);

  for (int i = 0; i < nobs; ++i)
    for (int j = 0; j < nobs; ++j) 
      fBLUECorrelationMatrix[i][j] = fBLUECovarianceMatrix[i][j] / fBLUEUncertainties[i] / fBLUEUncertainties[j];
}

// ---------------------------------------------------------
void MVCombination::PrintBLUEResults(std::string filename)
{
  // open file
  std::ofstream ofi(filename.c_str());
  
  // check if file is open
  if (!ofi.is_open()) {
    std::cerr << "Couldn't open file " << filename << std::endl;
    return;
  }
  
  int nobs = GetNObservables();
  int nmeas = GetNMeasurements();
  int nunc = GetNUncertainties();

  ofi << std::endl;
  ofi << "Summary of the combination" << std::endl;
  ofi << "==========================" << std::endl << std::endl;

  ofi << "* Observables:" << std::endl;
  ofi << "  Observable (range): " << std::endl;
  for (int i = 0; i < nobs; ++i) 
    ofi << "  " << std::setiosflags(std::ios::left) << GetParameter(i)->GetName()
	      << " (" << GetParameter(i)->GetLowerLimit() << " - " << GetParameter(i)->GetUpperLimit() << ")" << std::endl;
  ofi << std::endl;

  ofi << "* Measurements:" << std::endl;
  ofi << "  Measurement (observable): central value +- total uncertainty" << std::endl;
  for (int i = 0; i < nmeas; ++i) {
    MVMeasurement* m = GetMeasurement(i);
    ofi << "  " << std::setiosflags(std::ios::left) << m->GetName() 
	      << std::setiosflags(std::ios::left) << " (" << GetParameter(m->GetObservable())->GetName() << ")" 
	      << ": " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << m->GetCentralValue()
	      << " +- " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << m->GetTotalUncertainty() << std::endl;
  }
  ofi << std::endl;
	
  ofi << "* Uncertainties:" << std::endl;
  ofi << "  Measurement (observable): Uncertainty (";
  for (int j = 0; j < nunc-1; ++j )
    ofi << GetUncertainty(j)->GetName() << ", ";
  ofi << GetUncertainty(nunc-1)->GetName() << ")" << std::endl;
  for (int i = 0; i < nmeas; ++i) {
    MVMeasurement* m = GetMeasurement(i);
    ofi << "  " << std::setiosflags(std::ios::left) << m->GetName() 
	      << std::setiosflags(std::ios::left) << " (" << GetParameter(m->GetObservable())->GetName() << "): ";
    for (int j = 0; j < nunc; ++j )
      ofi << std::setiosflags(std::ios::left) << std::setw(7) << m->GetUncertainty(j);
    ofi << std::endl;
  }
  ofi << std::endl;

  for (int i = 0; i < nunc; ++i) {
    MVUncertainty* u = GetUncertainty(i);
    ofi << "  " << u->GetName() << " " << "(correlation matrix)" << std::endl;
    TMatrixD mat = u->GetCorrelationMatrix();

    for (int j = 0; j < nmeas; ++j) {
      ofi << "  ";
      for (int k = 0; k < nmeas; ++k) 
	ofi << std::setw(7) << std::showpos << mat[j][k] << " ";
      ofi << std::noshowpos << std::endl;
    }
    ofi << std::endl;
  }

  ofi << "* BLUE results: " << std::endl;
  ofi << "  Observable: estimate +- total uncertainty" << std::endl;
  for (int i = 0; i < nobs; ++i) {
    if (i < fBLUECentral.GetNoElements())
      ofi << "  " << GetParameter(i)->GetName() << ": " << fBLUECentral[i] << " +- " << std::setprecision(4) << fBLUEUncertainties[i] << std::endl;
  }
  ofi << std::endl;

  ofi << "  Observable: Uncertainty (";
  for (int j = 0; j < nunc-1; ++j )
    ofi << GetUncertainty(j)->GetName() << ", ";
  ofi << GetUncertainty(nunc-1)->GetName() << ")" << std::endl;
  for (int i = 0; i < nobs; ++i) {
    ofi << "  " << std::setiosflags(std::ios::left) << GetParameter(i)->GetName()<< ":";
    for (int j = 0; j < nunc; ++j )
      ofi << " " << std::setiosflags(std::ios::left) << std::setw(7) << std::setprecision(4) << GetBLUEUncertainties(j)[i];
    ofi << std::endl;
  }
  ofi << std::endl;

  if (nobs > 1) {
    ofi << "  Individual correlation matrices " << std::endl;
    for (int i = 0; i < nunc; ++i) {
      TMatrixD mat = GetBLUECorrelationMatrix(i);
      ofi << "  " << GetUncertainty(i)->GetName() << std::endl;
      for (int j = 0; j < nobs; ++j) {
	ofi << "  ";
	for (int k = 0; k < nobs; ++k) {
	  ofi << std::setw(7) << std::setprecision(4) << std::showpos << mat[j][k] << " ";
	}
	ofi << std::noshowpos << std::endl;
      }
      ofi << std::endl;
    }
  }
  
  if (nobs > 1) { 
    ofi << "  Overall correlation matrix" << std::endl;
    TMatrixD mat = fBLUECorrelationMatrix;
    for (int j = 0; j < nobs; ++j) {
      ofi << "  ";
      for (int k = 0; k < nobs; ++k) 
	ofi << std::setw(7) << std::setprecision(4) << std::showpos << mat[j][k] << " ";
      ofi << std::noshowpos << std::endl;
    }
    ofi << std::endl;
  }      
  
  ofi << "  Weights [%]:" <<std::endl;
  for (int j = 0; j < nmeas; ++j) {
    ofi << "  " << GetMeasurement(j)->GetName() << " : ";
    for (int k = 0; k < nobs; ++k) 
      ofi << std::setw(7) << std::setprecision(4) << std::showpos << fBLUEWeights[k][j]*100. << " ";
    ofi << std::endl;
  }
  ofi << std::endl;

  // close file
  ofi.close();
}

// ---------------------------------------------------------
int MVCombination::GetIndexObservable(std::string name)
{
  int n = GetNObservables();

  // go through the list of parameters and compare strings
  for (int i = 0; i < n; ++i) {
    if (name == std::string(GetParameter(i)->GetName()))
      return i;
  }
	
  // return -1 if not in the list
  return -1;
}

// ---------------------------------------------------------
void MVCombination::PrintMatrix(TMatrixD &matrix, std::string name)
{
  int nrows = matrix.GetNrows();
  int ncols = matrix.GetNcols();

  std::cout << std::endl;
  std::cout << name.c_str() << " (" << nrows << "x" << ncols << "):" << std::endl;

  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) 
      std::cout << std::setprecision(3) << std::setw(7) << matrix[i][j] << " ";
    std::cout << std::endl;
  }
}

// ---------------------------------------------------------
void MVCombination::PrintVector(TVectorD &vector, std::string name)
{
  int nrows = vector.GetNoElements();

  std::cout << std::endl;
  std::cout << name.c_str() << " (" << nrows << "):" << std::endl;

  for (int i = 0; i < nrows; ++i) {
    std::cout << std::setprecision(3) << std::setw(7) << vector[i] << " ";
    std::cout << std::endl;
  }
}

// ---------------------------------------------------------
