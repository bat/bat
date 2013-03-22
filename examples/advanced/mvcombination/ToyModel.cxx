#include "ToyModel.h"

#include <iomanip>

#include <TMath.h>
#include <TCanvas.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>

// ---------------------------------------------------------
ToyModel::ToyModel() : BCModel("ToyModel")
{
}

// ---------------------------------------------------------
ToyModel::~ToyModel()
{
}

// ---------------------------------------------------------
void ToyModel::SetNMeasurements(int n, double min, double max)
{ 
  for (int i = 0; i < n; ++i)
    AddParameter(Form("measurement_%i", i), min, max);
}

// ---------------------------------------------------------
void ToyModel::SetCovarianceMatrix(TMatrixD matrix) 
{ 
  fCovarianceMatrix.Clear();
  fCovarianceMatrix.ResizeTo(matrix);
  fCovarianceMatrix = matrix; 

  fInvCovarianceMatrix.Clear();
  fInvCovarianceMatrix.ResizeTo(fCovarianceMatrix);
  fInvCovarianceMatrix = fCovarianceMatrix;
  fInvCovarianceMatrix.Invert();

  fDetCovariance = fCovarianceMatrix.Determinant();
};

// ---------------------------------------------------------
void ToyModel::SetParameters(std::vector<double> parameters)
{ 
  fPars = parameters; 
  int nmeas = GetNParameters();

  for (int i = 0; i < nmeas; ++i) {
    fVectorObservables[i] = fPars[fVectorObservable[i]];
  }
};

// ---------------------------------------------------------
double ToyModel::LogLikelihood(const std::vector<double> &parameters)
{
  double logprob = 0.;

  int nmeas = GetNParameters();

  // copy parameters into a vector
  TVectorD observables(nmeas);
  TVectorD measurements(nmeas);

  for (int i = 0; i < nmeas; ++i) {
    observables[i] = fPars[fVectorObservable[i]];
    measurements[i] = parameters[i];
  }

  TVectorD prod1 = observables - measurements;
  TVectorD prod2 = fInvCovarianceMatrix * prod1;
  double prod = prod1 * prod2;

  logprob = -0.5 * prod - log(TMath::Power(2*TMath::Pi(), nmeas/2.) * sqrt(fDetCovariance));

  return logprob;
}

// ---------------------------------------------------------
double ToyModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  double logprob = 0.;

  return logprob;
}

// ---------------------------------------------------------
void ToyModel::MCMCIterationInterface()
{
  // get number of chains
  int nchains = MCMCGetNChains();

  // get number of parameters
  int npar = GetNParameters();

  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    // get the current values of the parameters x and y. These are
		
    // copy parameters into a vector
    TVectorD observables(npar);
    TVectorD measurements(npar);
		
    for (int j = 0; j < npar; ++j) {
      observables[j] = fPars[fVectorObservable[j]];
      measurements[j] = fMCMCx.at(i * npar + j);
    }

    double chi2 = Chi2(observables, measurements);

    // fill the histogram
    fHistChi2->Fill(chi2);
  }
}

// ---------------------------------------------------------
void ToyModel::PrintToys(std::string filename)
{
  // check if file extension is pdf
  if ( (filename.find_last_of(".") != std::string::npos) &&
       (filename.substr(filename.find_last_of(".")+1) == "pdf") ) {
    ; // it's a PDF file
    
  }
  else if ( (filename.find_last_of(".") != std::string::npos) &&
	    (filename.substr(filename.find_last_of(".")+1) == "ps") ) {
    ; // it's a PS file
  }
  else {
    ; // make it a PDF file
    filename += ".pdf";
  }
  
  TCanvas* c1 = new TCanvas("");
  c1->cd();

  // draw exoected chi2 distribution
  BCH1D* hist_chi2 = new BCH1D(fHistChi2);
  hist_chi2->Draw();

  // calculate observed chi2
  // debugKK: I am here
  //  std::cout << " chi2 : " << Chi2(fVectorObservables, fVectorMeasurements) << std::endl;

  c1->Print(filename.c_str());
}

// ---------------------------------------------------------
double ToyModel::Chi2(TVectorD observables, TVectorD measurements)
{
  TVectorD prod1 = observables - measurements;
  TVectorD prod2 = fInvCovarianceMatrix * prod1;
  double chi2 = prod1 * prod2;

  return chi2;
}

// ---------------------------------------------------------
