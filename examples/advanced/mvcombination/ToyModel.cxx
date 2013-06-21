#include "ToyModel.h"

#include <iomanip>

#include <TMath.h>
#include <TCanvas.h>
#include <TMarker.h>

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>

// ---------------------------------------------------------
ToyModel::ToyModel(MVCombination* mvc) : BCModel("ToyModel")
{
  SetNMeasurements(mvc->GetNMeasurements(),
		   mvc->GetParameter(0)->GetLowerLimit(),
		   mvc->GetParameter(0)->GetUpperLimit());
  SetVectorMeasurements(mvc->GetVectorMeasurements());
  SetVectorObservable(mvc->GetVectorObservable());
  SetCovarianceMatrix(mvc->GetCovarianceMatrix());
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

  fVectorObservables.Clear();
  fVectorObservables.ResizeTo(nmeas);

  for (int i = 0; i < nmeas; ++i) {
    fVectorObservables[i] = fPars[fVectorObservable[i]];
  }
};

// ---------------------------------------------------------
void ToyModel::SetMeasurementRanges(std::vector<double> min, std::vector<double> max)
{
  int npar = GetNParameters();

  if ( (int(min.size()) != npar) || (int(max.size()) != npar) ) {
    BCLog::OutWarning("ToyModel::SetMeasurementRanges. Size of ranges does not fit the number of measurements.");
    return;
  }

  // set the parameter ranges
  for (int i = 0; i < npar; ++i) {
      GetParameter(i)->SetLimits(min.at(i), max.at(i));
  }

}

// ---------------------------------------------------------
void ToyModel::SetMeasurementRanges(double min, double max)
{
  int npar = GetNParameters();

  std::vector<double> min_vec;
  std::vector<double> max_vec;

  // fill vector with ranges
  for (int i = 0; i < npar; ++i) {
    min_vec.push_back(min);
    max_vec.push_back(max);
  }

  SetMeasurementRanges(min_vec, max_vec);
}

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
    if (fHistChi2)
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

  int npars = GetNParameters();

  TCanvas* c1 = new TCanvas("");
  c1->cd();

  // calculate observed chi2
  double chi2 =  Chi2(fVectorObservables, fVectorMeasurements);

  // calculate expected chi2 distribution
  BCH1D* hist_chi2 = new BCH1D(fHistChi2);
  hist_chi2->GetHistogram()->Scale(1.0/hist_chi2->GetHistogram()->Integral("width"));

  // calculate p-value
  double pvalue = hist_chi2->GetHistogram()->Integral(hist_chi2->GetHistogram()->FindBin(chi2), hist_chi2->GetHistogram()->GetNbinsX(), "width");

  // draw expected chi2 distribution
  hist_chi2->Draw("BTllB1CS1L", 1-pvalue);

  c1->Print(std::string(filename+"(").c_str());

  // draw all 1D histograms indicating observed value
  for (int i = 0; i < npars; ++i) {
    double obs = fVectorMeasurements[i];
    const BCParameter* par = GetParameter(i);
    BCH1D* hist_par = GetMarginalized(par);
    double p = hist_par->GetHistogram()->Integral(hist_par->GetHistogram()->FindBin(obs), hist_par->GetHistogram()->GetNbinsX(), "width");
    hist_par->Draw("BTllB1CS1L", 1-p);
    c1->Print(filename.c_str());
  }

  // draw all 2D histograms indicating observed values
  for (int i = 0; i < npars; ++i) {
    for (int j = 0; j < i; ++j) {
      double obs_i = fVectorMeasurements[i];
      double obs_j = fVectorMeasurements[j];
      const BCParameter* par_i = GetParameter(i);
      const BCParameter* par_j = GetParameter(j);
      BCH2D* hist_par = GetMarginalized(par_i, par_j);
      hist_par->Draw("BTfB3CS1");
      TMarker* m = new TMarker();
      m->SetMarkerStyle(21);
      m->DrawMarker(obs_j, obs_i);
      if (i == npars - 1 && j == i-1)
	c1->Print(std::string(filename+")").c_str());
      else
	c1->Print(filename.c_str());
    }
  }

}

// ---------------------------------------------------------
void ToyModel::PrintSummary() {

  std::cout << " Goodness-of-fit test:" << std::endl << std::endl;

  // calculate observed chi2
  double chi2 =  Chi2(fVectorObservables, fVectorMeasurements);

  // calculate expected chi2 distribution
  BCH1D* hist_chi2 = new BCH1D(fHistChi2);
  hist_chi2->GetHistogram()->Scale(1.0/hist_chi2->GetHistogram()->Integral("width"));

  // calculate p-value
  double pvalue = hist_chi2->GetHistogram()->Integral(hist_chi2->GetHistogram()->FindBin(chi2), hist_chi2->GetHistogram()->GetNbinsX(), "width");

  std::cout << " chi2    : " << chi2 << std::endl;
  std::cout << " p-value : " << pvalue << std::endl;
  std::cout << std::endl;

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
