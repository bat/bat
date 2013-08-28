#include "RatioModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

#include <cmath>

// ---------------------------------------------------------
RatioModel::RatioModel()
 : BCModel()
 , fHistRatio(new BCH1D())
{
  // define the parameters x and y
  DefineParameters();

  // define the additional histogram for the ratio r
  DefineHistogram();
};

// ---------------------------------------------------------
RatioModel::RatioModel(const char * name)
 : BCModel(name)
 , fHistRatio(new BCH1D())
{
  // define the parameters x and y
  DefineParameters();

  // define the additional histogram for the ratio r
  DefineHistogram();
};

// ---------------------------------------------------------
RatioModel::~RatioModel()
{
  // free memory
  delete fHistRatio;
};

// ---------------------------------------------------------
void RatioModel::DefineParameters()
{
  // add two parameters, x and y.
  AddParameter("x", 0., 8.); // index 0
  AddParameter("y", 0., 16.); // index 1
}

// ---------------------------------------------------------
void RatioModel::DefineHistogram()
{
  // create a new ROOT and BAT histogram
  TH1D* hist = new TH1D("", ";r;p(r|x,y)", 100, 0.0, 2.0);
  hist->SetStats(kFALSE);
  fHistRatio->SetHistogram(hist);
}

// ---------------------------------------------------------
void RatioModel::PrintHistogram()
{
  // print the BAT histogram to an eps file
  fHistRatio->Print("ratio.pdf");
}

// ---------------------------------------------------------
double RatioModel::LogLikelihood(const std::vector<double> &parameters)
{
  // This methods returns the logarithm of the conditional probability
  // p(data|parameters). This is where you have to define your model.

  double logprob = 0.;

  double x = parameters.at(0);
  double y = parameters.at(1);

  // calculate the Gaussian probability densities
  logprob += BCMath::LogGaus(x, 4., 1.);
  logprob += BCMath::LogGaus(y, 8., 2.0);

  // return the log likelihood
  return logprob;
}

// ---------------------------------------------------------
double RatioModel::LogAPrioriProbability(const std::vector<double> &parameters)
{
  // This method returns the logarithm of the prior probability for the
  // parameters p(parameters).

  double logprob = 0.;

  // get width of the parameter ranges
  double dx = GetParameter(0)->GetRangeWidth();
  double dy = GetParameter(1)->GetRangeWidth();

  // add flat prior probabilities for x and y
  logprob += log(1./dx); // flat prior for x
  logprob += log(1./dy); // flat prior for y

  // return log prior
  return logprob;
}

// ---------------------------------------------------------
void RatioModel::MCMCUserIterationInterface()
{
  // get number of chains
  int nchains = MCMCGetNChains();

  // get number of parameters
  int npar = GetNParameters();

  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    // get the current values of the parameters x and y. These are
    // stored in fMCMCx.
    double x = fMCMCx.at(i * npar + 0);
    double y = fMCMCx.at(i * npar + 1);

    // fill the ratio histogram
    fHistRatio->GetHistogram()->Fill(x/y);
  }
}

// ---------------------------------------------------------
