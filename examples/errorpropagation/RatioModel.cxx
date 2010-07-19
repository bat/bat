#include "RatioModel.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
RatioModel::RatioModel() : BCModel()
			 , fHistRatio(new BCH1D())
{  
  DefineParameters();

  DefineHistogram(); 
};

// ---------------------------------------------------------
RatioModel::RatioModel(const char * name) : BCModel(name)
{ 
  DefineParameters();
};

// ---------------------------------------------------------
RatioModel::~RatioModel()
{
  delete fHistRatio;
};

// ---------------------------------------------------------
void RatioModel::DefineParameters()
{
  AddParameter("x", 1., 17.0); // index 0
  AddParameter("y", 6., 14.); // index 1
}

// ---------------------------------------------------------
void RatioModel::DefineHistogram()
{
  TH1D* hist = new TH1D("", ";r;p(r|x,y)", 100, 0.0, 2.0);
  fHistRatio->SetHistogram(hist);
}

// ---------------------------------------------------------
void RatioModel::PrintHistogram()
{
  fHistRatio->Print("ratio.eps");
}

// ---------------------------------------------------------
double RatioModel::LogLikelihood(std::vector <double> parameters)
{
  double logprob = 0.;

  double x = parameters.at(0);
  double y = parameters.at(1);

  logprob += BCMath::LogGaus(x, 9., 2.); 
  logprob += BCMath::LogGaus(y, 10., 1.0); 

  return logprob;
}

// ---------------------------------------------------------
double RatioModel::LogAPrioriProbability(std::vector <double> parameters)
{
  double logprob = 0.;

  double dx = GetParameter(0)->GetRangeWidth(); 
  double dy = GetParameter(1)->GetRangeWidth(); 

  logprob += log(1./dx);                    // flat prior for x
  logprob += log(1./dy);                    // flat prior for y

  return logprob;
}

// ---------------------------------------------------------
void RatioModel::MCMCIterationInterface()
{
  // get number of chains
  int nchains = MCMCGetNChains();

  // get number of parameters
  int npar = GetNParameters();

  // loop over all chains and fill histogram
  for (int i = 0; i < nchains; ++i) {
    double x = fMCMCx.at(i * npar + 0); 
    double y = fMCMCx.at(i * npar + 1); 
    
    fHistRatio->GetHistogram()->Fill(x/y);
  }
}

// ---------------------------------------------------------
