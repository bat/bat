// ***************************************************************
// This file was created using the ./CreateFitModel.sh script
// ./CreateFitModel.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <TemplateFit.h>
#include <TMath.h>
// ---------------------------------------------------------
TemplateFit::TemplateFit() : BinLikeModel()
{  
  DefineParameters();
};

// ---------------------------------------------------------
TemplateFit::~TemplateFit()
{}; 

// ---------------------------------------------------------
void TemplateFit::DefineParameters()
{
  // add model parameters 
	AddParameter("mean_offset",   100.0,  200.0); 
	AddParameter("mean_slope",      0.0,    0.4); 
	AddParameter("width_offset", - 20.0,   20.0); 
	AddParameter("width_slope",  -  0.1,    0.25); 
}

// ---------------------------------------------------------
double TemplateFit::LogAPrioriProbability(std::vector <double> parameters)
{
	double logprob = 0.;

	return logprob;
}

// ---------------------------------------------------------
double TemplateFit::Expectation(std::vector <double> parameters, double parvalue, double x)
{
  // initialize expectation
  double expectation = 0; 

 	// calculate parameters of fit function from model parameters and
	// leading parameter
	double mean = parameters.at(0) + parameters.at(1) * parvalue; 
	double rms = parameters.at(2) + parameters.at(3) * parvalue; 
 
	if (rms<=0)
		return 0; 

	expectation = TMath::Gaus(x, mean, rms, true); 

 // return expectation 
  return expectation;
}
// ---------------------------------------------------------

