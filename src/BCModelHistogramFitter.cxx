/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <iostream>
#include <fstream>

#include "BCModelHistogramFitter.h"

// ---------------------------------------------------------

BCModelHistogramFitter::BCModelHistogramFitter()
{

	fHistogram = 0; 
	fFitFunction = 0; 

}

// ---------------------------------------------------------

BCModelHistogramFitter::BCModelHistogramFitter(TH1D * hist, TF1 * func)
{
	this -> SetHistogram(hist); 
	this -> SetFitFunction(func);
}

// ---------------------------------------------------------

void BCModelHistogramFitter::SetHistogram(TH1D * hist)
{

	fHistogram = hist; 

	BCDataSet * ds = new BCDataSet(); 

	int nbins = fHistogram -> GetNbinsX(); 

	for (int i = 0; i < nbins; ++i)
		{
			BCDataPoint* dp = new BCDataPoint(2); 
			ds -> AddDataPoint(dp); 
		}

	this -> SetDataSet(ds); 

	double xmin = hist -> GetBinLowEdge(1); 
	double xmax = hist -> GetBinLowEdge(nbins) + hist -> GetBinWidth(1); 

	double histymin = hist -> GetMinimum(); 
	double histymax = hist -> GetMaximum(); 

	double ymin = TMath::Max(0., histymin - 5.*sqrt(histymin)); 
	double ymax = histymax + 5.*sqrt(histymax); 

	this -> SetDataBoundaries(0, xmin, xmax); 
	this -> SetDataBoundaries(1, ymin, ymax); 

	this -> SetFitFunctionIndices(0, 1);

}

// ---------------------------------------------------------

void BCModelHistogramFitter::SetFitFunction(TF1 * func)
{

	// set the function 
	fFitFunction = func; 

	// reset parameters 
	fParameterSet -> clear(); 

	// get the new number of parameters
	int n = func -> GetNpar(); 
	
	// add parameters
	for (int i = 0; i < n; ++i)
		{
			double xmin; 
			double xmax; 
			
			func -> GetParLimits(i, xmin, xmax); 

			this -> AddParameter(Form("parameter_%i", i), xmin, xmax);
		}
	
	return; 

}

// ---------------------------------------------------------

BCModelHistogramFitter::~BCModelHistogramFitter()
{

}

// ---------------------------------------------------------

double BCModelHistogramFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	
	return 0; 

}
  
// ---------------------------------------------------------

double BCModelHistogramFitter::LogLikelihood(std::vector <double> parameters)
{

	// initialize probability 
	double loglikelihood = 0; 

	// set the parameters of the function 
	for (int ipar = 0; ipar < int(parameters.size()); ++ipar)
		{
			fFitFunction -> SetParameter(ipar, parameters.at(ipar)); 
		}

	// get number of bins 
	int nbins = fHistogram -> GetNbinsX(); 

	// get bin width 
	double dx = fHistogram -> GetBinWidth(1); 

	// get function value of lower bin edge
	double fedgelow = fFitFunction -> Eval(fHistogram -> GetBinLowEdge(1)); 

	// loop over all bins 
	for (int ibin = 1; ibin <= nbins; ++ibin)
		{
			// get upper bin edge
			double xedgehi = fHistogram -> GetBinLowEdge(ibin) + dx; 

			// get function value at upper bin edge 
			double fedgehi = fFitFunction -> Eval(xedgehi); 

			// get the number of observed events 
			double y = fHistogram -> GetBinContent(ibin); 

			// get the number of expected events 
			// (integration via linear interpolation)
			double yexp = fedgelow * dx + (fedgehi - fedgelow) * dx / 2.; 

			// make the upper edge the lower edge for the next iteration
			fedgelow = fedgehi; 

			// get the value of the Poisson distribution 
			loglikelihood += BCMath::LogPoisson(y, yexp); 
		}

	// return the log likelihood 
	return loglikelihood; 

}

// --------------------------------------------------------- 

double BCModelHistogramFitter::FitFunction(std::vector <double> x, std::vector <double> parameters)
{

	// set the parameters of the function 
	for (int ipar = 0; ipar < int(parameters.size()); ++ipar)
		{
			fFitFunction -> SetParameter(ipar, parameters.at(ipar)); 
		}

	return fFitFunction -> Eval(x.at(0)) * fHistogram -> GetBinWidth(1); 

}

// // ---------------------------------------------------------

// void PrintHistogram(const char * options, const char * filename)
// {
	
// 	return; 

// }

// ---------------------------------------------------------
