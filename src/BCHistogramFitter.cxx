/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <iostream>
#include <fstream>

#include "BCHistogramFitter.h"

// ---------------------------------------------------------

BCHistogramFitter::BCHistogramFitter() : BCModel("HistogramFitter")
{
	fHistogram = 0;
	fFitFunction = 0;

	// set default options and values
	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

BCHistogramFitter::BCHistogramFitter(TH1D * hist, TF1 * func) : BCModel("HistogramFitter")
{
	this -> SetHistogram(hist);
	this -> SetFitFunction(func);

	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

int BCHistogramFitter::SetHistogram(TH1D * hist)
{
	// check if histogram exists
	if(!hist)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCHistogramFitter::SetHistogram() : TH1D not created.");
		return 0;
	}

	// set pointer to histogram
	fHistogram = hist;

	// create a data set. this is necessary in order to calculate the
	// error band. the data set contains as many data points as there
	// are bins. for now, the data points are empty. 
	BCDataSet * ds = new BCDataSet();

	// create data points and add them to the data set. 
	int nbins = fHistogram -> GetNbinsX();
	for (int i = 0; i < nbins; ++i)
	{
		BCDataPoint* dp = new BCDataPoint(2);
		ds -> AddDataPoint(dp);
	}

	// set the new data set. 
	this -> SetDataSet(ds);

	// calculate the lower and upper edge in x. 
	double xmin = hist -> GetBinLowEdge(1);
	double xmax = hist -> GetBinLowEdge(nbins+1);

	// calculate the minimum and maximum range in y. 
	double histymin = hist -> GetMinimum();
	double histymax = hist -> GetMaximum();

	// calculate the minimum and maximum of the function value based on
	// the minimum and maximum value in y. 
	double ymin = TMath::Max(0., histymin - 5.*sqrt(histymin));
	double ymax = histymax + 5.*sqrt(histymax);

	// set the data boundaries for x and y values. 
	this -> SetDataBoundaries(0, xmin, xmax);
	this -> SetDataBoundaries(1, ymin, ymax);

	// set the indeces for fitting.
	this -> SetFitFunctionIndices(0, 1);

	// no error 
	return 1;
}

// ---------------------------------------------------------

int BCHistogramFitter::SetFitFunction(TF1 * func)
{
	// check if function exists
	if(!func)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCHistogramFitter::SetFitFunction() : TF1 not created.");
		return 0;
	}

	// set the function
	fFitFunction = func;

	// update the model name to contain the function name
	this -> SetName(TString::Format("HistogramFitter with %s",fFitFunction->GetName()));

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

		this -> AddParameter(func->GetParName(i), xmin, xmax);
	}

	// no error 
	return 1;
}

// ---------------------------------------------------------

BCHistogramFitter::~BCHistogramFitter()
{}

// ---------------------------------------------------------

double BCHistogramFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	// using flat probability in all parameters
	return 0.;
}

// ---------------------------------------------------------

double BCHistogramFitter::LogLikelihood(std::vector <double> params)
{
	// initialize probability
	double loglikelihood = 0;

	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	// get the number of bins
	int nbins = fHistogram -> GetNbinsX();

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		// get bin boundaries
		double xmin = fHistogram -> GetBinLowEdge(ibin);
		double xmax = fHistogram -> GetBinLowEdge(ibin+1);

		// get the number of observed events
		double y = fHistogram -> GetBinContent(ibin);

		// get the number of expected events
		double yexp = fFitFunction -> Integral(xmin, xmax);

		// get the value of the Poisson distribution
		loglikelihood += BCMath::LogPoisson(y, yexp);
	}

	return loglikelihood;
}

// ---------------------------------------------------------

double BCHistogramFitter::FitFunction(std::vector <double> x, std::vector <double> params)
{
	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	return fFitFunction -> Eval(x[0]) * fHistogram -> GetBinWidth( fHistogram -> FindBin(x[0]) );
}

// ---------------------------------------------------------

int BCHistogramFitter::Fit(TH1D * hist, TF1 * func)
{
	// set histogram
	if (hist)
		this -> SetHistogram(hist);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::Fit() : Histogram not defined.");
		return 0;
	}

	// set function
	if (func)
		this -> SetFitFunction(func);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::Fit() : Fit function not defined.");
		return 0;
	}

	// perform marginalization
	this -> MarginalizeAll();

	// maximize posterior probability, using the best-fit values close
	// to the global maximum from the MCMC
	this -> FindModeMinuit( this -> GetBestFitParameters() , -1);

	// print summary to screen
	this -> PrintFitSummary(); 

	// no error 
	return 1;
}

// ---------------------------------------------------------

void BCHistogramFitter::DrawFit(const char * options)
{
	if (!fHistogram)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::DrawFit() : Histogram not defined.");
		return;
	}

	if (!fFitFunction)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::DrawFit() : Fit function not defined.");
		return;
	}

	// check wheather options contain "same"
	TString opt = options;
	opt.ToLower();

	// if not same, draw the histogram first to get the axes
	if(!opt.Contains("same"))
		fHistogram -> Draw(opt.Data());

	// draw the error band as central 68% probability interval
	this -> GetErrorBandGraph(0.16, 0.84) -> Draw("f same");

	// now draw the histogram again since it was covered by the band
	fHistogram -> Draw(TString::Format("%ssame",opt.Data()));

	// draw the fit function on top
	TGraph * gfit = this -> GetFitFunctionGraph( this->GetBestFitParameters() );
	gfit -> SetLineColor(kRed);
	gfit -> SetLineWidth(2);
	gfit -> Draw("l same");

	gPad -> RedrawAxis();
}

// ---------------------------------------------------------

void BCHistogramFitter::PrintFitSummary()
{
	std::cout << std::endl;
	std::cout << "Fit summary " << std::endl; 
	std::cout << "------------------------------------ " << std::endl;

	std::cout << "Number of parameters : " 
						<< this -> GetNParameters() << std::endl; 
	std::cout << std::endl; 

	std::cout << "Best fit parameters (global) : " << std::endl; 
	for (int i = 0; i < this -> GetNParameters(); ++i)
		std::cout << this -> GetParameter(i) -> GetName() << " : " 
							<< this -> GetBestFitParameter(i) << std::endl; 
	std::cout << std::endl; 

// 	std::cout << "Best fit parameters (marginalized) : " << std::endl; 
// 	for (int i = 0; i < this -> GetNParameters(); ++i)
// 		{
// 			BCH1D * bch1d = this -> GetMarginalized(fParameterSet -> at(i));
// 			std::cout << this -> GetParameter(i) -> GetName() << " : " 
// 								<< this -> GetBestFitParameterMarginalized(i) << std::endl; 
// 		}
// 	std::cout << std::endl; 
	
}

// ---------------------------------------------------------
