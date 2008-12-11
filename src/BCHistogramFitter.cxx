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

	this -> MCMCSetNIterationsRun(2000);

	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

BCHistogramFitter::BCHistogramFitter(TH1D * hist, TF1 * func) : BCModel("HistogramFitter")
{
	this -> SetHistogram(hist);
	this -> SetFitFunction(func);

	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

void BCHistogramFitter::SetHistogram(TH1D * hist)
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
	double xmax = hist -> GetBinLowEdge(nbins+1);

	double histymin = hist -> GetMinimum();
	double histymax = hist -> GetMaximum();

	double ymin = TMath::Max(0., histymin - 5.*sqrt(histymin));
	double ymax = histymax + 5.*sqrt(histymax);

	this -> SetDataBoundaries(0, xmin, xmax);
	this -> SetDataBoundaries(1, ymin, ymax);

	this -> SetFitFunctionIndices(0, 1);
}

// ---------------------------------------------------------

void BCHistogramFitter::SetFitFunction(TF1 * func)
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

		this -> AddParameter(func->GetParName(i), xmin, xmax);
	}

	return;
}

// ---------------------------------------------------------

BCHistogramFitter::~BCHistogramFitter()
{}

// ---------------------------------------------------------

double BCHistogramFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	return 0;
}

// ---------------------------------------------------------

double BCHistogramFitter::LogLikelihood(std::vector <double> params)
{
	// initialize probability
	double loglikelihood = 0;

	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

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
	this -> FindModeMinuit( this -> GetBestFitParameters() );

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
	TH1D * h0;
	if(!opt.Contains("same"))
	{
		// create a copy of the histogram and draw it
		h0 = new TH1D(*fHistogram);
		h0 -> Draw(opt.Data());
	}

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
