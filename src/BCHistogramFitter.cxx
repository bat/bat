/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TLegend.h> 

#include "BCLog.h"
#include "BCDataSet.h"
#include "BCDataPoint.h"
#include "BCMath.h"

#include "BCHistogramFitter.h"

// ---------------------------------------------------------

BCHistogramFitter::BCHistogramFitter() : BCModel("HistogramFitter")
{
	fHistogram = 0;
	fFitFunction = 0;

	// set default options and values
	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
	fFlagIntegration = true; 
}

// ---------------------------------------------------------

BCHistogramFitter::BCHistogramFitter(TH1D * hist, TF1 * func) : BCModel("HistogramFitter")
{
	fHistogram = 0;
	fFitFunction = 0;

	this -> SetHistogram(hist);
	this -> SetFitFunction(func);

	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
	fFlagIntegration = true; 
}

// ---------------------------------------------------------

int BCHistogramFitter::SetHistogram(TH1D * hist)
{
	// check if histogram exists
	if (!hist)
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

			double yexp = 0.; 

			// use ROOT's TH1D::Integral method
			if (fFlagIntegration)
				yexp = fFitFunction -> Integral(xedgehi-dx, xedgehi);

			// use linear interpolation 
			else
				{
					yexp = fedgelow * dx + (fedgehi - fedgelow) * dx / 2.; 

					// make the upper edge the lower edge for the next iteration
					fedgelow = fedgehi;
				}

			// get the value of the Poisson distribution
			loglikelihood += BCMath::LogPoisson(y, yexp);
		}

// 			// get bin boundaries
// 			double xmin = fHistogram -> GetBinLowEdge(ibin);
// 			double xmax = fHistogram -> GetBinLowEdge(ibin+1);
			
// 			// get the number of observed events
// 			double y = fHistogram -> GetBinContent(ibin);
			
// 			// get the number of expected events
// 			double yexp = fFitFunction -> Integral(xmin, xmax);
			
// 			// get the value of the Poisson distribution
// 			loglikelihood += BCMath::LogPoisson(y, yexp);

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

	// calculate the p-value using the fast MCMC algorithm 
	double pvalue; 
	if (this -> CalculatePValueFast(this -> GetBestFitParameters(), pvalue))
		{
			fPValue = pvalue; 
		}
	else
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::Fit() : Could not use the fast p-value evaluation.");
		}

	// print summary to screen
	this -> PrintFitSummary();

	// no error
	return 1;
}

// ---------------------------------------------------------

void BCHistogramFitter::DrawFit(const char * options, bool flaglegend)
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
	fErrorBand = this -> GetErrorBandGraph(0.16, 0.84);
	fErrorBand -> Draw("f same");

	// now draw the histogram again since it was covered by the band
	fHistogram -> Draw(TString::Format("%ssame",opt.Data()));

	// draw the fit function on top
	fGraphFitFunction = this -> GetFitFunctionGraph( this->GetBestFitParameters() );
	fGraphFitFunction -> SetLineColor(kRed);
	fGraphFitFunction -> SetLineWidth(2);
	fGraphFitFunction -> Draw("l same");

	// draw legend
	if (flaglegend)
		{
			TLegend * legend = new TLegend(0.25, 0.75, 0.55, 0.95); 
			legend -> SetBorderSize(0); 
			legend -> SetFillColor(kWhite); 
			legend -> AddEntry(fHistogram, "Data", "L"); 
			legend -> AddEntry(fGraphFitFunction, "Best fit", "L"); 
			legend -> AddEntry(fErrorBand, "Error band", "F"); 
			legend -> Draw(); 
		}

	gPad -> RedrawAxis();
}

// ---------------------------------------------------------

void BCHistogramFitter::PrintFitSummary()
{
	BCLog::Out(BCLog::summary, BCLog::summary, "-----------------------------------------"); 
	BCLog::Out(BCLog::summary, BCLog::summary, "Fit summary:");
	BCLog::Out(BCLog::summary, BCLog::summary, Form("Number of parameters = %i", this -> GetNParameters())); 

	BCLog::Out(BCLog::summary, BCLog::summary, "Best fit parameters (global):"); 
	for (unsigned int i = 0; i < this -> GetNParameters(); ++i)
		BCLog::Out(BCLog::summary, BCLog::summary, Form("%s = %.2lf", this -> GetParameter(i) -> GetName().data(), this -> GetBestFitParameter(i)));
	
	BCLog::Out(BCLog::summary, BCLog::summary, "Goodness-of-fit test:");
	BCLog::Out(BCLog::summary, BCLog::summary, Form("p-value = %.2lf", this -> GetPValue())); 
	BCLog::Out(BCLog::summary, BCLog::summary, "-----------------------------------------"); 
}

// ---------------------------------------------------------
int BCHistogramFitter::CalculatePValueFast(std::vector<double> par, double &pvalue)
{
	// check size of parameter vector 
	if (par.size() != this -> GetNParameters())
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::CalculatePValueFast() : Number of parameters is inconsistent.");
		return 0;
		}

	// check if histogram exists
	if (!fHistogram)
		{
			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramFitter::CalculatePValueFast() : Histogram not defined.");
		return 0;
		}

	// define temporary variables
	int nbins = fHistogram -> GetNbinsX(); 

	std::vector <int> histogram; 
	std::vector <double> expectation; 
	histogram.assign(nbins, 0); 
	expectation.assign(nbins, 0); 

	double logp = 0; 
	double logp_start = 0; 
	int counter_pvalue = 0; 

	// define starting distribution
	for (int ibin = 0; ibin < nbins; ++ibin)
		{
			// get bin boundaries
			double xmin = fHistogram -> GetBinLowEdge(ibin+1);
			double xmax = fHistogram -> GetBinLowEdge(ibin+2);
			
			// get the number of expected events
			double yexp = fFitFunction -> Integral(xmin, xmax);

			histogram[ibin]   = int(yexp); 
			expectation[ibin] = yexp; 

			// calculate p; 
			logp += BCMath::LogPoisson(double(int(yexp)), yexp); 
			logp_start += BCMath::LogPoisson(fHistogram -> GetBinContent(ibin+1), yexp); 
		}

	int niter = 100000; 
	
	// loop over iterations 
	for (int iiter = 0; iiter < niter; ++iiter)
		{
			
			// loop over bins 
			for (int ibin = 0; ibin < nbins; ++ibin)
				{
					// random step up or down in statistics for this bin 
					double ptest = fRandom -> Rndm() - 0.5; 

					// increase statistics by 1
					if (ptest > 0)
						{
							// calculate factor of probability 
							double r = expectation.at(ibin) / double(histogram.at(ibin) + 1); 

							// walk, or don't (this is the Metropolis part) 
							if (fRandom -> Rndm() < r)
								{
									histogram[ibin] = histogram.at(ibin) + 1; 
									logp += log(r); 
								}
						}

					// decrease statistics by 1 
					else
						{
							// calculate factor of probability 
							double r = double(histogram.at(ibin)) / expectation.at(ibin); 

							// walk, or don't (this is the Metropolis part) 
							if (fRandom -> Rndm() < r)
								{
									histogram[ibin] = histogram.at(ibin) - 1; 
									logp += log(r); 
								}
						} 					
				} // end of looping over bins 

			// increase counter 
			if (logp < logp_start)
				counter_pvalue++; 

		} // end of looping over iterations 

	// calculate p-value 
	pvalue = double(counter_pvalue) / double(niter); 

	// no error 
	return 1; 
}

// ---------------------------------------------------------
