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

#include "BCHistogramRatioFitter.h"

// ---------------------------------------------------------

BCHistogramRatioFitter::BCHistogramRatioFitter() : BCModel("HistogramFitter")
{
	fHistogram1 = 0;
	fHistogram2 = 0;
	fHistogramEfficiency = 0; 
	fFitFunction = 0;

	// set default options and values
	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
	fFlagIntegration = true; 
}

// ---------------------------------------------------------

BCHistogramRatioFitter::BCHistogramRatioFitter(TH1D * hist1, TH1D * hist2, TF1 * func) : BCModel("HistogramFitter")
{
	fHistogram1 = 0;
	fHistogram2 = 0;
	fHistogramEfficiency = 0; 
	fFitFunction = 0;

	this -> SetHistograms(hist1, hist2);
	this -> SetFitFunction(func);

	this -> MCMCSetNIterationsRun(2000);
	this -> SetFillErrorBand();
	fFlagIntegration = true; 
}

// ---------------------------------------------------------

int BCHistogramRatioFitter::SetHistograms(TH1D * hist1, TH1D * hist2)
{
	// check if histogram exists
	if (!hist1 || !hist2)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCHistogramRatioFitter::SetHistograms() : TH1D not created.");
		return 0;
	}

	// check compatibility of both histograms : number of bins 
	if (hist1 -> GetNbinsX() != hist2 -> GetNbinsX())
		{
		BCLog::Out(BCLog::error,BCLog::error,"BCHistogramRatioFitter::SetHistograms() : Histograms do not have the same number of bins.");
		return 0;
		}

	// check compatibility of both histograms : bin content 
	for (int i = 1; i <= hist1 -> GetNbinsX(); ++i)
		{
			if (hist1 -> GetBinContent(i) < hist2 -> GetBinContent(i))
				{
					BCLog::Out(BCLog::error,BCLog::error,"BCHistogramRatioFitter::SetHistograms() : Histogram 1 has fewer entries than histogram 2.");
					return 0;	
				}
		}
	
	// set pointer to histograms
	fHistogram1 = hist1;
	fHistogram2 = hist2;
	
	// create efficiency histogram 
	if (!fHistogramEfficiency)
		delete fHistogramEfficiency; 

	fHistogramEfficiency = new TH1D("hist_efficiency", 
																	"", 
																	hist1 -> GetNbinsX(), 
																	hist1 -> GetBinLowEdge(1),
																	hist1 -> GetBinLowEdge(hist1 -> GetNbinsX()+1)); 
	fHistogramEfficiency -> SetStats(kFALSE); 

	// create a data set. this is necessary in order to calculate the
	// error band. the data set contains as many data points as there
	// are bins. for now, the data points are empty.
	BCDataSet * ds = new BCDataSet();

	// create data points and add them to the data set.
	int nbins = fHistogram1 -> GetNbinsX();
	for (int i = 0; i < nbins; ++i)
	{
		BCDataPoint* dp = new BCDataPoint(2);
		ds -> AddDataPoint(dp);
	}

	// set the new data set.
	this -> SetDataSet(ds);

// 	// calculate the lower and upper edge in x.
// 	double xmin = hist1 -> GetBinLowEdge(1);
// 	double xmax = hist1 -> GetBinLowEdge(nbins+1);

// 	// calculate the minimum and maximum range in y.
// 	double histymin = hist2 -> GetMinimum();
// 	double histymax = hist1 -> GetMaximum();

// 	// calculate the minimum and maximum of the function value based on
// 	// the minimum and maximum value in y.
// 	double ymin = TMath::Max(0., histymin - 5.*sqrt(histymin));
// 	double ymax = histymax + 5.*sqrt(histymax);

// 	// set the data boundaries for x and y values.
// 	this -> SetDataBoundaries(0, xmin, xmax);
// 	this -> SetDataBoundaries(1, ymin, ymax);

	// set the indeces for fitting.
	this -> SetFitFunctionIndices(0, 1);

	// no error
	return 1;
}

// ---------------------------------------------------------

int BCHistogramRatioFitter::SetFitFunction(TF1 * func)
{
	// check if function exists
	if(!func)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCHistogramRatioFitter::SetFitFunction() : TF1 not created.");
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

BCHistogramRatioFitter::~BCHistogramRatioFitter()
{}

// ---------------------------------------------------------

double BCHistogramRatioFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	// using flat probability in all parameters
	return 0.;
}

// ---------------------------------------------------------

double BCHistogramRatioFitter::LogLikelihood(std::vector <double> params)
{
 	// initialize probability
 	double loglikelihood = 0;

	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	// get the number of bins
	int nbins = fHistogram1 -> GetNbinsX();

	// get bin width
	double dx = fHistogram1 -> GetBinWidth(1);

	// get function value of lower bin edge
	double fedgelow = fFitFunction -> Eval(fHistogram1 -> GetBinLowEdge(1));
	
	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
		{
			// get upper bin edge
			double xedgehi = fHistogram1 -> GetBinLowEdge(ibin) + dx;
			
			// get function value at upper bin edge
			double fedgehi = fFitFunction -> Eval(xedgehi);

			// get n 
			int n = int(fHistogram1 -> GetBinContent(ibin));

			// get k 
			int k = int(fHistogram2 -> GetBinContent(ibin));

			double eff = 0.; 

			// use ROOT's TH1D::Integral method
			if (fFlagIntegration)
				eff = fFitFunction -> Integral(xedgehi-dx, xedgehi);

			// use linear interpolation 
			else
				{
					eff = fedgelow * dx + (fedgehi - fedgelow) * dx / 2.; 
					
					// make the upper edge the lower edge for the next iteration
					fedgelow = fedgehi;
				}
			
			// get the value of the Poisson distribution
			loglikelihood += BCMath::LogApproxBinomial(n, k, eff);
		}
	
// 	// get bin boundaries
// 	double xmin = fHistogram -> GetBinLowEdge(ibin);
// 	double xmax = fHistogram -> GetBinLowEdge(ibin+1);
	
// 	// get the number of observed events
// 	double y = fHistogram -> GetBinContent(ibin);
	
// 	// get the number of expected events
// 	double yexp = fFitFunction -> Integral(xmin, xmax);
	
// 	// get the value of the Poisson distribution
// 	loglikelihood += BCMath::LogPoisson(y, yexp);

	return loglikelihood;
}

// ---------------------------------------------------------

double BCHistogramRatioFitter::FitFunction(std::vector <double> x, std::vector <double> params)
{
	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	return fFitFunction -> Eval(x[0]) * fHistogram1 -> GetBinWidth( fHistogram1 -> FindBin(x[0]) );
}

// ---------------------------------------------------------

int BCHistogramRatioFitter::Fit(TH1D * hist1, TH1D * hist2, TF1 * func)
{
	// set histogram
	if (hist1 && hist2)
		this -> SetHistograms(hist1, hist2);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::Fit() : Histogram(s) not defined.");
		return 0;
	}

	// set function
	if (func)
		this -> SetFitFunction(func);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::Fit() : Fit function not defined.");
		return 0;
	}

	// perform marginalization
	this -> MarginalizeAll();

	// maximize posterior probability, using the best-fit values close
	// to the global maximum from the MCMC
	this -> FindModeMinuit( this -> GetBestFitParameters() , -1);

	// calculate the p-value using the fast MCMC algorithm 
// 	double pvalue; 
// 	if (this -> CalculatePValueFast(this -> GetBestFitParameters(), pvalue))
// 		{
// 			fPValue = pvalue; 
// 		}
// 	else
// 		{
// 			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::Fit() : Could not use the fast p-value evaluation.");
// 		}

	// print summary to screen
	this -> PrintFitSummary();

	// no error
	return 1;
}

// ---------------------------------------------------------

void BCHistogramRatioFitter::DrawFit(const char * options, bool flaglegend)
{
// 	if (!fHistogram)
// 	{
// 		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::DrawFit() : Histogram not defined.");
// 		return;
// 	}

// 	if (!fFitFunction)
// 	{
// 		BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::DrawFit() : Fit function not defined.");
// 		return;
// 	}

// 	// check wheather options contain "same"
// 	TString opt = options;
// 	opt.ToLower();

// 	// if not same, draw the histogram first to get the axes
// 	if(!opt.Contains("same"))
// 		fHistogram -> Draw(opt.Data());

// 	// draw the error band as central 68% probability interval
// 	fErrorBand = this -> GetErrorBandGraph(0.16, 0.84);
// 	fErrorBand -> Draw("f same");

// 	// now draw the histogram again since it was covered by the band
// 	fHistogram -> Draw(TString::Format("%ssame",opt.Data()));

// 	// draw the fit function on top
// 	fGraphFitFunction = this -> GetFitFunctionGraph( this->GetBestFitParameters() );
// 	fGraphFitFunction -> SetLineColor(kRed);
// 	fGraphFitFunction -> SetLineWidth(2);
// 	fGraphFitFunction -> Draw("l same");

// 	// draw legend
// 	if (flaglegend)
// 		{
// 			TLegend * legend = new TLegend(0.25, 0.75, 0.55, 0.95); 
// 			legend -> SetBorderSize(0); 
// 			legend -> SetFillColor(kWhite); 
// 			legend -> AddEntry(fHistogram, "Data", "L"); 
// 			legend -> AddEntry(fGraphFitFunction, "Best fit", "L"); 
// 			legend -> AddEntry(fErrorBand, "Error band", "F"); 
// 			legend -> Draw(); 
// 		}

// 	gPad -> RedrawAxis();
}

// ---------------------------------------------------------

void BCHistogramRatioFitter::PrintFitSummary()
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
//int BCHistogramRatioFitter::CalculatePValueFast(std::vector<double> par, double &pvalue)
//{
// 	// check size of parameter vector 
// 	if (par.size() != this -> GetNParameters())
// 		{
// 			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::CalculatePValueFast() : Number of parameters is inconsistent.");
// 		return 0;
// 		}

// 	// check if histogram exists
// 	if (!fHistogram)
// 		{
// 			BCLog::Out(BCLog::warning, BCLog::warning,"BCHistogramRatioFitter::CalculatePValueFast() : Histogram not defined.");
// 		return 0;
// 		}

// 	// define temporary variables
// 	int nbins = fHistogram -> GetNbinsX(); 

// 	std::vector <int> histogram; 
// 	std::vector <double> expectation; 
// 	histogram.assign(nbins, 0); 
// 	expectation.assign(nbins, 0); 

// 	double logp = 0; 
// 	double logp_start = 0; 
// 	int counter_pvalue = 0; 

// 	// define starting distribution
// 	for (int ibin = 0; ibin < nbins; ++ibin)
// 		{
// 			// get bin boundaries
// 			double xmin = fHistogram -> GetBinLowEdge(ibin+1);
// 			double xmax = fHistogram -> GetBinLowEdge(ibin+2);
			
// 			// get the number of expected events
// 			double yexp = fFitFunction -> Integral(xmin, xmax);

// 			histogram[ibin]   = int(yexp); 
// 			expectation[ibin] = yexp; 

// 			// calculate p; 
// 			logp += BCMath::LogPoisson(double(int(yexp)), yexp); 
// 			logp_start += BCMath::LogPoisson(fHistogram -> GetBinContent(ibin+1), yexp); 
// 		}

// 	int niter = 100000; 
	
// 	// loop over iterations 
// 	for (int iiter = 0; iiter < niter; ++iiter)
// 		{
			
// 			// loop over bins 
// 			for (int ibin = 0; ibin < nbins; ++ibin)
// 				{
// 					// random step up or down in statistics for this bin 
// 					double ptest = fRandom -> Rndm() - 0.5; 

// 					// increase statistics by 1
// 					if (ptest > 0)
// 						{
// 							// calculate factor of probability 
// 							double r = expectation.at(ibin) / double(histogram.at(ibin) + 1); 

// 							// walk, or don't (this is the Metropolis part) 
// 							if (fRandom -> Rndm() < r)
// 								{
// 									histogram[ibin] = histogram.at(ibin) + 1; 
// 									logp += log(r); 
// 								}
// 						}

// 					// decrease statistics by 1 
// 					else
// 						{
// 							// calculate factor of probability 
// 							double r = double(histogram.at(ibin)) / expectation.at(ibin); 

// 							// walk, or don't (this is the Metropolis part) 
// 							if (fRandom -> Rndm() < r)
// 								{
// 									histogram[ibin] = histogram.at(ibin) - 1; 
// 									logp += log(r); 
// 								}
// 						} 					
// 				} // end of looping over bins 

// 			// increase counter 
// 			if (logp < logp_start)
// 				counter_pvalue++; 

// 		} // end of looping over iterations 

// 	// calculate p-value 
// 	pvalue = double(counter_pvalue) / double(niter); 

// 	// no error 
// 	return 1; 
//}

// ---------------------------------------------------------
