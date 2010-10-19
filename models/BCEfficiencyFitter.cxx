/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TMath.h>

#include "BAT/BCLog.h"
#include "BAT/BCDataSet.h"
#include "BAT/BCDataPoint.h"
#include "BAT/BCMath.h"

#include "BCEfficiencyFitter.h"

// ---------------------------------------------------------

BCEfficiencyFitter::BCEfficiencyFitter() : BCModel("EfficiencyFitter")
{
	fHistogram1 = 0;
	fHistogram2 = 0;
	fHistogramBinomial = 0;
	fHistogramRatio = 0;
	fFitFunction = 0;

	// set default options and values
	this -> MCMCSetNIterationsRun(2000);
	this -> MCMCSetRValueCriterion(0.01);
	this -> SetFillErrorBand();
	fFlagIntegration = false;
}

// ---------------------------------------------------------

BCEfficiencyFitter::BCEfficiencyFitter(TH1D * hist1, TH1D * hist2, TF1 * func) : BCModel("EfficiencyFitter")
{
	fHistogram1 = 0;
	fHistogram2 = 0;
	fHistogramBinomial = 0;
	fHistogramRatio = 0;
	fFitFunction = 0;

	this -> SetHistograms(hist1, hist2);
	this -> SetFitFunction(func);

	this -> MCMCSetNIterationsRun(2000);
	this -> MCMCSetRValueCriterion(0.01);
	this -> SetFillErrorBand();
	fFlagIntegration = false;
}

// ---------------------------------------------------------

int BCEfficiencyFitter::SetHistograms(TH1D * hist1, TH1D * hist2)
{
	// check if histogram exists
	if (!hist1 || !hist2)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCEfficiencyFitter::SetHistograms() : TH1D not created.");
		return 0;
	}

	// check compatibility of both histograms : number of bins
	if (hist1 -> GetNbinsX() != hist2 -> GetNbinsX())
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCEfficiencyFitter::SetHistograms() : Histograms do not have the same number of bins.");
		return 0;
	}

	// check compatibility of both histograms : bin content
	for (int i = 1; i <= hist1 -> GetNbinsX(); ++i)
	{
		if (hist1 -> GetBinContent(i) < hist2 -> GetBinContent(i))
		{
			BCLog::Out(BCLog::error,BCLog::error,"BCEfficiencyFitter::SetHistograms() : Histogram 1 has fewer entries than histogram 2.");
			return 0;
		}
	}

	// set pointer to histograms
	fHistogram1 = hist1;
	fHistogram2 = hist2;
	
	// create efficiency histogram
	if (fHistogramRatio)
		delete fHistogramRatio;

	fHistogramRatio = new TGraphAsymmErrors();
	fHistogramRatio -> SetMarkerStyle(20);
	fHistogramRatio -> SetMarkerSize(1.5);

	int npoints = 0;

	// set points
	for (int i = 1; i <= hist1 -> GetNbinsX(); ++i)
	{
		// calculate uncertainties
		double xmin;
		double xmax;
		int flag = this -> GetUncertainties(
				int(hist1 -> GetBinContent(i)),
				int(hist2 -> GetBinContent(i)),
				0.68, xmin, xmax);

		if (flag == 1)
		{
			fHistogramRatio -> SetPoint(
					npoints,
					hist1 -> GetBinCenter(i),
					hist2 -> GetBinContent(i) / hist1 -> GetBinContent(i));
			// set uncertainties
			fHistogramRatio -> SetPointEXhigh(npoints, 0.);
			fHistogramRatio -> SetPointEXlow(npoints, 0.);
			fHistogramRatio -> SetPointEYhigh(npoints, xmax - hist2 -> GetBinContent(i) / hist1 -> GetBinContent(i));
			fHistogramRatio -> SetPointEYlow(npoints++, hist2 -> GetBinContent(i) / hist1 -> GetBinContent(i) - xmin);
		}
		else if (flag == -2)
		{
			fHistogramRatio -> SetPoint(npoints, hist1 -> GetBinCenter(i), 0.);
			// set uncertainties
			fHistogramRatio -> SetPointEXhigh(npoints, 0.);
			fHistogramRatio -> SetPointEXlow(npoints, 0.);
			fHistogramRatio -> SetPointEYhigh(npoints, xmax);
			fHistogramRatio -> SetPointEYlow(npoints++, 0.);
		}
		else if (flag == -1)
		{
			fHistogramRatio -> SetPoint(npoints, hist1 -> GetBinCenter(i), 1.);
			// set uncertainties
			fHistogramRatio -> SetPointEXhigh(npoints, 0.);
			fHistogramRatio -> SetPointEXlow(npoints, 0.);
			fHistogramRatio -> SetPointEYhigh(npoints, 0.);
			fHistogramRatio -> SetPointEYlow(npoints++, 1.-xmin);
		}
	}

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

	// calculate the lower and upper edge in x.
	double xmin = hist1 -> GetBinLowEdge(1);
	double xmax = hist1 -> GetBinLowEdge(nbins+1);

//	// calculate the minimum and maximum range in y.
//	double histymin = hist2 -> GetMinimum();
//	double histymax = hist1 -> GetMaximum();

//	// calculate the minimum and maximum of the function value based on
//	// the minimum and maximum value in y.
//	double ymin = TMath::Max(0., histymin - 5.*sqrt(histymin));
//	double ymax = histymax + 5.*sqrt(histymax);

	// set the data boundaries for x and y values.
	this -> SetDataBoundaries(0, xmin, xmax);
	this -> SetDataBoundaries(1, 0.0, 1.0);

	// set the indeces for fitting.
	this -> SetFitFunctionIndices(0, 1);

	// no error
	return 1;
}

// ---------------------------------------------------------

int BCEfficiencyFitter::SetFitFunction(TF1 * func)
{
	// check if function exists
	if(!func)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCEfficiencyFitter::SetFitFunction() : TF1 not created.");
		return 0;
	}

	// set the function
	fFitFunction = func;

	// update the model name to contain the function name
	SetName(TString::Format("BCEfficiencyFitter with %s",fFitFunction->GetName()));

	// reset parameters
	fParameterSet->clear();

	// get the new number of parameters
	int n = func->GetNpar();

	// add parameters
	for (int i = 0; i < n; ++i)
	{
		double xmin;
		double xmax;

		func->GetParLimits(i, xmin, xmax);

		AddParameter(func->GetParName(i), xmin, xmax);
	}

	// set flat prior for all parameters by default
	SetPriorConstantAll();

	return GetNParameters();;
}

// ---------------------------------------------------------

BCEfficiencyFitter::~BCEfficiencyFitter()
{
	if (fHistogram1)
		delete fHistogram1;

	if (fHistogram2)
		delete fHistogram2;

	if (fHistogramBinomial)
		delete fHistogramBinomial;

	if (fHistogramRatio)
		delete fHistogramRatio;
}

// ---------------------------------------------------------

/*
double BCEfficiencyFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	// using flat probability in all parameters
	double logprob = 0.;
	for(unsigned int i=0; i < this -> GetNParameters(); i++)
		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}
*/


// ---------------------------------------------------------

double BCEfficiencyFitter::LogLikelihood(std::vector <double> params)
{

	// initialize probability
	double loglikelihood = 0;

	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	// get the number of bins
	int nbins = fHistogram1 -> GetNbinsX();

	// get bin width
	double dx = fHistogram1 -> GetXaxis() -> GetBinWidth(0);

	// loop over all bins
	for (int ibin = 1; ibin <= nbins; ++ibin)
	{
		// get n
		int n = int(fHistogram1 -> GetBinContent(ibin));

		// get k
		int k = int(fHistogram2 -> GetBinContent(ibin));

		// get x
		double x = fHistogram1 -> GetBinCenter(ibin);

		double eff = 0;

		// use ROOT's TH1D::Integral method
		if (fFlagIntegration)
			eff = fFitFunction -> Integral(x - dx/2., x + dx/2.) / dx;

		// use linear interpolation
		else
			eff = (fFitFunction -> Eval(x + dx/2.) + fFitFunction -> Eval(x - dx/2.)) / 2.;

		// get the value of the Poisson distribution
		loglikelihood += BCMath::LogApproxBinomial(n, k, eff);
	}

	return loglikelihood;
}

// ---------------------------------------------------------

double BCEfficiencyFitter::FitFunction(std::vector <double> x, std::vector <double> params)
{
	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	return fFitFunction -> Eval(x[0]);
}

// ---------------------------------------------------------

int BCEfficiencyFitter::Fit(TH1D * hist1, TH1D * hist2, TF1 * func)
{
	// set histogram
	if (hist1 && hist2)
		this -> SetHistograms(hist1, hist2);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::Fit() : Histogram(s) not defined.");
		return 0;
	}

	// set function
	if (func)
		this -> SetFitFunction(func);
	else
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::Fit() : Fit function not defined.");
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
			fPValue = pvalue;
	else
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::Fit() : Could not use the fast p-value evaluation.");

	// print summary to screen
	this -> PrintShortFitSummary();

	// no error
	return 1;
}

// ---------------------------------------------------------

void BCEfficiencyFitter::DrawFit(const char * options, bool flaglegend)
{
	if (!fHistogram1 || !fHistogram2 || !fHistogramRatio)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::DrawFit() : Histogram(s) not defined.");
		return;
	}

	if (!fFitFunction)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::DrawFit() : Fit function not defined.");
		return;
	}

	// check wheather options contain "same"
	TString opt = options;
	opt.ToLower();

	// if not same, draw the histogram first to get the axes
	if(!opt.Contains("same"))
	{
		// create new histogram
		TH2D * hist_axes = new TH2D("hist_axes",
				Form(";%s;ratio", fHistogram1 -> GetXaxis() -> GetTitle()),
				fHistogram1 -> GetNbinsX(),
				fHistogram1 -> GetXaxis() -> GetBinLowEdge(1),
				fHistogram1 -> GetXaxis() -> GetBinLowEdge(fHistogram1 -> GetNbinsX()+1),
				1, 0., 1.);
		hist_axes -> SetStats(kFALSE);
		hist_axes -> Draw();

		fHistogramRatio -> Draw(TString::Format("%sp",opt.Data()));
	}

	// draw the error band as central 68% probability interval
	fErrorBand = this -> GetErrorBandGraph(0.16, 0.84);
	fErrorBand -> Draw("f same");

	// now draw the histogram again since it was covered by the band
	fHistogramRatio -> Draw(TString::Format("%ssamep",opt.Data()));

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
		legend -> AddEntry(fHistogramRatio, "Data", "L");
		legend -> AddEntry(fGraphFitFunction, "Best fit", "L");
		legend -> AddEntry(fErrorBand, "Error band", "F");
		legend -> Draw();
	}
	gPad -> RedrawAxis();
}

// ---------------------------------------------------------
int BCEfficiencyFitter::CalculatePValueFast(std::vector<double> par, double &pvalue)
{
	// check size of parameter vector
	if (par.size() != this -> GetNParameters())
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::CalculatePValueFast() : Number of parameters is inconsistent.");
		return 0;
	}

	// check if histogram exists
	if (!fHistogram1 || !fHistogram2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCEfficiencyFitter::CalculatePValueFast() : Histogram not defined.");
		return 0;
	}

	// define temporary variables
	int nbins = fHistogram1 -> GetNbinsX();

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
		double xmin = fHistogram1 -> GetBinLowEdge(ibin+1);
		double xmax = fHistogram1 -> GetBinLowEdge(ibin+2);

		// get the number of expected events
		double yexp = fFitFunction -> Integral(xmin, xmax);

		// get n
		int n = int(fHistogram1 -> GetBinContent(ibin));

		// get k
		int k = int(fHistogram2 -> GetBinContent(ibin));

		histogram[ibin]   = k;
		expectation[ibin] = n * yexp;

		// calculate p;
		logp += BCMath::LogApproxBinomial(n, k, yexp);
	}
	logp_start = logp;

	int niter = 100000;

	// loop over iterations
	for (int iiter = 0; iiter < niter; ++iiter)
	{
		// loop over bins
		for (int ibin = 0; ibin < nbins; ++ibin)
		{
			// get n
			int n = int(fHistogram1 -> GetBinContent(ibin));

			// get k
			int k = histogram.at(ibin);

			// get expectation
			double yexp = 0;
			if (n > 0)
				yexp = expectation.at(ibin)/n;

			// random step up or down in statistics for this bin
			double ptest = fRandom -> Rndm() - 0.5;

			// continue, if efficiency is at limit
			if (!(yexp > 0. || yexp < 1.0))
				continue;

			// increase statistics by 1
			if (ptest > 0 && (k < n))
			{
				// calculate factor of probability
				double r = (double(n)-double(k))/(double(k)+1.) * yexp / (1. - yexp);

				// walk, or don't (this is the Metropolis part)
				if (fRandom -> Rndm() < r)
				{
					histogram[ibin] = k + 1;
					logp += log(r);
				}
			}

			// decrease statistics by 1
			else if (ptest <= 0 && (k > 0))
			{
				// calculate factor of probability
				double r = double(k) / (double(n)-(double(k)-1)) * (1-yexp)/yexp;

				// walk, or don't (this is the Metropolis part)
				if (fRandom -> Rndm() < r)
				{
					histogram[ibin] = k - 1;
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
int BCEfficiencyFitter::GetUncertainties(int n, int k, double p, double &xmin, double &xmax)
{
	// create a histogram with binomial distribution
	if (fHistogramBinomial)
		fHistogramBinomial -> Reset();
	else
		fHistogramBinomial = new TH1D("hist_binomial", "", 1000, 0., 1.);

	// loop over bins and fill histogram
	for (int i = 1; i <= 1000; ++i)
	{
		double x   = fHistogramBinomial -> GetBinCenter(i);
		double val = TMath::Binomial(n, k) * TMath::Power(x, double(k)) * TMath::Power(1-x, double(n-k));
		fHistogramBinomial -> SetBinContent(i, val);
	}

	// normalize
	fHistogramBinomial -> Scale(1.0 / fHistogramBinomial -> Integral());

	// calculate quantiles
	int nprobSum = 4;
	double q[4];
	double probSum[4];
	probSum[0] = (1.-p)/2.;
	probSum[1] = 1.-(1.-p)/2.;
	probSum[2] = 0.05;
	probSum[3] = 0.95;

	fHistogramBinomial -> GetQuantiles(nprobSum, q, probSum);

	double xexp = double(k)/double(n);

	// calculate uncertainties
	if (n == 0)
	{
		xmin = 0.0;
		xmax = 0.0;
		return -3;
	}
	else if (xexp < q[0])
	{
		xmin = 0;
		xmax = q[3];
		return -2;
	}

	else if (xexp > q[1])
	{
		xmin = q[2];
		xmax = 1.0;
		return -1;
	}
	else
	{
		xmin = q[0];
		xmax = q[1];
		return 1;
	}

}

// ---------------------------------------------------------
