/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include <TGraphErrors.h>
#include <TF1.h>
#include <TString.h>
#include <TPad.h>
#include <TLegend.h> 
#include <TLegendEntry.h> 

#include "BCLog.h"
#include "BCDataSet.h"
#include "BCDataPoint.h"
#include "BCMath.h"

#include "BCGraphFitter.h"

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter() : BCModel("GraphFitter")
{
	fGraph = 0;
	fFitFunction = 0;
	fErrorBand = 0; 
	fGraphFitFunction = 0; 

	this -> MCMCSetNIterationsRun(2000);

	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

BCGraphFitter::BCGraphFitter(TGraphErrors * graph, TF1 * func) : BCModel("GraphFitter")
{
	fGraph = 0;
	fFitFunction = 0;
	fErrorBand = 0; 
	fGraphFitFunction = 0; 

	this -> MCMCSetNIterationsRun(2000);

	this -> SetGraph(graph);
	this -> SetFitFunction(func);

	this -> SetFillErrorBand();
}

// ---------------------------------------------------------

int BCGraphFitter::SetGraph(TGraphErrors * graph)
{
	if(!graph)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors not created.");
		return 0;
	}

	int npoints = graph -> GetN();
	if(!npoints)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors is empty.");
		return 0;
	}
	else if(npoints==1)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors has only one point. Not able to fit.");
		return 0;
	}
	else if(npoints==2)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors has only two points. Fit has no meaning.");
		return 0;
	}

	fGraph = graph;

	double * x = fGraph -> GetX();
	double * y = fGraph -> GetY();
	double * ex = fGraph -> GetEX();
	double * ey = fGraph -> GetEY();

	if(!ey)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetGraph() : TGraphErrors has no errors set on Y. Not able to fit.");
		return 0;
	}

	BCDataSet * ds = new BCDataSet();

	// fill the dataset
	// find x and y boundaries for the error band calculation
	double xmin=x[0];
	double xmax=x[0];
	double ymin=y[0];
	double ymax=y[0];
	for (int i = 0; i < npoints; ++i)
	{
		// if x errors are not set, set them to zero
		double errx = ex ? ex[i] : 0.;

		// create the data point
		BCDataPoint * dp = new BCDataPoint(4);
		dp -> SetValue(0, x[i]);
		dp -> SetValue(1, y[i]);
		dp -> SetValue(2, errx);
		dp -> SetValue(3, ey[i]);
		ds -> AddDataPoint(dp);

		if(x[i]-errx < xmin)
			xmin = x[i]-errx;
		else if(x[i]+errx > xmax)
			xmax = x[i]+errx;

		if(y[i] - 5.*ey[i] < ymin)
			ymin = y[i] - 5.*ey[i];
		else if(y[i] + 5.*ey[i] > ymax)
			ymax = y[i] + 5.*ey[i];
	}

	this -> SetDataSet(ds);

	// set boundaries for the error band calculation
	this -> SetDataBoundaries(0, xmin, xmax);
	this -> SetDataBoundaries(1, ymin, ymax);

	this -> SetFitFunctionIndices(0, 1);

	return this -> GetNDataPoints();
}

// ---------------------------------------------------------

int BCGraphFitter::SetFitFunction(TF1 * func)
{
	if(!func)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetFitFunction() : TF1 not created.");
		return 0;
	}

	// get the new number of parameters
	int npar = func -> GetNpar();
	if(!npar)
	{
		BCLog::Out(BCLog::error,BCLog::error,"BCGraphFitter::SetFitFunction() : TF1 has zero parameters. Not able to fit.");
		return 0;
	}

	// set the function
	fFitFunction = func;

	// update the model name to contain the function name
	this -> SetName(TString::Format("GraphFitter with %s",fFitFunction->GetName()));

	// reset parameters
	fParameterSet -> clear();

	// add parameters
	for (int i = 0; i < npar; ++i)
	{
		double xmin;
		double xmax;
		fFitFunction -> GetParLimits(i, xmin, xmax);

		this -> AddParameter(fFitFunction->GetParName(i), xmin, xmax);
	}

	return this -> GetNParameters();
}

// ---------------------------------------------------------

BCGraphFitter::~BCGraphFitter()
{}

// ---------------------------------------------------------

double BCGraphFitter::LogAPrioriProbability(std::vector <double> parameters)
{
	// using flat probability in all parameters
	return 0.;
}

// ---------------------------------------------------------

double BCGraphFitter::LogLikelihood(std::vector <double> params)
{
	// initialize probability
	double logl = 0.;

	// set the parameters of the function
	// passing the pointer to first element of the vector is
	// not completely safe as there might be an implementation where
	// the vector elements are not stored consecutively in memory.
	// however it is much faster than copying the contents, especially
	// for large number of parameters
	fFitFunction -> SetParameters(&params[0]);

	// loop over all data points
	for (int i = 0; i < this -> GetNDataPoints(); i++)
	{
		std::vector <double> x = GetDataPoint(i) -> GetValues();

		// her we ignore the errors on x even when they're available
		// i.e. we treat them just as the region specifiers
		double y = x[1];
		double yerr = x[3];
		double yexp = this -> FitFunction(x,params);

		// calculate log of probability assuming
		// a Gaussian distribution for each point
		logl += BCMath::LogGaus(y, yexp, yerr, true);
	}

	return logl;
}

// ---------------------------------------------------------

double BCGraphFitter::FitFunction(std::vector <double> x, std::vector <double> params)
{
	// set the parameters of the function
	fFitFunction -> SetParameters(&params[0]);

	return fFitFunction -> Eval(x[0]);
}

// ---------------------------------------------------------

int BCGraphFitter::Fit(TGraphErrors * graph, TF1 * func)
{
	// set graph
	if (!this -> SetGraph(graph))
		return 0;

	// set function
	if (!this -> SetFitFunction(func))
		return 0;

	// perform marginalization
	this -> MarginalizeAll();

	// maximize posterior probability, using the best-fit values close
	// to the global maximum from the MCMC
	this -> FindModeMinuit( this -> GetBestFitParameters(), -1);

	// calculate p-value from the chi2 probability
	// this is only valid for a product of gaussiang which is the case for
	// the BCGraphFitter
	this -> GetPvalueFromChi2(this -> GetBestFitParameters(), 3);

	// print summary to screen
	this -> PrintFitSummary();

	return 1;
}

// ---------------------------------------------------------

void BCGraphFitter::DrawFit(const char * options)
{
	if (!fGraph)
	{
		BCLog::Out(BCLog::error, BCLog::error,"BCGraphFitter::DrawFit() : TGraphErrors not defined.");
		return;
	}

	if (!fFitFunction)
	{
		BCLog::Out(BCLog::error, BCLog::error,"BCGraphFitter::DrawFit() : Fit function not defined.");
		return;
	}

	// check wheather options contain "same"
	TString opt = options;
	opt.ToLower();

	// if not same, draw the histogram first to get the axes
	if(!opt.Contains("same"))
		fGraph -> Draw("ap");

	// draw the error band as central 68% probability interval
	fErrorBand = this -> GetErrorBandGraph(0.16, 0.84);
	fErrorBand -> Draw("f same");

	// draw the fit function on top
	fGraphFitFunction = this -> GetFitFunctionGraph( this->GetBestFitParameters() );
	fGraphFitFunction -> SetLineColor(kRed);
	fGraphFitFunction -> SetLineWidth(2);
	fGraphFitFunction -> Draw("l same");

	// now draw the histogram again since it was covered by the band and
	// the best fit
	fGraph -> Draw("p same");

	// draw legend
	if (opt.Contains("leg"))
		{
			TLegend * legend = new TLegend(0.25, 0.75, 0.55, 0.95); 
			legend -> SetBorderSize(0); 
			legend -> SetFillColor(kWhite); 
			legend -> AddEntry(fGraph, "Data", "P"); 
			legend -> AddEntry(fGraphFitFunction, "Best fit", "L"); 
			legend -> AddEntry(fErrorBand, "Error band", "F"); 
			legend -> Draw(); 
		}

	gPad -> RedrawAxis();
}

// ---------------------------------------------------------

void BCGraphFitter::PrintFitSummary()
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


// std::cout << std::endl;
// 	std::cout << "Fit summary " << std::endl;
// 	std::cout << "------------------------------------ " << std::endl;

// 	std::cout
// 			<< "Number of parameters : "
// 			<< this -> GetNParameters() << std::endl;
// 	std::cout << std::endl;

// 	std::cout << "Best fit parameters (global) : " << std::endl;
// 	for (unsigned int i = 0; i < this -> GetNParameters(); ++i)
// 		std::cout
// 				<< this -> GetParameter(i) -> GetName() << " : "
// 				<< this -> GetBestFitParameter(i) << std::endl;
// 	std::cout << std::endl;

// 	std::cout << "Goodness-of-fit test : " << std::endl;
// 	std::cout << " p-value = " << this -> GetPValue() << std::endl; 

//	std::cout << "Best fit parameters (marginalized) : " << std::endl;
//	for (int i = 0; i < this -> GetNParameters(); ++i)
//	{
//		BCH1D * bch1d = this -> GetMarginalized(fParameterSet -> at(i));
//		std::cout
//				<< this -> GetParameter(i) -> GetName() << " : "
//				<< this -> GetBestFitParameterMarginalized(i) << std::endl;
//	}
//	std::cout << std::endl;
}

// ---------------------------------------------------------
