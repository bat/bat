/*
 * Copyright (C) 2008, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "config.h"

#include "BAT/BCIntegrate.h"
#include "BAT/BCLog.h"
#include "BAT/BCMath.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TTree.h>

#ifdef HAVE_CUBA_H
#include "cuba.h"
#endif

// ---------------------------------------------------------

class BCIntegrate * global_this;

// *********************************************
BCIntegrate::BCIntegrate() : BCEngineMCMC()
{
	fNvar=0;
	fNiterPerDimension = 100;
	fNSamplesPer2DBin = 100;
	fRandom = new TRandom3(0);

	fNIterationsMax    = 1000000;
	fRelativePrecision = 1e-3;

#ifdef HAVE_CUBA_H
	fIntegrationMethod = BCIntegrate::kIntCuba;
#else
	fIntegrationMethod = BCIntegrate::kIntMonteCarlo;
#endif
	fMarginalizationMethod = BCIntegrate::kMargMetropolis;
	fOptimizationMethod = BCIntegrate::kOptMinuit;

	fNbins=100;

	fDataPointLowerBoundaries = 0;
	fDataPointUpperBoundaries = 0;

	fFillErrorBand = false;

	fFitFunctionIndexX = -1;
	fFitFunctionIndexY = -1;

	fErrorBandNbinsX = 100;
	fErrorBandNbinsY = 500;

	fErrorBandContinuous = true;

	fMinuit = 0;
	fMinuitArglist[0] = 20000;
	fMinuitArglist[1] = 0.01;

	fFlagWriteMarkovChain = false;
	fMarkovChainTree = 0;

	fMarkovChainStepSize = 0.1;

	fMarkovChainAutoN = true;

	fSAT0 = 100.0;
	fSATmin = 0.1;
	fSASchedule = BCIntegrate::kSACauchy;
}

// *********************************************
BCIntegrate::BCIntegrate(BCParameterSet * par) : BCEngineMCMC(1)
{
	fNvar=0;
	fNiterPerDimension = 100;
	fNSamplesPer2DBin = 100;
	fRandom = new TRandom3(0);

	fNbins=100;

	this->SetParameters(par);

	fDataPointLowerBoundaries = 0;
	fDataPointUpperBoundaries = 0;

	fFillErrorBand = false;

	fFitFunctionIndexX = -1;
	fFitFunctionIndexY = -1;

	fErrorBandNbinsX = 100;
	fErrorBandNbinsY = 200;

	fErrorBandContinuous = true;

	fMinuit = 0;
	fMinuitArglist[0] = 20000;
	fMinuitArglist[1] = 0.01;

	fFlagWriteMarkovChain = false;
	fMarkovChainTree = 0;

	fMarkovChainStepSize = 0.1;

	fMarkovChainAutoN = true;
	
	fSAT0 = 100.0;
	fSATmin = 0.1;
}

// *********************************************
BCIntegrate::~BCIntegrate()
{
	DeleteVarList();

	fx=0;

	delete fRandom;
	fRandom=0;

	if (fMinuit)
		delete fMinuit;

	int n1 = fHProb1D.size();
	if(n1>0)
	{
		for (int i=0;i<n1;i++)
			delete fHProb1D.at(i);
	}

	int n2 = fHProb2D.size();
	if(n2>0)
	{
		for (int i=0;i<n2;i++)
			delete fHProb2D.at(i);
	}
}

// *********************************************
void BCIntegrate::SetParameters(BCParameterSet * par)
{
	DeleteVarList();

	fx = par;
	fNvar = fx->size();
	fMin = new double[fNvar];
	fMax = new double[fNvar];
	fVarlist = new int[fNvar];

	this->ResetVarlist(1);

	for(int i=0;i<fNvar;i++)
	{
		fMin[i]=(fx->at(i))->GetLowerLimit();
		fMax[i]=(fx->at(i))->GetUpperLimit();
	}

	fXmetro0.clear();
	for(int i=0;i<fNvar;i++)
		fXmetro0.push_back((fMin[i]+fMax[i])/2.0);

	fXmetro1.clear();
	fXmetro1.assign(fNvar,0.);

	fMCMCBoundaryMin.clear();
	fMCMCBoundaryMax.clear();

	for(int i=0;i<fNvar;i++)
	{
		fMCMCBoundaryMin.push_back(fMin[i]);
		fMCMCBoundaryMax.push_back(fMax[i]);
	}

	for (int i = int(fMCMCH1NBins.size()); i<fNvar; ++i)
		fMCMCH1NBins.push_back(100);

	fMCMCNParameters = fNvar;
}

// *********************************************
void BCIntegrate::SetMarkovChainInitialPosition(std::vector<double> position)
{
	int n = position.size();

	fXmetro0.clear();

	for (int i = 0; i < n; ++i)
		fXmetro0.push_back(position.at(i));
}

// *********************************************
void BCIntegrate::DeleteVarList()
{
	if(fNvar)
	{
		delete[] fVarlist;
		fVarlist=0;

		delete[] fMin;
		fMin=0;

		delete[] fMax;
		fMax=0;

		fx=0;
		fNvar=0;
	}
}

// *********************************************
void BCIntegrate::SetNbins(int nbins, int index)
{
	if (fNvar == 0)
		return;

	// check if index is in range
	if (index >= fNvar)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,"BCIntegrate::SetNbins : Index out of range.");
		return;
	}
	// set for all parameters at once
	else if (index < 0)
	{
		for (int i = 0; i < fNvar; ++i)
			this -> SetNbins(nbins, i);
		return;
	}

	// sanity check for nbins
	if (nbins <= 0)
		nbins = 10;

	fMCMCH1NBins[index] = nbins;

	return;

// 	if(n<2)
// 	{
// 		BCLog::Out(BCLog::warning, BCLog::warning,"BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
// 		n=2;
// 	}
// 	this -> MCMCSetH1NBins(n, -1);

	//	fNbins=n;

	//	fMCMCH1NBins = n;
	//	fMCMCH2NBinsX = n;
	//	fMCMCH2NBinsY = n;
}

// *********************************************
// void BCIntegrate::SetNbinsX(int n)
// {
// 	if(n<2)
// 	{
// 		BCLog::Out(BCLog::warning, BCLog::warning,"BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
// 		n=2;
// 	}
// 	fMCMCH2NBinsX = n;
// }

// *********************************************
// void BCIntegrate::SetNbinsY(int n)
// {
// 	if(n<2)
// 	{
// 		BCLog::Out(BCLog::warning, BCLog::warning,"BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
// 		n=2;
// 	}
// 	fNbins=n;

// 	fMCMCH2NBinsY = n;
// }

// *********************************************
void BCIntegrate::SetVarList(int * varlist)
{
	for(int i=0;i<fNvar;i++)
		fVarlist[i]=varlist[i];
}

// *********************************************
void BCIntegrate::ResetVarlist(int v)
{
	for(int i=0;i<fNvar;i++)
		fVarlist[i]=v;
}

// *********************************************
double BCIntegrate::Eval(std::vector <double> x)
{
	BCLog::Out(BCLog::warning, BCLog::warning, "BCIntegrate::Eval. No function. Function needs to be overloaded.");
	return 0;
}

// *********************************************
double BCIntegrate::LogEval(std::vector <double> x)
{
	// this method should better also be overloaded
	return log(this->Eval(x));
}

// *********************************************
double BCIntegrate::EvalSampling(std::vector <double> x)
{
	BCLog::Out(BCLog::warning, BCLog::warning, "BCIntegrate::EvalSampling. No function. Function needs to be overloaded.");
	return 0;
}

// *********************************************
double BCIntegrate::LogEvalSampling(std::vector <double> x)
{
	return log(this->EvalSampling(x));
}

// *********************************************
double BCIntegrate::EvalPrint(std::vector <double> x)
{
	double val=this->Eval(x);
	BCLog::Out(BCLog::detail, BCLog::detail, Form("BCIntegrate::EvalPrint. Value: %d.", val));

	return val;
}

// *********************************************
void BCIntegrate::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
#ifdef HAVE_CUBA_H
	fIntegrationMethod = method;
#else
	BCLog::Out(BCLog::warning,BCLog::warning,"!!! This version of BAT is compiled without Cuba.");
	BCLog::Out(BCLog::warning,BCLog::warning,"    Monte Carlo Sampled Mean integration method will be used.");
	BCLog::Out(BCLog::warning,BCLog::warning,"    To be able to use Cuba integration, install Cuba and recompile BAT.");
#endif
}

// *********************************************
double BCIntegrate::Integrate()
{
	std::vector <double> parameter;
	parameter.assign(fNvar, 0.0);

	switch(fIntegrationMethod)
	{
		case BCIntegrate::kIntMonteCarlo:
			return IntegralMC(parameter);

		case BCIntegrate::kIntMetropolis:
			return this -> IntegralMetro(parameter);

		case BCIntegrate::kIntImportance:
			return this -> IntegralImportance(parameter);

		case BCIntegrate::kIntCuba:
#ifdef HAVE_CUBA_H
			return this -> CubaIntegrate();
#else
			BCLog::Out(BCLog::error,BCLog::error,"!!! This version of BAT is compiled without Cuba.");
			BCLog::Out(BCLog::error,BCLog::error,"    Use other integration methods or install Cuba and recompile BAT.");
#endif
	}

	BCLog::Out(BCLog::error, BCLog::error,
			Form("BCIntegrate::Integrate : Invalid integration method: %d", fIntegrationMethod));

	return 0;
}

// *********************************************
double BCIntegrate::IntegralMC(std::vector <double> x, int * varlist)
{
	this->SetVarList(varlist);
	return IntegralMC(x);
}

// *********************************************
double BCIntegrate::IntegralMC(std::vector <double> x)
{
	// count the variables to integrate over
	int NvarNow=0;

	for(int i = 0; i < fNvar; i++)
		if(fVarlist[i])
			NvarNow++;

	// print to log
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Running MC integation over %i dimensions.", NvarNow));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Maximum number of iterations: %i", this->GetNIterationsMax()));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Aimed relative precision:     %e", this->GetRelativePrecision()));

	// calculate D (the integrated volume)
	double D = 1.0;
	for(int i = 0; i < fNvar; i++)
		if (fVarlist[i])
			D = D * (fMax[i] - fMin[i]);

	// reset variables
	double pmax = 0.0;
	double sumW  = 0.0;
	double sumW2 = 0.0;
	double integral = 0.0;
	double variance = 0.0;
	double relprecision = 1.0;

	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	// reset number of iterations
	fNIterations = 0;

	// iterate while precision is not reached and the number of iterations is lower than maximum number of iterations
	while ((fRelativePrecision < relprecision && fNIterationsMax > fNIterations) || fNIterations < 10)
	{
		// increase number of iterations
		fNIterations++;

		// get random numbers
		this -> GetRandomVector(randx);

		// scale random numbers
		for(int i = 0; i < fNvar; i++)
		{
			if(fVarlist[i])
				randx[i]=fMin[i]+randx[i]*(fMax[i]-fMin[i]);
			else
				randx[i]=x[i];
		}

		// evaluate function at sampled point
		double value = this->Eval(randx);

		// add value to sum and sum of squares
		sumW  += value;
		sumW2 += value * value;

		// search for maximum probability
		if (value > pmax)
		{
			// set new maximum value
			pmax = value;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// calculate integral and variance
		integral = D * sumW / fNIterations;

		if (fNIterations%10000 == 0)
		{
			variance = (1.0 / double(fNIterations)) * (D * D * sumW2 / double(fNIterations) - integral * integral);
			double error = sqrt(variance);
			relprecision = error / integral;

			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralMC. Iteration %i, integral: %e +- %e.", fNIterations, integral, error));
		}
	}

	fError = variance / fNIterations;

	// print to log
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Result of integration:        %e +- %e   in %i iterations.", integral, sqrt(variance), fNIterations));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Obtained relative precision:  %e. ", sqrt(variance) / integral));

	return integral;
}


// *********************************************
double BCIntegrate::IntegralMetro(std::vector <double> x)
{
	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralMetro. Integate over %i dimensions.", fNvar));

	// get total number of iterations
	double Niter = pow(fNiterPerDimension, fNvar);

	// print if more than 100,000 iterations
	if(Niter>1e5)
		BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Total number of iterations in Metropolis: %d.", Niter));

	// reset sum
	double sumI = 0;

	// prepare Metropolis algorithm
	std::vector <double> randx;
	randx.assign(fNvar,0.);
	InitMetro();

	// prepare maximum value
	double pmax = 0.0;

	// loop over iterations
	for(int i = 0; i <= Niter; i++)
	{
		// get random point from sampling distribution
		this -> GetRandomPointSamplingMetro(randx);

		// calculate probability at random point
		double val_f = this -> Eval(randx);

		// calculate sampling distributions at that point
		double val_g = this -> EvalSampling(randx);

		// add ratio to sum
		if (val_g > 0)
			sumI += val_f / val_g;

		// search for maximum probability
		if (val_f > pmax)
		{
			// set new maximum value
			pmax = val_f;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// write intermediate results
		if((i+1)%100000 == 0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralMetro. Iteration %d, integral: %d.", i+1, sumI/(i+1)));
	}

	// calculate integral
	double result=sumI/Niter;

	// print debug information
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Integral is %f in %i iterations.", result, Niter));

	return result;
}

// *********************************************
double BCIntegrate::IntegralImportance(std::vector <double> x)
{
	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralImportance. Integate over %i dimensions.", fNvar));

	// get total number of iterations
	double Niter = pow(fNiterPerDimension, fNvar);

	// print if more than 100,000 iterations
	if(Niter>1e5)
		BCLog::Out(BCLog::detail, BCLog::detail, Form("BCIntegrate::IntegralImportance. Total number of iterations: %d.", Niter));

	// reset sum
	double sumI = 0;

	std::vector <double> randx;
	randx.assign(fNvar,0.);

	// prepare maximum value
	double pmax = 0.0;

	// loop over iterations
	for(int i = 0; i <= Niter; i++)
	{
		// get random point from sampling distribution
		this -> GetRandomPointImportance(randx);

		// calculate probability at random point
		double val_f = this -> Eval(randx);

		// calculate sampling distributions at that point
		double val_g = this -> EvalSampling(randx);

		// add ratio to sum
		if (val_g > 0)
			sumI += val_f / val_g;

		// search for maximum probability
		if (val_f > pmax)
		{
			// set new maximum value
			pmax = val_f;

			// delete old best fit parameter values
			fBestFitParameters.clear();

			// write best fit parameters
			for(int i = 0; i < fNvar; i++)
				fBestFitParameters.push_back(randx.at(i));
		}

		// write intermediate results
		if((i+1)%100000 == 0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::IntegralImportance. Iteration %d, integral: %d.", i+1, sumI/(i+1)));
	}

	// calculate integral
	double result=sumI/Niter;

	// print debug information
	BCLog::Out(BCLog::debug, BCLog::debug, Form("BCIntegrate::IntegralImportance. Integral %f in %i iterations.", result, Niter));

	return result;
}

// *********************************************
TH1D* BCIntegrate::Marginalize(BCParameter * parameter)
{
	BCLog::Out(BCLog::detail, BCLog::detail,
						 Form(" --> Marginalizing model wrt. parameter %s using method %d.", parameter->GetName().data(), fMarginalizationMethod));

	switch(fMarginalizationMethod)
	{
		case BCIntegrate::kMargMonteCarlo:
			return MarginalizeByIntegrate(parameter);

		case BCIntegrate::kMargMetropolis:
			return MarginalizeByMetro(parameter);
	}

	BCLog::Out(BCLog::warning, BCLog::warning,
		Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizationMethod));

	return 0;
}

// *********************************************
TH2D * BCIntegrate::Marginalize(BCParameter * parameter1, BCParameter * parameter2)
{
	switch(fMarginalizationMethod)
	{
		case BCIntegrate::kMargMonteCarlo:
			return MarginalizeByIntegrate(parameter1,parameter2);

		case BCIntegrate::kMargMetropolis:
			return MarginalizeByMetro(parameter1,parameter2);
	}

	BCLog::Out(BCLog::warning, BCLog::warning,
		Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizationMethod));

	return 0;
}

// *********************************************
TH1D* BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter)
{
	// set parameter to marginalize
	this->ResetVarlist(1);
	int index = parameter->GetIndex();
	this->UnsetVar(index);

	// define histogram
	double hmin = parameter -> GetLowerLimit();
	double hmax = parameter -> GetUpperLimit();
	double hdx  = (hmax - hmin) / double(fNbins);
	TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

	// integrate
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	for(int i=0;i<=fNbins;i++)
	{
		randx[index] = hmin + (double)i * hdx;

		double val = IntegralMC(randx);
		hist->Fill(randx[index], val);
	}

	// normalize
	hist -> Scale( 1./hist->Integral("width") );

	return hist;
}

// *********************************************
TH2D * BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter1, BCParameter * parameter2)
{
	// set parameter to marginalize
	this->ResetVarlist(1);
	int index1 = parameter1->GetIndex();
	this->UnsetVar(index1);
	int index2 = parameter2->GetIndex();
	this->UnsetVar(index2);

	// define histogram
	double hmin1 = parameter1 -> GetLowerLimit();
	double hmax1 = parameter1 -> GetUpperLimit();
	double hdx1  = (hmax1 - hmin1) / double(fNbins);

	double hmin2 = parameter2 -> GetLowerLimit();
	double hmax2 = parameter2 -> GetUpperLimit();
	double hdx2  = (hmax2 - hmin2) / double(fNbins);

	TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1 -> GetName().data(), parameter2 -> GetName().data()),"",
			fNbins, hmin1, hmax1,
			fNbins, hmin2, hmax2);

	// integrate
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);

	for(int i=0;i<=fNbins;i++)
	{
		randx[index1] = hmin1 + (double)i * hdx1;
		for(int j=0;j<=fNbins;j++)
		{
			randx[index2] = hmin2 + (double)j * hdx2;

			double val = IntegralMC(randx);
			hist->Fill(randx[index1],randx[index2], val);
		}
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist;
}

// *********************************************
TH1D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter)
{
	int niter = fMarkovChainNIterations;

	if (fMarkovChainAutoN == true)
		niter = fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	BCLog::Out(BCLog::detail, BCLog::detail,
		Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

	// set parameter to marginalize
	int index = parameter->GetIndex();

	// define histogram
	double hmin = parameter -> GetLowerLimit();
	double hmax = parameter -> GetUpperLimit();
	TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

	// prepare Metro
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);
	InitMetro();

	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);
		hist->Fill(randx[index]);
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist;
}

// *********************************************
TH2D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter1, BCParameter * parameter2)
{
	int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	// set parameter to marginalize
	int index1 = parameter1->GetIndex();
	int index2 = parameter2->GetIndex();

	// define histogram
	double hmin1 = parameter1 -> GetLowerLimit();
	double hmax1 = parameter1 -> GetUpperLimit();

	double hmin2 = parameter2 -> GetLowerLimit();
	double hmax2 = parameter2 -> GetUpperLimit();

	TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1 -> GetName().data(), parameter2 -> GetName().data()),"",
			fNbins, hmin1, hmax1,
			fNbins, hmin2, hmax2);

	// prepare Metro
	std::vector <double> randx;
	randx.assign(fNvar, 0.0);
	InitMetro();

	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);
		hist->Fill(randx[index1],randx[index2]);
	}

	// normalize
	hist -> Scale(1.0/hist->Integral("width"));

	return hist;
}

// *********************************************
int BCIntegrate::MarginalizeAllByMetro(const char * name="")
{
	int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

	BCLog::Out(BCLog::detail, BCLog::detail,
			Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

	// define 1D histograms
	for(int i=0;i<fNvar;i++)
	{
		double hmin1 = fx->at(i) -> GetLowerLimit();
		double hmax1 = fx->at(i) -> GetUpperLimit();

		TH1D * h1 = new TH1D(
				TString::Format("h%s_%s", name, fx->at(i) -> GetName().data()),"",
				fNbins, hmin1, hmax1);

		fHProb1D.push_back(h1);
	}

	// define 2D histograms
	for(int i=0;i<fNvar-1;i++)
		for(int j=i+1;j<fNvar;j++)
		{
			double hmin1 = fx->at(i) -> GetLowerLimit();
			double hmax1 = fx->at(i) -> GetUpperLimit();

			double hmin2 = fx->at(j) -> GetLowerLimit();
			double hmax2 = fx->at(j) -> GetUpperLimit();

			TH2D * h2 = new TH2D(
				TString::Format("h%s_%s_%s", name, fx->at(i) -> GetName().data(), fx->at(j) -> GetName().data()),"",
				fNbins, hmin1, hmax1,
				fNbins, hmin2, hmax2);

			fHProb2D.push_back(h2);
		}

	// get number of 2d distributions
	int nh2d = fHProb2D.size();

	BCLog::Out(BCLog::detail, BCLog::detail,
			Form(" --> Marginalizing %d 1D distributions and %d 2D distributions.", fNvar, nh2d));

	// prepare function fitting
	double dx = 0.0;
	double dy = 0.0;

	if (fFitFunctionIndexX >= 0)
	{
		dx = (fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexX) -
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexX))
				/ double(fErrorBandNbinsX);

		dx = (fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexY) -
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexY))
				/ double(fErrorBandNbinsY);

		fErrorBandXY = new TH2D(
				TString::Format("errorbandxy_%d",BCLog::GetHIndex()), "",
				fErrorBandNbinsX,
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexX) - 0.5 * dx,
				fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexX) + 0.5 * dx,
				fErrorBandNbinsY,
				fDataPointLowerBoundaries -> GetValue(fFitFunctionIndexY) - 0.5 * dy,
				fDataPointUpperBoundaries -> GetValue(fFitFunctionIndexY) + 0.5 * dy);
		fErrorBandXY -> SetStats(kFALSE);

		for (int ix = 1; ix <= fErrorBandNbinsX; ++ix)
			for (int iy = 1; iy <= fErrorBandNbinsX; ++iy)
				fErrorBandXY -> SetBinContent(ix, iy, 0.0);
	}

	// prepare Metro
	std::vector <double> randx;

	randx.assign(fNvar, 0.0);
	InitMetro();

	// reset counter
	fNacceptedMCMC = 0;

	// run Metro and fill histograms
	for(int i=0;i<=niter;i++)
	{
		GetRandomPointMetro(randx);

		// save this point to the markov chain in the ROOT file
		if (fFlagWriteMarkovChain)
			{
				fMCMCIteration = i;
				fMarkovChainTree -> Fill();
			}

		for(int j=0;j<fNvar;j++)
			fHProb1D[j] -> Fill( randx[j] );

		int ih=0;
		for(int j=0;j<fNvar-1;j++)
			for(int k=j+1;k<fNvar;k++)
			{
				fHProb2D[ih] -> Fill(randx[j],randx[k]);
				ih++;
			}

		if((i+1)%100000==0)
			BCLog::Out(BCLog::debug, BCLog::debug,
				Form("BCIntegrate::MarginalizeAllByMetro. %d samples done.",i+1));

		// function fitting

		if (fFitFunctionIndexX >= 0)
		{
			// loop over all possible x values ...
			if (fErrorBandContinuous)
			{
				double x = 0;

				for (int ix = 0; ix < fErrorBandNbinsX; ix++)
				{
					// calculate x
					x = fErrorBandXY -> GetXaxis() -> GetBinCenter(ix + 1);

					// calculate y
					std::vector <double> xvec;
					xvec.push_back(x);
					double y = this -> FitFunction(xvec, randx);
					xvec.clear();

					// fill histogram
					fErrorBandXY -> Fill(x, y);
				}
			}

			// ... or evaluate at the data point x-values
			else
			{
				int ndatapoints = int(fErrorBandX.size());
				double x = 0;

				for (int ix = 0; ix < ndatapoints; ++ix)
				{
					// calculate x
					x = fErrorBandX.at(ix);

					// calculate y
					std::vector <double> xvec;
					xvec.push_back(x);
					double y = this -> FitFunction(xvec, randx);
					xvec.clear();

					// fill histogram
					fErrorBandXY -> Fill(x, y);
				}
			}
		}
	}

	// normalize histograms
	for(int i=0;i<fNvar;i++)
		fHProb1D[i] -> Scale( 1./fHProb1D[i]->Integral("width") );

	for (int i=0;i<nh2d;i++)
		fHProb2D[i] -> Scale( 1./fHProb2D[i]->Integral("width") );

	if (fFitFunctionIndexX >= 0)
		fErrorBandXY -> Scale(1.0/fErrorBandXY -> Integral() * fErrorBandXY -> GetNbinsX());

	if (fFitFunctionIndexX >= 0)
	{
		for (int ix = 1; ix <= fErrorBandNbinsX; ix++)
		{
			double sum = 0;

			for (int iy = 1; iy <= fErrorBandNbinsY; iy++)
				sum += fErrorBandXY -> GetBinContent(ix, iy);

			for (int iy = 1; iy <= fErrorBandNbinsY; iy++)
			{
				double newvalue = fErrorBandXY -> GetBinContent(ix, iy) / sum;
				fErrorBandXY -> SetBinContent(ix, iy, newvalue);
			}
		}
	}

	BCLog::Out(BCLog::detail, BCLog::detail,
		   Form("BCIntegrate::MarginalizeAllByMetro done with %i trials and %i accepted trials. Efficiency is %f",fNmetro, fNacceptedMCMC, double(fNacceptedMCMC)/double(fNmetro)));

	return fNvar+nh2d;
}

// *********************************************
TH1D * BCIntegrate::GetH1D(int parIndex)
{
	if(fHProb1D.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCModel::GetH1D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	if(parIndex<0 || parIndex>fNvar)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH1D. Parameter index %d is invalid.",parIndex));
		return 0;
	}

	return fHProb1D[parIndex];
}

// *********************************************
int BCIntegrate::GetH2DIndex(int parIndex1, int parIndex2)
{
	if(parIndex1>fNvar-1 || parIndex1<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex1));
		return -1;
	}

	if(parIndex2>fNvar-1 || parIndex2<0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex2));
		return -1;
	}

	if(parIndex1==parIndex2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			Form("BCIntegrate::GetH2DIndex. Parameters have equal indeces: %d , %d", parIndex1, parIndex2));
		return -1;
	}

	if(parIndex1>parIndex2)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCIntegrate::GetH2DIndex. First parameters must be smaller than second (sorry).");
		return -1;
	}

	int index=0;
	for(int i=0;i<fNvar-1;i++)
		for(int j=i+1;j<fNvar;j++)
		{
			if(i==parIndex1 && j==parIndex2)
				return index;
			index++;
		}

	BCLog::Out(BCLog::warning, BCLog::warning,
		"BCIntegrate::GetH2DIndex. Invalid index combination.");

	return -1;
}

// *********************************************
TH2D * BCIntegrate::GetH2D(int parIndex1, int parIndex2)
{
	if(fHProb2D.size()==0)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCModel::GetH2D. MarginalizeAll() has to be run prior to this to fill the distributions.");
		return 0;
	}

	int hindex = this -> GetH2DIndex(parIndex1, parIndex2);
	if(hindex==-1)
		return 0;

	if(hindex>(int)fHProb2D.size()-1)
	{
		BCLog::Out(BCLog::warning, BCLog::warning,
			"BCIntegrate::GetH2D. Got invalid index from GetH2DIndex(). Something went wrong.");
		return 0;
	}

	return fHProb2D[hindex];
}

// *********************************************
double BCIntegrate::GetRandomPoint(std::vector <double> &x)
{
	this->GetRandomVector(x);

	for(int i=0;i<fNvar;i++)
		x[i]=fMin[i]+x[i]*(fMax[i]-fMin[i]);

	return this->Eval(x);
}

// *********************************************
double BCIntegrate::GetRandomPointImportance(std::vector <double> &x)
{
	double p = 1.1;
	double g = 1.0;

	while (p>g)
	{
		this->GetRandomVector(x);

		for(int i=0;i<fNvar;i++)
			x[i] = fMin[i] + x[i] * (fMax[i]-fMin[i]);

		p = fRandom->Rndm();

		g = this -> EvalSampling(x);
	}

	return this->Eval(x);
}

// *********************************************
void BCIntegrate::InitMetro()
{
	fNmetro=0;

	if (fXmetro0.size() <= 0)
	{
		// start in the center of the phase space
		for(int i=0;i<fNvar;i++)
			fXmetro0.push_back((fMin[i]+fMax[i])/2.0);
	}

	fMarkovChainValue = this ->  LogEval(fXmetro0);

	// run metropolis for a few times and dump the result... (to forget the initial position)
	std::vector <double> x;
	x.assign(fNvar,0.);

	for(int i=0;i<1000;i++)
		GetRandomPointMetro(x);

	fNmetro = 0;
}

// *********************************************
void BCIntegrate::GetRandomPointMetro(std::vector <double> &x)
{
	// get new point
	this->GetRandomVectorMetro(fXmetro1);

	// scale the point to the allowed region and stepsize
	int in=1;
	for(int i=0;i<fNvar;i++)
	{
		fXmetro1[i] = fXmetro0[i] + fXmetro1[i] * (fMax[i]-fMin[i]);

		// check whether the generated point is inside the allowed region
		if( fXmetro1[i]<fMin[i] || fXmetro1[i]>fMax[i] )
			in=0;
	}

	// calculate the log probabilities and compare old and new point

	double p0 = fMarkovChainValue; // old point
	double p1 = 0; // new point
	int accept=0;

	if(in)
	{
		p1 = this -> LogEval(fXmetro1);

		if(p1>=p0)
			accept=1;
		else
		{
			double r=log(fRandom->Rndm());
			if(r<p1-p0)
				accept=1;
		}
	}

	// fill the return point after the decision
	if(accept)
	{
		// increase counter
		fNacceptedMCMC++;

		for(int i=0;i<fNvar;i++)
		{
			fXmetro0[i]=fXmetro1[i];
			x[i]=fXmetro1[i];
			fMarkovChainValue = p1;
		}
	}
	else
		for(int i=0;i<fNvar;i++)
		{
			x[i]=fXmetro0[i];
			fXmetro1[i] = x[i];
			fMarkovChainValue = p0;
		}

	fNmetro++;
}

// *********************************************
void BCIntegrate::GetRandomPointSamplingMetro(std::vector <double> &x)
{
	// get new point
	this->GetRandomVectorMetro(fXmetro1);

	// scale the point to the allowed region and stepsize
	int in=1;
	for(int i=0;i<fNvar;i++)
	{
		fXmetro1[i] = fXmetro0[i] + fXmetro1[i] * (fMax[i]-fMin[i]);

		// check whether the generated point is inside the allowed region
		if( fXmetro1[i]<fMin[i] || fXmetro1[i]>fMax[i] )
			in=0;
	}

	// calculate the log probabilities and compare old and new point
	double p0 = this -> LogEvalSampling(fXmetro0); // old point
	double p1=0; // new point
	if(in)
		p1= this -> LogEvalSampling(fXmetro1);

	// compare log probabilities
	int accept=0;
	if(in)
	{
		if(p1>=p0)
			accept=1;
		else
		{
			double r=log(fRandom->Rndm());
			if(r<p1-p0)
				accept=1;
		}
	}

	// fill the return point after the decision
	if(accept)
		for(int i=0;i<fNvar;i++)
		{
			fXmetro0[i]=fXmetro1[i];
			x[i]=fXmetro1[i];
		}
	else
		for(int i=0;i<fNvar;i++)
			x[i]=fXmetro0[i];

	fNmetro++;
}

// *********************************************
void BCIntegrate::GetRandomVector(std::vector <double> &x)
{
	double * randx = new double[fNvar];

	fRandom -> RndmArray(fNvar, randx);

	for(int i=0;i<fNvar;i++)
		x[i] = randx[i];

	delete[] randx;
	randx = 0;
}

// *********************************************
void BCIntegrate::GetRandomVectorMetro(std::vector <double> &x)
{
	double * randx = new double[fNvar];

	fRandom -> RndmArray(fNvar, randx);

	for(int i=0;i<fNvar;i++)
		x[i] = 2.0 * (randx[i] - 0.5) * fMarkovChainStepSize;

	delete[] randx;
	randx = 0;
}

// *********************************************
TMinuit * BCIntegrate::GetMinuit()
{
	if (!fMinuit)
		fMinuit = new TMinuit();

	return fMinuit;
}

// *********************************************

void BCIntegrate::FindModeMinuit(std::vector<double> start, int printlevel)
{
	bool have_start = true;

	// check start values
	if (int(start.size()) != fNvar)
		have_start = false;

	// set global this
	global_this = this;

	// define minuit
	if (fMinuit)
		delete fMinuit;
	fMinuit = new TMinuit(fNvar);

	// set function
	fMinuit -> SetFCN(&BCIntegrate::FCNLikelihood);

	// set print level
	fMinuit -> SetPrintLevel(printlevel);

	// set parameters
	int flag;
	for (int i = 0; i < fNvar; i++)
	{
		double starting_point = (fMin[i]+fMax[i])/2.;
		if(have_start)
			starting_point = start[i];
		fMinuit -> mnparm(i,
				fx -> at(i) -> GetName().data(),
				starting_point,
				(fMax[i]-fMin[i])/100.0,
				fMin[i],
				fMax[i],
				flag);
	}

	// do minimization
	fMinuit -> mnexcm("MIGRAD", fMinuitArglist, 2, flag);

	// copy flag
	fMinuitErrorFlag = flag;

	// set best fit parameters
	fBestFitParameters.clear();

	for (int i = 0; i < fNvar; i++)
	{
		double par;
		double parerr;
		fMinuit -> GetParameter(i, par, parerr);
		fBestFitParameters.push_back(par);
	}

	// delete minuit
		delete fMinuit;
		fMinuit = 0;

	return;
}

// *********************************************

void BCIntegrate::FindModeSA(std::vector<double> start)
{
	// note: if f(x) is the function to be minimized, then
	// f(x) := - this->LogEval(parameters)

	bool have_start = true;
	std::vector<double> x, y, best_fit; // vectors for current state, new proposed state and best fit up to now
	double fval_x, fval_y, fval_best_fit; // function values at points x, y and best_fit (we save them rather than to re-calculate them every time)
	int t = 1; // time iterator

	// check start values
	if (int(start.size()) != fNvar)
		have_start = false;

	// if no starting point is given, set to center of parameter space
	if ( !have_start )
	{
		start.clear();
		for (int i = 0; i < fNvar; i++)
			start.push_back((fMin[i]+fMax[i])/2.);
	}

	// set current state and best fit to starting point
	x.clear();
	best_fit.clear();
	for (int i = 0; i < fNvar; i++)
	{
		x.push_back(start[i]);
		best_fit.push_back(start[i]);
	}
	// calculate function value at starting point
	fval_x = fval_best_fit = this->LogEval(x);

	// run while still "hot enough"
	while ( this->SATemperature(t) > fSATmin )
	{

		// generate new state
		y = this->GetProposalPointSA(x, t);

		// check if the proposed point is inside the phase space
		// if not, reject it
		bool is_in_ranges = true;

		for (int i = 0; i < fNvar; i++)
			if (y[i] > fMax[i] || y[i] < fMin[i])
				is_in_ranges = false;
		
		if ( !is_in_ranges )
			; // do nothing...
		else
		{
			// calculate function value at new point
			fval_y = this->LogEval(y);

			// is it better than the last one?
			// if so, update state and chef if it is the new best fit...
			if (fval_y >= fval_x)
			{
				x.clear();
				for (int i = 0; i < fNvar; i++)
					x.push_back(y[i]);

				fval_x = fval_y;

				if (fval_y > fval_best_fit)
				{
					best_fit.clear();
					for (int i = 0; i < fNvar; i++)
						best_fit.push_back(y[i]);
					
					fval_best_fit = fval_y;
				}
			}
			// ...else, only accept new state w/ certain probability
			else
			{
				if (fRandom->Rndm() <= exp( (fval_y - fval_x)
						/ this->SATemperature(t) ))
				{
					x.clear();
					for (int i = 0; i < fNvar; i++)
						x.push_back(y[i]);
					
					fval_x = fval_y;
				}
			}
		}
		t++;
	}

	// set best fit parameters
	fBestFitParameters.clear();

	for (int i = 0; i < fNvar; i++)
		fBestFitParameters.push_back(best_fit[i]);

	return;
}

// *********************************************

double BCIntegrate::SATemperature(double t)
{
	// do we have Cauchy (default) or Boltzmann annealing schedule?
	if (this->fSASchedule == BCIntegrate::kSABoltzmann)
		return this->SATemperatureBoltzmann(t);
	else
		return this->SATemperatureCauchy(t);
}

// *********************************************

double BCIntegrate::SATemperatureBoltzmann(double t)
{
	return fSAT0 / log((double)(t + 1));
}

// *********************************************

double BCIntegrate::SATemperatureCauchy(double t)
{
	return fSAT0 / (double)t;
}

// *********************************************

std::vector<double> BCIntegrate::GetProposalPointSA(std::vector<double> x, int t)
{
	// do we have Cauchy (default) or Boltzmann annealing schedule?
	if (this->fSASchedule == BCIntegrate::kSABoltzmann)
		return this->GetProposalPointSABoltzmann(x, t);
	else
		return this->GetProposalPointSACauchy(x, t);
}

// *********************************************

std::vector<double> BCIntegrate::GetProposalPointSABoltzmann(std::vector<double> x, int t)
{
	std::vector<double> y;
	y.clear();
	double new_val, norm;

	for (int i = 0; i < fNvar; i++)
	{
		norm = (fMax[i] - fMin[i]) * this->SATemperature(t) / 2.;
		new_val = x[i] + norm * fRandom->Gaus();
		y.push_back(new_val);
	}

	return y;
}

// *********************************************

std::vector<double> BCIntegrate::GetProposalPointSACauchy(std::vector<double> x, int t)
{
	std::vector<double> y;
	y.clear();

	if (fNvar == 1)
	{
		double cauchy, new_val, norm;

		norm = (fMax[0] - fMin[0]) * this->SATemperature(t) / 2.;
		cauchy = tan(3.14159 * (fRandom->Rndm() - 0.5));
		new_val = x[0] + norm * cauchy;
		y.push_back(new_val);
	}
	else
	{
		// use sampling to get radial n-dim Cauchy distribution

		// first generate a random point uniformly distributed on a
		// fNvar-dimensional hypersphere
		y = this->SAHelperGetRandomPointOnHypersphere();

		// scale the vector by a random factor determined by the radial
		// part of the fNvar-dimensional Cauchy distribution
		double radial = this->SATemperature(t) * this->SAHelperGetRadialCauchy();

		// scale y by radial part and the size of dimension i in phase space
		// afterwards, move by x
		for (int i = 0; i < fNvar; i++)
			y[i] = (fMax[i] - fMin[i]) * y[i] * radial / 2. + x[i];
	}

	return y;
}

// *********************************************

std::vector<double> BCIntegrate::SAHelperGetRandomPointOnHypersphere()
{
	std::vector<double> rand_point, gauss_array;
	double s = 0.,
		gauss_num;

	for (int i = 0; i < fNvar; i++)
	{
		gauss_num = fRandom->Gaus();
		gauss_array.push_back(gauss_num);
		s += gauss_num * gauss_num;
	}
	s = sqrt(s);

	for (int i = 0; i < fNvar; i++)
		rand_point.push_back(gauss_array[i] / s);

	return rand_point;
}

// *********************************************

double BCIntegrate::SAHelperGetRadialCauchy()
{
	// theta is sampled from a rather complicated distribution,
	// so first we create a lookup table with 10000 random numbers
	// once and then, each time we need a new random number,
	// we just look it up in the table.
	double theta;

	// static vectors for theta-sampling-map
	static std::vector<double> map_u (10001);
	static std::vector<double> map_theta (10001);
	static bool initialized = false;
	static int map_dimension = 0;

	// is the lookup-table already initialized? if not, do it!
	if (!initialized || map_dimension != fNvar)
	{
		double init_theta;
		double init_cdf;
		double beta = this->SAHelperSinusToNIntegral(fNvar - 1, 1.57079632679);

		for (int i = 0; i <= 10000; i++)
		{
			init_theta = 3.14159265 * (double)i / 5000.;
			map_theta.push_back(init_theta);

			init_cdf = this->SAHelperSinusToNIntegral(fNvar - 1, init_theta) / beta;
			map_u.push_back(init_cdf);
		}

		map_dimension = fNvar;
		initialized = true;
	} // initializing is done.

	// generate uniform random number for sampling
	double u = fRandom->Uniform();

	// Find the two elements just greater than and less than u
	// using a binary search (O(log(N))).
	int lo = 0;
	int up = map_u.size() - 1;
	int mid;

	while (up != lo)
	{
		mid = ((up - lo + 1) / 2) + lo;

		if (u >= map_u[mid])
			lo = mid;
		else
			up = mid - 1;
	}
	up++;

	// perform linear interpolation:
	theta = map_theta[lo] + (u - map_u[lo]) / (map_u[up] - map_u[lo])
		* (map_theta[up] - map_theta[lo]);

	return tan(theta);
}

// *********************************************

double BCIntegrate::SAHelperSinusToNIntegral(int dim, double theta)
{
	if (dim < 1)
		return theta;
	else if (dim == 1)
		return (1. - cos(theta));
	else if (dim == 2)
		return 0.5 * (theta - sin(theta) * cos(theta));
	else if (dim == 3)
		return (2. - sin(theta) * sin(theta) * cos(theta) - 2. * cos(theta)) / 3.;
	else
		return - pow(sin(theta), (double)(dim - 1)) * cos(theta) / (double)dim
			+ (double)(dim - 1) / (double)dim
			* this->SAHelperSinusToNIntegral(dim - 2, theta);
}

// *********************************************

void BCIntegrate::SetMode(std::vector <double> mode)
{
	if((int)mode.size() == fNvar)
	{
		fBestFitParameters.clear();
		for (int i = 0; i < fNvar; i++)
			fBestFitParameters.push_back(mode[i]);
	}
}

// *********************************************

void BCIntegrate::FCNLikelihood(int &npar, double * grad, double &fval, double * par, int flag)
{
	// copy parameters
	std::vector <double> parameters;

	int n = global_this -> GetNvar();

	for (int i = 0; i < n; i++)
		parameters.push_back(par[i]);

	fval = - global_this -> LogEval(parameters);

	// delete parameters
	parameters.clear();
}

// *********************************************

//void BCIntegrate::FindModeMCMC(int flag_run)
void BCIntegrate::FindModeMCMC()
{
	// call PreRun
//	if (flag_run == 0)
//	if (!fMCMCFlagPreRun)
//		this -> MCMCMetropolisPreRun();

	// find global maximum
	double probmax = (this -> MCMCGetMaximumLogProb()).at(0);
	fBestFitParameters = this -> MCMCGetMaximumPoint(0);

	// loop over all chains and find the maximum point
	for (int i = 1; i < fMCMCNChains; ++i)
	{
		double prob = (this -> MCMCGetMaximumLogProb()).at(i);

		// copy the point into the vector
		if (prob > probmax)
		{
			probmax = prob;

			fBestFitParameters.clear();

			fBestFitParameters = this -> MCMCGetMaximumPoint(i);
		}
	}
}

// *********************************************
void BCIntegrate::CubaIntegrand(const int *ndim, const double xx[],
		const int *ncomp, double ff[])
{
#ifdef HAVE_CUBA_H
	// scale variables
	double jacobian = 1.0;

	std::vector<double> scaled_parameters;

	for (int i = 0; i < *ndim; i++)
	{
		double range = global_this -> fx -> at(i) -> GetUpperLimit() -  global_this -> fx -> at(i) -> GetLowerLimit();

		// multiply range to jacobian
		jacobian *= range;

		// get the scaled parameter value
		scaled_parameters.push_back(global_this -> fx -> at(i) -> GetLowerLimit() + xx[i] * range);
	}

	// call function to integrate
	ff[0] = global_this -> Eval(scaled_parameters);

	// multiply jacobian
	ff[0] *= jacobian;

	// multiply fudge factor
	ff[0] *= 1e99;

	// remove parameter vector
	scaled_parameters.clear();
#else
	BCLog::Out(BCLog::error,BCLog::error,"!!! This version of BAT is compiled without Cuba.");
	BCLog::Out(BCLog::error,BCLog::error,"    Use other integration methods or install Cuba and recompile BAT.");
#endif
}

// *********************************************
double BCIntegrate::CubaIntegrate()
{
#ifdef HAVE_CUBA_H
	double EPSREL = 1e-3;
	double EPSABS = 1e-12;
	double VERBOSE   = 0;
	double MINEVAL   = 0;
	double MAXEVAL   = 2000000;
	double NSTART    = 25000;
	double NINCREASE = 25000;

	std::vector<double> parameters_double;
	std::vector<double>    parameters_int;

	parameters_double.push_back(EPSREL);
	parameters_double.push_back(EPSABS);

	parameters_int.push_back(VERBOSE);
	parameters_int.push_back(MINEVAL);
	parameters_int.push_back(MAXEVAL);
	parameters_int.push_back(NSTART);
	parameters_int.push_back(NINCREASE);

	// print to log
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Running Cuba/Vegas integation over %i dimensions.", fNvar));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Maximum number of iterations: %i", (int)MAXEVAL));
	BCLog::Out(BCLog::detail, BCLog::detail, Form(" --> Aimed relative precision:     %e", EPSREL));

	return this -> CubaIntegrate(0, parameters_double, parameters_int);
#else
	BCLog::Out(BCLog::error,BCLog::error,"!!! This version of BAT is compiled without Cuba.");
	BCLog::Out(BCLog::error,BCLog::error,"    Use other integration methods or install Cuba and recompile BAT.");
	return 0.;
#endif
}

// *********************************************
double BCIntegrate::CubaIntegrate(int method, std::vector<double> parameters_double, std::vector<double> parameters_int)
{
#ifdef HAVE_CUBA_H
	const int NDIM      = int(fx ->size());
	const int NCOMP     = 1;

	const double EPSREL = parameters_double[0];
	const double EPSABS = parameters_double[1];
	const int VERBOSE   = int(parameters_int[0]);
	const int MINEVAL   = int(parameters_int[1]);
	const int MAXEVAL   = int(parameters_int[2]);

	int neval;
	int fail;
	int nregions;
	double integral[NCOMP];
	double error[NCOMP];
	double prob[NCOMP];

	global_this = this;

	integrand_t an_integrand = &BCIntegrate::CubaIntegrand;

	if (method == 0)
	{
		// set VEGAS specific parameters
		const int NSTART    = int(parameters_int[3]);
		const int NINCREASE = int(parameters_int[4]);

		// call VEGAS integration method
		Vegas(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
			NSTART, NINCREASE,
			&neval, &fail, integral, error, prob);
	}
	else if (method == 1)
	{
		// set SUAVE specific parameters
		const int LAST     = int(parameters_int[3]);
		const int NNEW     = int(parameters_int[4]);
		const int FLATNESS = int(parameters_int[5]);

		// call SUAVE integration method
		Suave(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
			NNEW, FLATNESS,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else if (method == 2)
	{
		// set DIVONNE specific parameters
		const int KEY1         = int(parameters_int[3]);
		const int KEY2         = int(parameters_int[4]);
		const int KEY3         = int(parameters_int[5]);
		const int MAXPASS      = int(parameters_int[6]);
		const int BORDER       = int(parameters_int[7]);
		const int MAXCHISQ     = int(parameters_int[8]);
		const int MINDEVIATION = int(parameters_int[9]);
		const int NGIVEN       = int(parameters_int[10]);
		const int LDXGIVEN     = int(parameters_int[11]);
		const int NEXTRA       = int(parameters_int[12]);

		// call DIVONNE integration method
		Divonne(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
			KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
			NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else if (method == 3)
	{
		// set CUHRE specific parameters
		const int LAST = int(parameters_int[3]);
		const int KEY  = int(parameters_int[4]);

		// call CUHRE integration method
		Cuhre(NDIM, NCOMP, an_integrand,
			EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL, KEY,
			&nregions, &neval, &fail, integral, error, prob);
	}
	else
	{
		std::cout << " Integration method not available. " << std::endl;
		integral[0] = -1e99;
	}

	if (fail != 0)
	{
		std::cout << " Warning, integral did not converge with the given set of parameters. "<< std::endl;
		std::cout << " neval    = " << neval       << std::endl;
		std::cout << " fail     = " << fail        << std::endl;
		std::cout << " integral = " << integral[0] << std::endl;
		std::cout << " error    = " << error[0]    << std::endl;
		std::cout << " prob     = " << prob[0]     << std::endl;
	}

	return integral[0] / 1e99;
#else
	BCLog::Out(BCLog::error,BCLog::error,"!!! This version of BAT is compiled without Cuba.");
	BCLog::Out(BCLog::error,BCLog::error,"    Use other integration methods or install Cuba and recompile BAT.");
	return 0.;
#endif
}

// *********************************************
void BCIntegrate::MCMCIterationInterface()
{
	// what's within this method will be executed
	// for every iteration of the MCMC

	// fill error band
	this -> MCMCFillErrorBand();

	// do user defined stuff
	this -> MCMCUserIterationInterface();
}

// *********************************************
void BCIntegrate::MCMCFillErrorBand()
{
	if (!fFillErrorBand)
		return;

	// function fitting
	if (fFitFunctionIndexX < 0)
		return;

	// loop over all possible x values ...
	if (fErrorBandContinuous)
	{
		double x = 0;
		for (int ix = 0; ix < fErrorBandNbinsX; ix++)
		{
			// calculate x
			x = fErrorBandXY -> GetXaxis() -> GetBinCenter(ix + 1);

			// calculate y
			std::vector <double> xvec;
			xvec.push_back(x);

			// loop over all chains
			for (int ichain = 0; ichain < this -> MCMCGetNChains(); ++ichain)
			{
				// calculate y
				double y = this -> FitFunction(xvec, this -> MCMCGetx(ichain));

				// fill histogram
				fErrorBandXY -> Fill(x, y);
			}

			xvec.clear();
		}
	}
	// ... or evaluate at the data point x-values
	else
	{
		int ndatapoints = int(fErrorBandX.size());
		double x = 0;

		for (int ix = 0; ix < ndatapoints; ++ix)
		{
			// calculate x
			x = fErrorBandX.at(ix);

			// calculate y
			std::vector <double> xvec;
			xvec.push_back(x);

			// loop over all chains
			for (int ichain = 0; ichain < this -> MCMCGetNChains(); ++ichain)
			{
				// calculate y
				double y = this -> FitFunction(xvec, this -> MCMCGetx(ichain));

				// fill histogram
				fErrorBandXY -> Fill(x, y);
			}

			xvec.clear();
		}
	}
}

// *********************************************

