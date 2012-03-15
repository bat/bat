/*
 * Copyright (C) 2008-2011, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "config.h"

#include "BCIntegrate.h"
#include "BCLog.h"
#include "BCMath.h"

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

// ---------------------------------------------------------
BCIntegrate::BCIntegrate()
   : BCEngineMCMC()
   , fNvar(0)
   , fNbins(100)
   , fNSamplesPer2DBin(100)
   , fMarkovChainStepSize(0.1)
   , fMarkovChainAutoN(true)
   , fDataPointLowerBoundaries(0)
   , fDataPointUpperBoundaries(0)
   , fFillErrorBand(false)
   , fFitFunctionIndexX(-1)
   , fFitFunctionIndexY(-1)
   , fErrorBandContinuous(true)
   , fErrorBandNbinsX(100)
   , fErrorBandNbinsY(500)
   , fMinuit(0)
   , fFlagIgnorePrevOptimization(false)
   , fFlagWriteMarkovChain(false)
   , fMarkovChainTree(0)
   , fSAT0(100)
   , fSATmin(0.1)
   , fTreeSA(0)
   , fFlagWriteSAToFile(false)
	 , fIntegrationVolume(0)
#ifdef HAVE_CUBA_H
   , fIntegrationMethod(BCIntegrate::kIntCuba)
#else
   , fIntegrationMethod(BCIntegrate::kIntMonteCarlo)
#endif
   , fMarginalizationMethod(BCIntegrate::kMargMetropolis)
   , fOptimizationMethod(BCIntegrate::kOptMinuit)
   , fOptimizationMethodMode(BCIntegrate::kOptMinuit)
   , fSASchedule(BCIntegrate::kSACauchy)
   , fNIterationsMin(0)
   , fNIterationsMax(1000000)
   , fNIterationsPrecisionCheck(1000)
   , fNIterationsOutput(-1)
   , fNIterations(0)
   , fRelativePrecision(1e-3)
   , fAbsolutePrecision(1e-12)
   , fCubaIntegrationMethod(BCIntegrate::kCubaVegas)
   , fCubaVerbosity(0)
   , fCubaVegasNStart(25000)
   , fCubaVegasNIncrease(25000)
   , fCubaSuaveNNew(10000)
   , fCubaSuaveFlatness(50)
   , fError(-999.)
{
   fMinuitArglist[0] = 20000;
   fMinuitArglist[1] = 0.01;
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(BCParameterSet * par)
   : BCEngineMCMC()
   , fNvar(0)
   , fNbins(100)
   , fNSamplesPer2DBin(100)
   , fMarkovChainStepSize(0.1)
   , fMarkovChainAutoN(true)
   , fDataPointLowerBoundaries(0)
   , fDataPointUpperBoundaries(0)
   , fFillErrorBand(false)
   , fFitFunctionIndexX(-1)
   , fFitFunctionIndexY(-1)
   , fErrorBandContinuous(true)
   , fErrorBandNbinsX(100)
   , fErrorBandNbinsY(500)
   , fMinuit(0)
   , fFlagIgnorePrevOptimization(false)
   , fFlagWriteMarkovChain(false)
   , fMarkovChainTree(0)
   , fSAT0(100)
   , fSATmin(0.1)
   , fTreeSA(0)
   , fFlagWriteSAToFile(false)
	 , fIntegrationVolume(0)
#ifdef HAVE_CUBA_H
   , fIntegrationMethod(BCIntegrate::kIntCuba)
#else
   , fIntegrationMethod(BCIntegrate::kIntMonteCarlo)
#endif
   , fMarginalizationMethod(BCIntegrate::kMargMetropolis)
   , fOptimizationMethod(BCIntegrate::kOptMinuit)
   , fOptimizationMethodMode(BCIntegrate::kOptMinuit)
   , fSASchedule(BCIntegrate::kSACauchy)
   , fNIterationsMin(0)
   , fNIterationsMax(1000000)
	 , fNIterationsPrecisionCheck(1000)
	 , fNIterationsOutput(-1)
   , fNIterations(0)
   , fRelativePrecision(1e-3)
   , fAbsolutePrecision(1e-12)
   , fCubaIntegrationMethod(BCIntegrate::kCubaVegas)
   , fCubaVerbosity(0)
   , fCubaVegasNStart(25000)
   , fCubaVegasNIncrease(25000)
   , fCubaSuaveNNew(10000)
   , fCubaSuaveFlatness(50)
   , fError(-999.)
{
   SetParameters(par);

   fMinuitArglist[0] = 20000;
   fMinuitArglist[1] = 0.01;
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(const BCIntegrate & bcintegrate) : BCEngineMCMC(bcintegrate)
 {
    fNvar                     = bcintegrate.fNvar;
    fNbins                    = bcintegrate.fNbins;
    fNSamplesPer2DBin         = bcintegrate.fNSamplesPer2DBin;
    fMarkovChainStepSize      = bcintegrate.fMarkovChainStepSize;
    fMarkovChainNIterations   = bcintegrate.fMarkovChainNIterations;
    fMarkovChainAutoN         = bcintegrate.fMarkovChainAutoN;
    if (bcintegrate.fDataPointLowerBoundaries)
       fDataPointLowerBoundaries = new BCDataPoint(*bcintegrate.fDataPointLowerBoundaries);
    else
       fDataPointLowerBoundaries = 0;
    if (bcintegrate.fDataPointUpperBoundaries)
       fDataPointUpperBoundaries = new BCDataPoint(*bcintegrate.fDataPointUpperBoundaries);
    else
       fDataPointUpperBoundaries = 0;
    fDataFixedValues          = bcintegrate.fDataFixedValues;
    fBestFitParameters        = bcintegrate.fBestFitParameters;
    fBestFitParameterErrors   = bcintegrate.fBestFitParameterErrors;
    fBestFitParametersMarginalized = bcintegrate.fBestFitParametersMarginalized;
    for (int i = 0; i < int(bcintegrate.fHProb1D.size()); ++i) {
       if (bcintegrate.fHProb1D.at(i))
          fHProb1D.push_back(new TH1D(*(bcintegrate.fHProb1D.at(i))));
       else
          fHProb1D.push_back(0);
    }
    for (int i = 0; i < int(bcintegrate.fHProb2D.size()); ++i) {
       if (bcintegrate.fHProb2D.at(i))
          fHProb2D.push_back(new TH2D(*(fHProb2D.at(i))));
       else
          fHProb2D.push_back(0);
    }
    fFillErrorBand            = bcintegrate.fFillErrorBand;
    fFitFunctionIndexX        = bcintegrate.fFitFunctionIndexX;
    fFitFunctionIndexY        = bcintegrate.fFitFunctionIndexY;
    fErrorBandX               = bcintegrate.fErrorBandX;
    if (bcintegrate.fErrorBandXY)
       fErrorBandXY = new TH2D(*(bcintegrate.fErrorBandXY));
    else
       fErrorBandXY = 0;
    fErrorBandNbinsX          = bcintegrate.fErrorBandNbinsX;
    fErrorBandNbinsY          = bcintegrate.fErrorBandNbinsY;
    fMinuit                   = new TMinuit();
    // debugKK
    //    *fMinuit = *(bcintegrate.fMinuit);
    fMinuitArglist[0]         = bcintegrate.fMinuitArglist[0];
    fMinuitArglist[1]         = bcintegrate.fMinuitArglist[1];
    fMinuitErrorFlag          = bcintegrate.fMinuitErrorFlag;
    fFlagIgnorePrevOptimization = bcintegrate.fFlagIgnorePrevOptimization;
    fFlagWriteMarkovChain     = bcintegrate.fFlagWriteMarkovChain;
    fMarkovChainTree          = bcintegrate.fMarkovChainTree;
    fMCMCIteration            = bcintegrate.fMCMCIteration;
    fSAT0                     = bcintegrate.fSAT0;
    fSATmin                   = bcintegrate.fSATmin;
    // debugKK
    fTreeSA = 0;
    fFlagWriteSAToFile        = bcintegrate.fFlagWriteSAToFile;
    fSANIterations            = bcintegrate.fSANIterations;
    fSATemperature            = bcintegrate.fSATemperature;
    fSALogProb                = bcintegrate.fSALogProb;
    fSAx                      = bcintegrate.fSAx;
    if (bcintegrate.fx)
       fx = new BCParameterSet(*(bcintegrate.fx));
    else
       fx = 0;
    fMin                      = new double[fNvar];
    fMax                      = new double[fNvar];
    fVarlist                  = new int[fNvar];
    fMin                      = bcintegrate.fMin;
    fMax                      = bcintegrate.fMax;
    fVarlist                  = bcintegrate.fVarlist;
		fIntegrationVolume        = bcintegrate.fIntegrationVolume;
    fIntegrationMethod        = bcintegrate.fIntegrationMethod;
    fMarginalizationMethod    = bcintegrate.fMarginalizationMethod;
    fOptimizationMethod       = bcintegrate.fOptimizationMethod;
    fOptimizationMethodMode   = bcintegrate.fOptimizationMethodMode;
    fSASchedule               = bcintegrate.fSASchedule;
    fNIterationsMin           = bcintegrate.fNIterationsMin;
    fNIterationsMax           = bcintegrate.fNIterationsMax;
		fNIterationsPrecisionCheck= bcintegrate.fNIterationsPrecisionCheck;
		fNIterationsOutput        = bcintegrate.fNIterationsOutput;
    fNIterations              = bcintegrate.fNIterations;
    fRelativePrecision        = bcintegrate.fRelativePrecision;
    fAbsolutePrecision        = bcintegrate.fAbsolutePrecision;
    fCubaIntegrationMethod    = bcintegrate.fCubaIntegrationMethod;
    fCubaVerbosity            = bcintegrate.fCubaVerbosity;
    fCubaVegasNStart          = bcintegrate.fCubaVegasNStart;
    fCubaVegasNIncrease       = bcintegrate.fCubaVegasNIncrease;
    fCubaSuaveNNew            = bcintegrate.fCubaSuaveNNew;
    fCubaSuaveFlatness        = bcintegrate.fCubaSuaveFlatness;
    fError                    = bcintegrate.fError;
    fNmetro                   = bcintegrate.fNmetro;
    fNacceptedMCMC            = bcintegrate.fNacceptedMCMC;
    fXmetro0                  = bcintegrate.fXmetro0;
    fXmetro1                  = bcintegrate.fXmetro1;
    fMarkovChainValue         = bcintegrate.fMarkovChainValue;
}

// ---------------------------------------------------------
BCIntegrate & BCIntegrate::operator = (const BCIntegrate & bcintegrate)
{
    BCEngineMCMC::operator=(bcintegrate);

    fNvar                     = bcintegrate.fNvar;
    fNbins                    = bcintegrate.fNbins;
    fNSamplesPer2DBin         = bcintegrate.fNSamplesPer2DBin;
    fMarkovChainStepSize      = bcintegrate.fMarkovChainStepSize;
    fMarkovChainNIterations   = bcintegrate.fMarkovChainNIterations;
    fMarkovChainAutoN         = bcintegrate.fMarkovChainAutoN;
    if (bcintegrate.fDataPointLowerBoundaries)
       fDataPointLowerBoundaries = new BCDataPoint(*bcintegrate.fDataPointLowerBoundaries);
    else
       fDataPointLowerBoundaries = 0;
    if (bcintegrate.fDataPointUpperBoundaries)
       fDataPointUpperBoundaries = new BCDataPoint(*bcintegrate.fDataPointUpperBoundaries);
    else
       fDataPointUpperBoundaries = 0;
    fDataFixedValues          = bcintegrate.fDataFixedValues;
    fBestFitParameters        = bcintegrate.fBestFitParameters;
    fBestFitParameterErrors   = bcintegrate.fBestFitParameterErrors;
    fBestFitParametersMarginalized = bcintegrate.fBestFitParametersMarginalized;
    for (int i = 0; i < int(bcintegrate.fHProb1D.size()); ++i) {
       if (bcintegrate.fHProb1D.at(i))
          fHProb1D.push_back(new TH1D(*(bcintegrate.fHProb1D.at(i))));
       else
          fHProb1D.push_back(0);
    }
    for (int i = 0; i < int(bcintegrate.fHProb2D.size()); ++i) {
       if (bcintegrate.fHProb2D.at(i))
          fHProb2D.push_back(new TH2D(*(fHProb2D.at(i))));
       else
          fHProb2D.push_back(0);
    }
    fFillErrorBand            = bcintegrate.fFillErrorBand;
    fFitFunctionIndexX        = bcintegrate.fFitFunctionIndexX;
    fFitFunctionIndexY        = bcintegrate.fFitFunctionIndexY;
    fErrorBandX               = bcintegrate.fErrorBandX;
    if (bcintegrate.fErrorBandXY)
       fErrorBandXY = new TH2D(*(bcintegrate.fErrorBandXY));
    else
       fErrorBandXY = 0;
    fErrorBandNbinsX          = bcintegrate.fErrorBandNbinsX;
    fErrorBandNbinsY          = bcintegrate.fErrorBandNbinsY;
    fMinuit                   = new TMinuit();
    // debugKK
    //    *fMinuit = *(bcintegrate.fMinuit);
    fMinuitArglist[0]         = bcintegrate.fMinuitArglist[0];
    fMinuitArglist[1]         = bcintegrate.fMinuitArglist[1];
    fMinuitErrorFlag          = bcintegrate.fMinuitErrorFlag;
    fFlagIgnorePrevOptimization = bcintegrate.fFlagIgnorePrevOptimization;
    fFlagWriteMarkovChain     = bcintegrate.fFlagWriteMarkovChain;
    fMarkovChainTree          = bcintegrate.fMarkovChainTree;
    fMCMCIteration            = bcintegrate.fMCMCIteration;
    fSAT0                     = bcintegrate.fSAT0;
    fSATmin                   = bcintegrate.fSATmin;
    // debugKK
    fTreeSA = 0;
    fFlagWriteSAToFile        = bcintegrate.fFlagWriteSAToFile;
    fSANIterations            = bcintegrate.fSANIterations;
    fSATemperature            = bcintegrate.fSATemperature;
    fSALogProb                = bcintegrate.fSALogProb;
    fSAx                      = bcintegrate.fSAx;
    if (bcintegrate.fx)
       fx = new BCParameterSet(*(bcintegrate.fx));
    else
       fx = 0;
    fMin                      = new double[fNvar];
    fMax                      = new double[fNvar];
    fVarlist                  = new int[fNvar];
		fIntegrationVolume        = bcintegrate.fIntegrationVolume;
    fIntegrationMethod        = bcintegrate.fIntegrationMethod;
    fMarginalizationMethod    = bcintegrate.fMarginalizationMethod;
    fOptimizationMethod       = bcintegrate.fOptimizationMethod;
    fOptimizationMethodMode   = bcintegrate.fOptimizationMethodMode;
    fSASchedule               = bcintegrate.fSASchedule;

    fNIterationsMin           = bcintegrate.fNIterationsMin;
    fNIterationsMax           = bcintegrate.fNIterationsMax;
		fNIterationsPrecisionCheck= bcintegrate.fNIterationsPrecisionCheck;
		fNIterationsOutput        = bcintegrate.fNIterationsOutput;
    fNIterations              = bcintegrate.fNIterations;
    fRelativePrecision        = bcintegrate.fRelativePrecision;
    fAbsolutePrecision        = bcintegrate.fAbsolutePrecision;
    fCubaIntegrationMethod    = bcintegrate.fCubaIntegrationMethod;
    fCubaVerbosity            = bcintegrate.fCubaVerbosity;
    fCubaVegasNStart          = bcintegrate.fCubaVegasNStart;
    fCubaVegasNIncrease       = bcintegrate.fCubaVegasNIncrease;
    fCubaSuaveNNew            = bcintegrate.fCubaSuaveNNew;
    fCubaSuaveFlatness        = bcintegrate.fCubaSuaveFlatness;
    fError                    = bcintegrate.fError;
    fNmetro                   = bcintegrate.fNmetro;
    fNacceptedMCMC            = bcintegrate.fNacceptedMCMC;
    fXmetro0                  = bcintegrate.fXmetro0;
    fXmetro1                  = bcintegrate.fXmetro1;
    fMarkovChainValue         = bcintegrate.fMarkovChainValue;

   // return this
   return *this;
}

// ---------------------------------------------------------
BCIntegrate::~BCIntegrate()
{
   DeleteVarList();

   fx=0;

   if (fMinuit)
      delete fMinuit;

   int n1 = fHProb1D.size();
   if(n1>0) {
      for (int i=0;i<n1;i++)
         delete fHProb1D.at(i);
   }

   int n2 = fHProb2D.size();
   if(n2>0) {
      for (int i=0;i<n2;i++)
         delete fHProb2D.at(i);
   }
}

// ---------------------------------------------------------
void BCIntegrate::SetParameters(BCParameterSet * par)
{
   DeleteVarList();

   fx = par;
   fNvar = fx->size();
   fMin = new double[fNvar];
   fMax = new double[fNvar];
   fVarlist = new int[fNvar];

   ResetVarlist(1);

   for(int i=0;i<fNvar;i++) {
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

   for(int i=0;i<fNvar;i++) {
      fMCMCBoundaryMin.push_back(fMin[i]);
      fMCMCBoundaryMax.push_back(fMax[i]);
      fMCMCFlagsFillHistograms.push_back(true);
   }

   for (int i = int(fMCMCH1NBins.size()); i<fNvar; ++i)
      fMCMCH1NBins.push_back(100);

   fMCMCNParameters = fNvar;
}

// ---------------------------------------------------------
void BCIntegrate::SetMarkovChainInitialPosition(std::vector<double> position)
{
   int n = position.size();

   fXmetro0.clear();

   for (int i = 0; i < n; ++i)
      fXmetro0.push_back(position.at(i));
}

// ---------------------------------------------------------
void BCIntegrate::DeleteVarList()
{
   if(fNvar) {
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

// ---------------------------------------------------------
void BCIntegrate::SetNbins(int nbins, int index)
{
   if (fNvar == 0)
      return;

   // check if index is in range
   if (index >= fNvar) {
      BCLog::OutWarning("BCIntegrate::SetNbins : Index out of range.");
      return;
   }
   // set for all parameters at once
   else if (index < 0) {
      for (int i = 0; i < fNvar; ++i)
         SetNbins(nbins, i);
      return;
   }

   // sanity check for nbins
   if (nbins <= 0)
      nbins = 10;

   fMCMCH1NBins[index] = nbins;

   return;

//    if(n<2) {
//       BCLog::OutWarning("BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
//       n=2;
//    }
//    MCMCSetH1NBins(n, -1);

   //   fNbins=n;

   //   fMCMCH1NBins = n;
   //   fMCMCH2NBinsX = n;
   //   fMCMCH2NBinsY = n;
}

// ---------------------------------------------------------
// void BCIntegrate::SetNbinsX(int n)
// {
//    if(n<2) {
//       BCLog::OutWarning("BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
//       n=2;
//    }
//    fMCMCH2NBinsX = n;
// }

// ---------------------------------------------------------
// void BCIntegrate::SetNbinsY(int n)
// {
//    if(n<2) {
//       BCLog::OutWarning("BCIntegrate::SetNbins. Number of bins less than 2 makes no sense. Setting to 2.");
//       n=2;
//    }
//    fNbins=n;

//    fMCMCH2NBinsY = n;
// }

// ---------------------------------------------------------
void BCIntegrate::SetVarList(int * varlist)
{
   for(int i=0;i<fNvar;i++)
      fVarlist[i]=varlist[i];
}

// ---------------------------------------------------------
void BCIntegrate::ResetVarlist(int v)
{
   for(int i=0;i<fNvar;i++)
      fVarlist[i]=v;
}

// ---------------------------------------------------------
double BCIntegrate::Eval(const std::vector<double> & /*x*/)
{
   BCLog::OutWarning( "BCIntegrate::Eval. No function. Function needs to be overloaded.");
   return 0;
}

// ---------------------------------------------------------
double BCIntegrate::LogEval(const std::vector<double> &x)
{
   // this method should better also be overloaded
   return log(Eval(x));
}

// ---------------------------------------------------------
double BCIntegrate::EvalSampling(const std::vector<double> & /*x*/)
{
   BCLog::OutWarning( "BCIntegrate::EvalSampling. No function. Function needs to be overloaded.");
   return 0;
}

// ---------------------------------------------------------
double BCIntegrate::LogEvalSampling(const std::vector<double> &x)
{
   return log(EvalSampling(x));
}

// ---------------------------------------------------------
void BCIntegrate::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
   fIntegrationMethod = method;
}

// ---------------------------------------------------------
int BCIntegrate::IntegrateResetResults()
{
   fBestFitParameters.clear();
   fBestFitParameterErrors.clear();
   fBestFitParametersMarginalized.clear();

   // no errors
   return 1;
}

// ---------------------------------------------------------
int BCIntegrate::GetNIntegrationVariables() {
	 int n = 0;
   for(int i = 0; i < fNvar; i++)
			if (fVarlist[i])
				 n++;
	 return n;
}

// ---------------------------------------------------------
double BCIntegrate::CalculateIntegrationVolume() {
	fIntegrationVolume=-1.;

   for(int i = 0; i < fNvar; i++)
		 if(fVarlist[i]) {
			 if (fIntegrationVolume<0)
				 fIntegrationVolume = 1;
			 fIntegrationVolume *= fMax[i] - fMin[i];
		 }

	 if (fIntegrationVolume<0)
		 fIntegrationVolume = 0;

	 return fIntegrationVolume;
}

// ---------------------------------------------------------
double BCIntegrate::Integrate()
{
	return Integrate(fIntegrationMethod);
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type)
{
   std::vector<double> parameter;
   parameter.assign(fNvar, 0.0);
	 return Integrate(type,parameter);
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type, const std::vector<double> &x)
{
   switch(type)
		 {

			 // Monte Carlo Integration
		 case BCIntegrate::kIntMonteCarlo:

			 CalculateIntegrationVolume();
			 return Integrate(kIntMonteCarlo,
												&BCIntegrate::GetRandomVectorInParameterSpaceFixed,
												&BCIntegrate::EvaluatorMC,
												&BCIntegrate::IntegralUpdaterMC,
												std::vector<double>(2,0.0),x);
			 
			 
			 // Metropolis Integration
		 case BCIntegrate::kIntMetropolis:

			 InitMetroFixed(x);
			 return Integrate(kIntMetropolis,
												&BCIntegrate::GetRandomPointSamplingMetroFixed,
												&BCIntegrate::EvaluatorMetro,
												&BCIntegrate::IntegralUpdaterMetro,
												std::vector<double>(2,0.0),x);


			 // Importance Sampling Integration
		 case BCIntegrate::kIntImportance:

			 return Integrate(kIntImportance,
												&BCIntegrate::GetRandomPointImportanceFixed,
												&BCIntegrate::EvaluatorMetro,	 // use same evaluator as for metropolis
												&BCIntegrate::IntegralUpdaterMetro,	 // use same updater as for metropolis
												std::vector<double>(2,0.),x);
			 
			 
		 case BCIntegrate::kIntCuba:
#ifdef HAVE_CUBA_H
			 return IntegrateCuba(x);
#else
			 BCLog::OutError("!!! This version of BAT is compiled without Cuba.");
			 BCLog::OutError("    Use other integration methods or install Cuba and recompile BAT.");
			 break;
#endif
		 default:
			 BCLog::OutError(
											 Form("BCIntegrate::Integrate : Invalid integration method: %d", fIntegrationMethod));
			 break;
		 }
	 
   return 0;
}

// ---------------------------------------------------------
void BCIntegrate::SetBestFitParameters(const std::vector<double> &x) {
			fBestFitParameters.clear();
			for (int i=0; i<fNvar; i++)
						fBestFitParameters.push_back(x.at(i));
}

// ---------------------------------------------------------
void BCIntegrate::SetBestFitParameters(const std::vector<double> &x, const double &new_value, double &old_value) {
			if (new_value < old_value)
						return;
			old_value = new_value;
			SetBestFitParameters(x);
}

// ---------------------------------------------------------
int BCIntegrate::IntegrationOutputFrequency()
{
	if (fNIterationsOutput>0)
		return fNIterationsOutput;
	int nwrite = fNIterationsMax/10;
	if(nwrite < 10000)
		return 1000;
	if(nwrite < 100000)
		return 10000;
	if(nwrite < 1000000)
		return 100000;
	return 1000000;
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtStartOfIntegration(BCIntegrationMethod type, BCCubaMethod cubatype, const std::vector<double> &x) {
	 int NVarNow = GetNIntegrationVariables();

   BCLog::LogLevel level=BCLog::summary;

   if(fNvar!=NVarNow) {

      level=BCLog::detail;

			if (type==kIntCuba)
				 BCLog::OutDetail(Form("Running %s (%s) integration over %i dimensions out of %i.", 
															 DumpIntegrationMethod(type).c_str(),
															 DumpCubaIntegrationMethod(cubatype).c_str(),
															 NVarNow, fNvar));
			else
				 BCLog::OutDetail(Form("Running %s integration over %i dimensions out of %i.", 
															 DumpIntegrationMethod(type).c_str(),
															 NVarNow, fNvar));

			BCLog::OutDetail(" --> Fixed parameters:");
      for(int i = 0; i < fNvar; i++)
         if(!fVarlist[i])
            BCLog::OutDetail(Form("      %3i :  %g", i, x[i]));
   }
   else {
 			if (type==kIntCuba)
				 BCLog::OutDetail(Form("Running %s (%s) integration over %i dimensions.", 
															 DumpIntegrationMethod(type).c_str(),
															 DumpCubaIntegrationMethod(cubatype).c_str(),
															 NVarNow));
			else
				 BCLog::OutDetail(Form("Running %s integration over %i dimensions.", 
															 DumpIntegrationMethod(type).c_str(),
															 NVarNow));
	 }

   BCLog::Out(level, Form(" --> Minimum number of iterations: %i", GetNIterationsMin()));
   BCLog::Out(level, Form(" --> Maximum number of iterations: %i", GetNIterationsMax()));
   BCLog::Out(level, Form(" --> Aimed relative precision:     %e", GetRelativePrecision()));
	 BCLog::Out(level, Form(" --> Aimed absolute precision:     %e", GetAbsolutePrecision()));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtEndOfIntegration(double integral, double absprecision, double relprecision, int nIterations) {
   BCLog::OutSummary(Form(" --> Result of integration:        %e +- %e   in %i iterations.", integral, absprecision, nIterations));
   BCLog::OutSummary(Form(" --> Obtained relative precision:  %e. ", relprecision));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtIntegrationStatusUpdate(BCIntegrationMethod type, double integral, double absprecision, int nIterations) {
	BCLog::OutDetail(Form("%s. Iteration %i, integral: %e +- %e.", DumpIntegrationMethod(type).c_str(), nIterations, integral, absprecision));
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type, tRandomizer randomizer, tEvaluator evaluator, tIntegralUpdater updater,
													 std::vector<double> sums, const std::vector<double> &x) {
	 LogOutputAtStartOfIntegration(type, NCubaMethods, x);

   // reset variables
   double pmax = 0.;
   double integral = 0.;
   double absprecision = 2.*fAbsolutePrecision;
   double relprecision = 2.*fRelativePrecision;

   std::vector<double> randx (fNvar, 0.);

   // how often to print out the info line to screen
   int nwrite = IntegrationOutputFrequency();

   // reset number of iterations
   fNIterations = 0;

   // iterate while number of iterations is lower than minimum number of iterations
	 // or precision is not reached and the number of iterations is lower than maximum number of iterations 
   while ( ( fRelativePrecision < relprecision && fAbsolutePrecision < absprecision && fNIterationsMax > fNIterations )
					 || fNIterations < fNIterationsMin ) {
			
			// increase number of iterations
			fNIterations++;

			// get random numbers
			(this->*randomizer)(randx,x);

      // evaluate function at sampled point
			// updating sums & checking for maximum probablity
			SetBestFitParameters(randx, (this->*evaluator)(sums,randx), pmax);

			// update precisions
      if (fNIterations%fNIterationsPrecisionCheck == 0) {
				(this->*updater)(sums,fNIterations,integral,absprecision);
         relprecision = absprecision / integral;
			}
			
			// write status
			if (fNIterations%nwrite == 0) {
         double temp_integral;
				 double temp_absprecision;
				 (this->*updater)(sums,fNIterations,temp_integral,temp_absprecision);
				 LogOutputAtIntegrationStatusUpdate(type,temp_integral,temp_absprecision,fNIterations);
      }
   }

	 // calculate integral
	 (this->*updater)(sums,fNIterations,integral,absprecision);
	 relprecision = absprecision / integral;
	 
   // print to log
	 LogOutputAtEndOfIntegration(integral,absprecision,relprecision,fNIterations);

	 fError = absprecision;
   return integral;
}

// ---------------------------------------------------------
double BCIntegrate::EvaluatorMC(std::vector<double> &sums, const std::vector<double> &point) {
	double value = Eval(point);

	// add value to sum and sum of squares
	sums[0] += value;
	sums[1] += value * value;

	return value;
}

// ---------------------------------------------------------
void BCIntegrate::IntegralUpdaterMC(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision) {
	integral = fIntegrationVolume * sums[0] / nIterations;
	absprecision = sqrt((1.0 / double(nIterations)) * (fIntegrationVolume * fIntegrationVolume * sums[1] / double(nIterations) - integral * integral));
}

double BCIntegrate::EvaluatorMetro(std::vector<double> &sums, const std::vector<double> &point)
{
	// calculate probability at random point
	double val_f = Eval(point);

	// calculate sampling distributions at that point
	double val_g = EvalSampling(point);

	// add ratio to sum and sum of squares
	if (val_g > 0) {
		sums[0] += val_f / val_g;
		sums[1] += val_f * val_f / val_g / val_g;
	}

	return val_f;
}

void BCIntegrate::IntegralUpdaterMetro(const std::vector<double> &sums, const int &nIterations, double &integral, double &absprecision) {
	integral = sums[0] / nIterations;
	absprecision = sqrt((1.0 / double(nIterations)) * (sums[1] / double(nIterations) - integral * integral));
}

// ---------------------------------------------------------
TH1D* BCIntegrate::Marginalize(BCParameter * parameter)
{
   BCLog::OutDebug(Form(" --> Marginalizing model wrt. parameter %s using %s.", parameter->GetName().data(), DumpMarginalizationMethod().c_str()));

   switch(fMarginalizationMethod)
   {
      case BCIntegrate::kMargMonteCarlo:
         return MarginalizeByIntegrate(parameter);

      case BCIntegrate::kMargMetropolis:
         return MarginalizeByMetro(parameter);

      default:
         BCLog::OutError(
            Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizationMethod));
         break;

   }

   return 0;
}

// ---------------------------------------------------------
TH2D * BCIntegrate::Marginalize(BCParameter * parameter1, BCParameter * parameter2)
{
   switch(fMarginalizationMethod)
   {
      case BCIntegrate::kMargMonteCarlo:
         return MarginalizeByIntegrate(parameter1,parameter2);

      case BCIntegrate::kMargMetropolis:
         return MarginalizeByMetro(parameter1,parameter2);

      default:
         BCLog::OutError(
            Form("BCIntegrate::Marginalize. Invalid marginalization method: %d. Return 0.", fMarginalizationMethod));
         break;
   }

   return 0;
}

// ---------------------------------------------------------
TH1D* BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter)
{
   // set parameter to marginalize
   ResetVarlist(1);
   int index = parameter->GetIndex();
   UnsetVar(index);

   // define histogram
   double hmin = parameter->GetLowerLimit();
   double hmax = parameter->GetUpperLimit();
   double hdx  = (hmax - hmin) / double(fNbins);
   TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

   // integrate
   std::vector<double> randx;
   randx.assign(fNvar, 0.);

   for(int i=0;i<=fNbins;i++) {
      randx[index] = hmin + (double)i * hdx;

      double val = Integrate(kIntMonteCarlo,randx);
      hist->Fill(randx[index], val);
   }

   // normalize
   hist->Scale( 1./hist->Integral("width") );

   return hist;
}

// ---------------------------------------------------------
TH2D * BCIntegrate::MarginalizeByIntegrate(BCParameter * parameter1, BCParameter * parameter2)
{
   // set parameter to marginalize
   ResetVarlist(1);
   int index1 = parameter1->GetIndex();
   UnsetVar(index1);
   int index2 = parameter2->GetIndex();
   UnsetVar(index2);

   // define histogram
   double hmin1 = parameter1->GetLowerLimit();
   double hmax1 = parameter1->GetUpperLimit();
   double hdx1  = (hmax1 - hmin1) / double(fNbins);

   double hmin2 = parameter2->GetLowerLimit();
   double hmax2 = parameter2->GetUpperLimit();
   double hdx2  = (hmax2 - hmin2) / double(fNbins);

   TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1->GetName().data(), parameter2->GetName().data()),"",
         fNbins, hmin1, hmax1,
         fNbins, hmin2, hmax2);

   // integrate
   std::vector<double> randx;
   randx.assign(fNvar, 0.0);

   for(int i=0;i<=fNbins;i++) {
      randx[index1] = hmin1 + (double)i * hdx1;
      for(int j=0;j<=fNbins;j++) {
         randx[index2] = hmin2 + (double)j * hdx2;

         double val = Integrate(kIntMonteCarlo,randx);
         hist->Fill(randx[index1],randx[index2], val);
      }
   }

   // normalize
   hist->Scale(1.0/hist->Integral("width"));

   return hist;
}

// ---------------------------------------------------------
TH1D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter)
{
   int niter = fMarkovChainNIterations;

   if (fMarkovChainAutoN == true)
      niter = fNbins*fNbins*fNSamplesPer2DBin*fNvar;

   BCLog::OutDetail(Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

   // set parameter to marginalize
   int index = parameter->GetIndex();

   // define histogram
   double hmin = parameter->GetLowerLimit();
   double hmax = parameter->GetUpperLimit();
   TH1D * hist = new TH1D("hist","", fNbins, hmin, hmax);

   // prepare Metro
   std::vector<double> randx;
   randx.assign(fNvar, 0.0);
   InitMetro();

   for(int i=0;i<=niter;i++) {
      GetRandomPointMetro(randx);
      hist->Fill(randx[index]);
   }

   // normalize
   hist->Scale(1.0/hist->Integral("width"));

   return hist;
}

// ---------------------------------------------------------
TH2D * BCIntegrate::MarginalizeByMetro(BCParameter * parameter1, BCParameter * parameter2)
{
   int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

   // set parameter to marginalize
   int index1 = parameter1->GetIndex();
   int index2 = parameter2->GetIndex();

   // define histogram
   double hmin1 = parameter1->GetLowerLimit();
   double hmax1 = parameter1->GetUpperLimit();

   double hmin2 = parameter2->GetLowerLimit();
   double hmax2 = parameter2->GetUpperLimit();

   TH2D * hist = new TH2D(Form("hist_%s_%s", parameter1->GetName().data(), parameter2->GetName().data()),"",
         fNbins, hmin1, hmax1,
         fNbins, hmin2, hmax2);

   // prepare Metro
   std::vector<double> randx;
   randx.assign(fNvar, 0.0);
   InitMetro();

   for(int i=0;i<=niter;i++) {
      GetRandomPointMetro(randx);
      hist->Fill(randx[index1],randx[index2]);
   }

   // normalize
   hist->Scale(1.0/hist->Integral("width"));

   return hist;
}

// ---------------------------------------------------------
int BCIntegrate::MarginalizeAllByMetro(const char * name="")
{
   int niter=fNbins*fNbins*fNSamplesPer2DBin*fNvar;

   BCLog::OutDetail(Form(" --> Number of samples in Metropolis marginalization: %d.", niter));

   // define 1D histograms
   for(int i=0;i<fNvar;i++)
   {
      double hmin1 = fx->at(i)->GetLowerLimit();
      double hmax1 = fx->at(i)->GetUpperLimit();

      TH1D * h1 = new TH1D(
            TString::Format("h%s_%s", name, fx->at(i)->GetName().data()),"",
            fNbins, hmin1, hmax1);

      fHProb1D.push_back(h1);
   }

   // define 2D histograms
   for(int i=0;i<fNvar-1;i++)
      for(int j=i+1;j<fNvar;j++) {
         double hmin1 = fx->at(i)->GetLowerLimit();
         double hmax1 = fx->at(i)->GetUpperLimit();

         double hmin2 = fx->at(j)->GetLowerLimit();
         double hmax2 = fx->at(j)->GetUpperLimit();

         TH2D * h2 = new TH2D(
            TString::Format("h%s_%s_%s", name, fx->at(i)->GetName().data(), fx->at(j)->GetName().data()),"",
            fNbins, hmin1, hmax1,
            fNbins, hmin2, hmax2);

         fHProb2D.push_back(h2);
      }

   // get number of 2d distributions
   int nh2d = fHProb2D.size();

   BCLog::OutDetail(Form(" --> Marginalizing %d 1D distributions and %d 2D distributions.", fNvar, nh2d));

   // prepare function fitting
   double dx = 0.;
   double dy = 0.;

   if (fFitFunctionIndexX >= 0) {
      dx = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX) -
            fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX))
            / double(fErrorBandNbinsX);

      dx = (fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY) -
            fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY))
            / double(fErrorBandNbinsY);

      fErrorBandXY = new TH2D(
            TString::Format("errorbandxy_%d",BCLog::GetHIndex()), "",
            fErrorBandNbinsX,
            fDataPointLowerBoundaries->GetValue(fFitFunctionIndexX) - 0.5 * dx,
            fDataPointUpperBoundaries->GetValue(fFitFunctionIndexX) + 0.5 * dx,
            fErrorBandNbinsY,
            fDataPointLowerBoundaries->GetValue(fFitFunctionIndexY) - 0.5 * dy,
            fDataPointUpperBoundaries->GetValue(fFitFunctionIndexY) + 0.5 * dy);
      fErrorBandXY->SetStats(kFALSE);

      for (int ix = 1; ix <= fErrorBandNbinsX; ++ix)
         for (int iy = 1; iy <= fErrorBandNbinsX; ++iy)
            fErrorBandXY->SetBinContent(ix, iy, 0.0);
   }

   // prepare Metro
   std::vector<double> randx;

   randx.assign(fNvar, 0.0);
   InitMetro();

   // reset counter
   fNacceptedMCMC = 0;

   // run Metro and fill histograms
   for(int i=0;i<=niter;i++) {
      GetRandomPointMetro(randx);

      // save this point to the markov chain in the ROOT file
      if (fFlagWriteMarkovChain) {
         fMCMCIteration = i;
         fMarkovChainTree->Fill();
      }

      for(int j=0;j<fNvar;j++)
         fHProb1D[j]->Fill( randx[j] );

      int ih=0;
      for(int j=0;j<fNvar-1;j++)
         for(int k=j+1;k<fNvar;k++) {
            fHProb2D[ih]->Fill(randx[j],randx[k]);
            ih++;
         }

      if((i+1)%100000==0)
         BCLog::OutDebug(Form("BCIntegrate::MarginalizeAllByMetro. %d samples done.",i+1));

      // function fitting

      if (fFitFunctionIndexX >= 0) {
         // loop over all possible x values ...
         if (fErrorBandContinuous) {
            double x = 0;

            for (int ix = 0; ix < fErrorBandNbinsX; ix++) {
               // calculate x
               x = fErrorBandXY->GetXaxis()->GetBinCenter(ix + 1);

               // calculate y
               std::vector<double> xvec;
               xvec.push_back(x);
               double y = FitFunction(xvec, randx);
               xvec.clear();

               // fill histogram
               fErrorBandXY->Fill(x, y);
            }
         }

         // ... or evaluate at the data point x-values
         else {
            int ndatapoints = int(fErrorBandX.size());
            double x = 0;

            for (int ix = 0; ix < ndatapoints; ++ix) {
               // calculate x
               x = fErrorBandX.at(ix);

               // calculate y
               std::vector<double> xvec;
               xvec.push_back(x);
               double y = FitFunction(xvec, randx);
               xvec.clear();

               // fill histogram
               fErrorBandXY->Fill(x, y);
            }
         }
      }
   }

   // normalize histograms
   for(int i=0;i<fNvar;i++)
      fHProb1D[i]->Scale( 1./fHProb1D[i]->Integral("width") );

   for (int i=0;i<nh2d;i++)
      fHProb2D[i]->Scale( 1./fHProb2D[i]->Integral("width") );

   if (fFitFunctionIndexX >= 0)
      fErrorBandXY->Scale(1.0/fErrorBandXY->Integral() * fErrorBandXY->GetNbinsX());

   if (fFitFunctionIndexX >= 0) {
      for (int ix = 1; ix <= fErrorBandNbinsX; ix++) {
         double sum = 0;

         for (int iy = 1; iy <= fErrorBandNbinsY; iy++)
            sum += fErrorBandXY->GetBinContent(ix, iy);

         for (int iy = 1; iy <= fErrorBandNbinsY; iy++) {
            double newvalue = fErrorBandXY->GetBinContent(ix, iy) / sum;
            fErrorBandXY->SetBinContent(ix, iy, newvalue);
         }
      }
   }

   BCLog::OutDebug(Form("BCIntegrate::MarginalizeAllByMetro done with %i trials and %i accepted trials. Efficiency is %f",fNmetro, fNacceptedMCMC, double(fNacceptedMCMC)/double(fNmetro)));

   return fNvar+nh2d;
}

// ---------------------------------------------------------
TH1D * BCIntegrate::GetH1D(int parIndex)
{
   if(fHProb1D.size()==0) {
      BCLog::OutWarning("BCModel::GetH1D. MarginalizeAll() has to be run prior to this to fill the distributions.");
      return 0;
   }

   if(parIndex<0 || parIndex>fNvar) {
      BCLog::OutWarning(Form("BCIntegrate::GetH1D. Parameter index %d is invalid.",parIndex));
      return 0;
   }

   return fHProb1D[parIndex];
}

// ---------------------------------------------------------
int BCIntegrate::GetH2DIndex(int parIndex1, int parIndex2)
{
   if(parIndex1>fNvar-1 || parIndex1<0) {
      BCLog::OutWarning(Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex1));
      return -1;
   }

   if(parIndex2>fNvar-1 || parIndex2<0) {
      BCLog::OutWarning(Form("BCIntegrate::GetH2DIndex. Parameter index %d is invalid", parIndex2));
      return -1;
   }

   if(parIndex1==parIndex2) {
      BCLog::OutWarning(Form("BCIntegrate::GetH2DIndex. Parameters have equal indeces: %d , %d", parIndex1, parIndex2));
      return -1;
   }

   if(parIndex1>parIndex2) {
      BCLog::OutWarning("BCIntegrate::GetH2DIndex. First parameters must be smaller than second (sorry).");
      return -1;
   }

   int index=0;
   for(int i=0;i<fNvar-1;i++)
      for(int j=i+1;j<fNvar;j++) {
         if(i==parIndex1 && j==parIndex2)
            return index;
         index++;
      }

   BCLog::OutWarning("BCIntegrate::GetH2DIndex. Invalid index combination.");

   return -1;
}

// ---------------------------------------------------------
TH2D * BCIntegrate::GetH2D(int parIndex1, int parIndex2)
{
   if(fHProb2D.size()==0) {
      BCLog::OutWarning("BCModel::GetH2D. MarginalizeAll() has to be run prior to this to fill the distributions.");
      return 0;
   }

   int hindex = GetH2DIndex(parIndex1, parIndex2);
   if(hindex==-1)
      return 0;

   if(hindex>(int)fHProb2D.size()-1) {
      BCLog::OutWarning("BCIntegrate::GetH2D. Got invalid index from GetH2DIndex(). Something went wrong.");
      return 0;
   }

   return fHProb2D[hindex];
}

// ---------------------------------------------------------
double BCIntegrate::GetRandomPoint(std::vector<double> &x)
{
   GetRandomVectorInParameterSpace(x);
   return Eval(x);
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomPointImportance(std::vector<double> &x)
{
   double p = 1.1;
   double g = 1.0;

   while (p>g) {
      GetRandomVectorInParameterSpace(x);

      p = fRandom->Rndm();

      g = EvalSampling(x);
   }
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomPointImportanceFixed(std::vector<double> &x, const std::vector<double> &fix_x)
{
   double p = 1.1;
   double g = 1.0;

   while (p>g) {
		 GetRandomVectorInParameterSpaceFixed(x,fix_x);

      p = fRandom->Rndm();

      g = EvalSampling(x);
   }
}

// ---------------------------------------------------------
void BCIntegrate::InitMetro()
{
	fNmetro=0;
	
	if (fXmetro0.size() <= 0)
		// start in the center of the phase space
		for(int i=0;i<fNvar;i++)
			fXmetro0.push_back((fMin[i]+fMax[i])/2.0);
	
	fMarkovChainValue =  LogEval(fXmetro0);

	// run metropolis for a few times and dump the result... (to forget the initial position)
	std::vector<double> x (fNvar, 0.0);
												 
	for(int i=0;i<1000;i++)
		GetRandomPointMetro(x);

	fNmetro = 0;
}

// ---------------------------------------------------------
void BCIntegrate::InitMetroFixed(const std::vector<double> &fix_x)
{
	fNmetro=0;

	if (fXmetro0.size() <= 0)
		// start in the center of the phase space
		for(int i=0;i<fNvar;i++)
			if(fVarlist[i])
				fXmetro0.push_back((fMin[i]+fMax[i])/2.0);
			else
				fXmetro0.push_back(fix_x[i]);
	
	fMarkovChainValue =  LogEval(fXmetro0);
	
	// run metropolis for a few times and dump the result... (to forget the initial position)
	std::vector<double> x (fNvar, 0.0);
												 
	for(int i=0;i<1000;i++)
		GetRandomPointMetroFixed(x,fix_x);

	fNmetro = 0;
}

// ---------------------------------------------------------
bool BCIntegrate::UpdatePointsMetro(std::vector<double> &old_point, double p_old,
																		std::vector<double> &new_point, double p_new,
																		std::vector<double> &x, bool in)
{
	bool accept = false;

	if (in)
		// compare log probabilities
		if(p_new>=p_old)
			accept = true;
		else {
			double r = log(fRandom->Rndm());
			if(r<p_new-p_old)
				accept = true;
		}

	// update points after the decision
	if(accept)
		for(int i=0;i<fNvar;i++) {
			old_point[i] = new_point[i];
			x[i]         = new_point[i];
		}
	else
		for(int i=0;i<fNvar;i++) {
			new_point[i] = old_point[i];
			x[i]         = old_point[i];
		}

	return accept;
}



// ---------------------------------------------------------
void BCIntegrate::GetRandomPointMetro(std::vector<double> &x)
{
	// get new point and check if in allowed region
	bool in = GetRandomVectorMetroScaled(fXmetro1,fXmetro0);

	// calculate the log probability
	double p1  = 0;
	if(in)
		p1 = LogEval(fXmetro1);

	if(UpdatePointsMetro(fXmetro0,fMarkovChainValue,fXmetro1,p1,x,in)) {
		// increase counter
		fNacceptedMCMC++;
		fMarkovChainValue = p1;
	}
	else
		fMarkovChainValue = fMarkovChainValue;

	fNmetro++;
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomPointSamplingMetro(std::vector<double> &x)
{
	// get new point and check if in allowed region
	bool in = GetRandomVectorMetroScaled(fXmetro1,fXmetro0);

	// calculate the log probability
	double p1  = 0;
	if(in)
		p1 = LogEvalSampling(fXmetro1);

	UpdatePointsMetro(fXmetro0,LogEvalSampling(fXmetro0),fXmetro1,p1,x,in);

	fNmetro++;
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomPointMetroFixed(std::vector<double> &x, const std::vector<double> &/*fix_x*/)
{
	// get new point and check if in allowed region
	bool in = GetRandomVectorMetroScaledFixed(fXmetro1,fXmetro0);

	// calculate the log probability
	double p1  = 0;
	if(in)
		p1 = LogEval(fXmetro1);

	if(UpdatePointsMetro(fXmetro0,fMarkovChainValue,fXmetro1,p1,x,in)) {
		// increase counter
		fNacceptedMCMC++;
		fMarkovChainValue = p1;
	}
	else
		fMarkovChainValue = fMarkovChainValue;

	fNmetro++;
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomPointSamplingMetroFixed(std::vector<double> &x, const std::vector<double> &/*fix_x*/)
{
	// get new point and check if in allowed region
	bool in = GetRandomVectorMetroScaledFixed(fXmetro1,fXmetro0);

	// calculate the log probability
	double p1  = 0;
	if(in)
		p1 = LogEvalSampling(fXmetro1);

	UpdatePointsMetro(fXmetro0,LogEvalSampling(fXmetro0),fXmetro1,p1,x,in);

	fNmetro++;
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorUnitHypercube(std::vector<double> &x)
{
   double * randx = new double[fNvar];

   fRandom->RndmArray(fNvar, randx);

   for(int i=0;i<fNvar;i++)
      x[i] = randx[i];

   delete[] randx;
   randx = 0;
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorInParameterSpace(std::vector<double> &x) {
	 GetRandomVectorUnitHypercube(x);
	 
	 for (int i=0; i<fNvar; i++)
			x[i] = fMin[i] + x[i]*(fMax[i]-fMin[i]);
}


// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorInParameterSpaceFixed(std::vector<double> &x, const std::vector<double> &fix_x) {
	 GetRandomVectorUnitHypercube(x);
	 
	 for (int i=0; i<fNvar; i++)
			if(fVarlist[i])
				 x[i] = fMin[i] + x[i]*(fMax[i]-fMin[i]);
			else
				 x[i] = fix_x[i];
}


// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorMetro(std::vector<double> &x)
{
   double * randx = new double[fNvar];

   fRandom->RndmArray(fNvar, randx);

   for(int i=0;i<fNvar;i++)
      x[i] = 2.0 * (randx[i] - 0.5) * fMarkovChainStepSize;

   delete[] randx;
   randx = 0;
}
 
// ---------------------------------------------------------
bool BCIntegrate::GetRandomVectorMetroScaled(std::vector<double> &x, const std::vector<double> &prev_x)
{
	GetRandomVectorMetro(x);

	bool in = true;

	for(int i=0;i<fNvar;i++) {	
		x[i] = prev_x[i] + x[i] * (fMax[i]-fMin[i]);
		// check whether the generated point is inside the allowed region
		if( x[i]<fMin[i] || x[i]>fMax[i] )
			in=false;
	}
	return in;
}

// ---------------------------------------------------------
bool BCIntegrate::GetRandomVectorMetroScaledFixed(std::vector<double> &x, const std::vector<double> &prev_x)
{
	GetRandomVectorMetro(x);

	bool in = true;

	for(int i=0;i<fNvar;i++) {
		 if(fVarlist[i])
			 x[i] = prev_x[i] + x[i] * (fMax[i]-fMin[i]);
		 else
			 x[i] = prev_x[i];
		 // check whether the generated point is inside the allowed region
		 if( x[i]<fMin[i] || x[i]>fMax[i] )
			 in=false;
	}

	return in;
}


// ---------------------------------------------------------
TMinuit * BCIntegrate::GetMinuit()
{
   if (!fMinuit)
      fMinuit = new TMinuit();

   return fMinuit;
}

// ---------------------------------------------------------

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

   // set print level
   fMinuit->SetPrintLevel(printlevel);

   // set function
   fMinuit->SetFCN(&BCIntegrate::FCNLikelihood);

   // set UP for likelihood
   fMinuit->SetErrorDef(0.5);

   // set parameters
   int flag;
   for (int i = 0; i < fNvar; i++) {
      double starting_point = (fMin[i]+fMax[i])/2.;
      if(have_start)
         starting_point = start[i];
         fMinuit->mnparm(i,
                         fx->at(i)->GetName().data(),
                         starting_point,
                         (fMax[i]-fMin[i])/100.,
                         fMin[i],
                         fMax[i],
                         flag);
   }

   // do mcmc minimization
   //   fMinuit->mnseek();

   // do minimization
   fMinuit->mnexcm("MIGRAD", fMinuitArglist, 2, flag);

   // improve search for local minimum
   //   fMinuit->mnimpr();

   // copy flag
   fMinuitErrorFlag = flag;

   // check if mode found by minuit is better than previous estimation
   double probmax = 0;
   bool valid = false;

   if ( int(fBestFitParameters.size()) == fNvar) {
      valid = true;
      for (int i = 0; i < fNvar; ++i)
         if (fBestFitParameters.at(i) < fMin[i] || fBestFitParameters.at(i) > fMax[i])
            valid= false;

      if (valid)
         probmax = Eval( fBestFitParameters );
   }

   std::vector<double> tempvec;
   for (int i = 0; i < fNvar; i++) {
      double par;
      double parerr;
      fMinuit->GetParameter(i, par, parerr);
      tempvec.push_back(par);
   }
   double probmaxminuit = Eval( tempvec );

   // set best fit parameters
   if ( (probmaxminuit > probmax && !fFlagIgnorePrevOptimization) || !valid || fFlagIgnorePrevOptimization) {
      fBestFitParameters.clear();
      fBestFitParameterErrors.clear();

      for (int i = 0; i < fNvar; i++) {
         double par;
         double parerr;
         fMinuit->GetParameter(i, par, parerr);
         fBestFitParameters.push_back(par);
         fBestFitParameterErrors.push_back(parerr);
      }

      // set optimization method used to find the mode
      fOptimizationMethodMode = BCIntegrate::kOptMinuit;
   }

   // delete minuit
   delete fMinuit;
   fMinuit = 0;

   return;
}

// ---------------------------------------------------------
void BCIntegrate::InitializeSATree()
{
  if (fTreeSA)
    delete fTreeSA;
  fTreeSA = new TTree("SATree", "SATree");

   fTreeSA->Branch("Iteration",      &fSANIterations,   "iteration/I");
   fTreeSA->Branch("NParameters",    &fNvar,            "parameters/I");
   fTreeSA->Branch("Temperature",    &fSATemperature,   "temperature/D");
   fTreeSA->Branch("LogProbability", &fSALogProb,       "log(probability)/D");

   for (int i = 0; i < fNvar; ++i)
      fTreeSA->Branch(TString::Format("Parameter%i", i), &fSAx[i], TString::Format("parameter %i/D", i));
}

// ---------------------------------------------------------
void BCIntegrate::FindModeSA(std::vector<double> start)
{
   // note: if f(x) is the function to be minimized, then
   // f(x) := - LogEval(parameters)

   bool have_start = true;
   std::vector<double> x, y, best_fit; // vectors for current state, new proposed state and best fit up to now
   double fval_x, fval_y, fval_best_fit; // function values at points x, y and best_fit (we save them rather than to re-calculate them every time)
   int t = 1; // time iterator

   // check start values
   if (int(start.size()) != fNvar)
      have_start = false;

   // if no starting point is given, set to center of parameter space
   if ( !have_start ) {
      start.clear();
      for (int i = 0; i < fNvar; i++)
         start.push_back((fMin[i]+fMax[i])/2.);
   }

   // set current state and best fit to starting point
   x.clear();
   best_fit.clear();
   for (int i = 0; i < fNvar; i++) {
      x.push_back(start[i]);
      best_fit.push_back(start[i]);
   }
   // calculate function value at starting point
   fval_x = fval_best_fit = LogEval(x);

   // run while still "hot enough"
   while ( SATemperature(t) > fSATmin ) {
      // generate new state
      y = GetProposalPointSA(x, t);

      // check if the proposed point is inside the phase space
      // if not, reject it
      bool is_in_ranges = true;

      for (int i = 0; i < fNvar; i++)
         if (y[i] > fMax[i] || y[i] < fMin[i])
            is_in_ranges = false;

      if ( !is_in_ranges )
         ; // do nothing...
      else {
         // calculate function value at new point
         fval_y = LogEval(y);

         // is it better than the last one?
         // if so, update state and chef if it is the new best fit...
         if (fval_y >= fval_x) {
            x.clear();
            for (int i = 0; i < fNvar; i++)
               x.push_back(y[i]);

            fval_x = fval_y;

            if (fval_y > fval_best_fit) {
               best_fit.clear();
               for (int i = 0; i < fNvar; i++)
                  best_fit.push_back(y[i]);

               fval_best_fit = fval_y;
            }
         }
         // ...else, only accept new state w/ certain probability
         else {
            if (fRandom->Rndm() <= exp( (fval_y - fval_x) / SATemperature(t) )) {
               x.clear();
               for (int i = 0; i < fNvar; i++)
                  x.push_back(y[i]);

               fval_x = fval_y;
            }
         }
      }

      // update tree variables
      fSANIterations = t;
      fSATemperature = SATemperature(t);
      fSALogProb = fval_x;
      fSAx.clear();
      for (int i = 0; i < fNvar; ++i)
         fSAx.push_back(x[i]);

      // fill tree
      if (fFlagWriteSAToFile)
         fTreeSA->Fill();

      // increate t
      t++;
   }

   // check if mode found by minuit is better than previous estimation
   double probmax = 0;
   bool valid = false;

   if ( int(fBestFitParameters.size()) == fNvar) {
      valid = true;
      for (int i = 0; i < fNvar; ++i)
         if (fBestFitParameters.at(i) < fMin[i] || fBestFitParameters.at(i) > fMax[i])
            valid= false;

      if (valid)
         probmax = Eval( fBestFitParameters );
   }

   double probmaxsa = Eval( best_fit);

   if ( (probmaxsa > probmax  && !fFlagIgnorePrevOptimization) || !valid || fFlagIgnorePrevOptimization) {
      // set best fit parameters
      fBestFitParameters.clear();
      fBestFitParameterErrors.clear();

      for (int i = 0; i < fNvar; i++) {
         fBestFitParameters.push_back(best_fit[i]);
         fBestFitParameterErrors.push_back(-1.); // error undefined
      }

      // set optimization moethod used to find the mode
      fOptimizationMethodMode = BCIntegrate::kOptSA;
   }

   return;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperature(double t)
{
   // do we have Cauchy (default), Boltzmann or custom annealing schedule?
   if (fSASchedule == BCIntegrate::kSABoltzmann)
      return SATemperatureBoltzmann(t);
   else if (fSASchedule == BCIntegrate::kSACauchy)
      return SATemperatureCauchy(t);
   else
      return SATemperatureCustom(t);
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureBoltzmann(double t)
{
   return fSAT0 / log((double)(t + 1));
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCauchy(double t)
{
   return fSAT0 / (double)t;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCustom(double /*t*/)
{
   BCLog::OutError("BCIntegrate::SATemperatureCustom : No custom temperature schedule defined");
   return 0.;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSA(const std::vector<double> &x, int t)
{
   // do we have Cauchy (default), Boltzmann or custom annealing schedule?
   if (fSASchedule == BCIntegrate::kSABoltzmann)
      return GetProposalPointSABoltzmann(x, t);
   else if (fSASchedule == BCIntegrate::kSACauchy)
      return GetProposalPointSACauchy(x, t);
   else
      return GetProposalPointSACustom(x, t);
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSABoltzmann(const std::vector<double> &x, int t)
{
   std::vector<double> y;
   y.clear();
   double new_val, norm;

   for (int i = 0; i < fNvar; i++) {
      norm = (fMax[i] - fMin[i]) * SATemperature(t) / 2.;
      new_val = x[i] + norm * fRandom->Gaus();
      y.push_back(new_val);
   }

   return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACauchy(const std::vector<double> &x, int t)
{
   std::vector<double> y;
   y.clear();

   if (fNvar == 1) {
      double cauchy, new_val, norm;

      norm = (fMax[0] - fMin[0]) * SATemperature(t) / 2.;
      cauchy = tan(3.14159 * (fRandom->Rndm() - 0.5));
      new_val = x[0] + norm * cauchy;
      y.push_back(new_val);
   }
   else {
      // use sampling to get radial n-dim Cauchy distribution

      // first generate a random point uniformly distributed on a
      // fNvar-dimensional hypersphere
      y = SAHelperGetRandomPointOnHypersphere();

      // scale the vector by a random factor determined by the radial
      // part of the fNvar-dimensional Cauchy distribution
      double radial = SATemperature(t) * SAHelperGetRadialCauchy();

      // scale y by radial part and the size of dimension i in phase space
      // afterwards, move by x
      for (int i = 0; i < fNvar; i++)
         y[i] = (fMax[i] - fMin[i]) * y[i] * radial / 2. + x[i];
   }

   return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACustom(const std::vector<double> & /*x*/, int /*t*/)
{
   BCLog::OutError("BCIntegrate::GetProposalPointSACustom : No custom proposal function defined");
   return std::vector<double>(fNvar);
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::SAHelperGetRandomPointOnHypersphere()
{
   std::vector<double> rand_point(fNvar);

   // This method can only be called with fNvar >= 2 since the 1-dim case
   // is already hard wired into the Cauchy annealing proposal function.
   // To speed things up, hard-code fast method for 2 and dimensions.
   // The algorithm for 2D can be found at
   // http://mathworld.wolfram.com/CirclePointPicking.html
   // For 3D just using ROOT's algorithm.
   if (fNvar == 2) {
      double x1, x2, s;
      do {
         x1 = fRandom->Rndm() * 2. - 1.;
         x2 = fRandom->Rndm() * 2. - 1.;
         s = x1*x1 + x2*x2;
      }
      while (s >= 1);

      rand_point[0] = (x1*x1 - x2*x2) / s;
      rand_point[1] = (2.*x1*x2) / s;
   }
   else if (fNvar == 3) {
      fRandom->Sphere(rand_point[0], rand_point[1], rand_point[2], 1.0);
   }
   else {
      double s = 0.,
         gauss_num;

      for (int i = 0; i < fNvar; i++) {
         gauss_num = fRandom->Gaus();
         rand_point[i] = gauss_num;
         s += gauss_num * gauss_num;
      }
      s = sqrt(s);

      for (int i = 0; i < fNvar; i++)
         rand_point[i] = rand_point[i] / s;
   }

   return rand_point;
}

// ---------------------------------------------------------
double BCIntegrate::SAHelperGetRadialCauchy()
{
   // theta is sampled from a rather complicated distribution,
   // so first we create a lookup table with 10000 random numbers
   // once and then, each time we need a new random number,
   // we just look it up in the table.
   double theta;

   // static vectors for theta-sampling-map
   static double map_u[10001];
   static double map_theta[10001];
   static bool initialized = false;
   static int map_dimension = 0;

   // is the lookup-table already initialized? if not, do it!
   if (!initialized || map_dimension != fNvar) {
      double init_theta;
      double beta = SAHelperSinusToNIntegral(fNvar - 1, 1.57079632679);

      for (int i = 0; i <= 10000; i++) {
         init_theta = 3.14159265 * (double)i / 5000.;
         map_theta[i] = init_theta;

         map_u[i] = SAHelperSinusToNIntegral(fNvar - 1, init_theta) / beta;
      }

      map_dimension = fNvar;
      initialized = true;
   } // initializing is done.

   // generate uniform random number for sampling
   double u = fRandom->Rndm();

   // Find the two elements just greater than and less than u
   // using a binary search (O(log(N))).
   int lo = 0;
   int up = 10000;
   int mid;

   while (up != lo) {
      mid = ((up - lo + 1) / 2) + lo;

      if (u >= map_u[mid])
         lo = mid;
      else
         up = mid - 1;
   }
   up++;

   // perform linear interpolation:
   theta = map_theta[lo] + (u - map_u[lo]) / (map_u[up] - map_u[lo]) * (map_theta[up] - map_theta[lo]);

   return tan(theta);
}

// ---------------------------------------------------------
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
         * SAHelperSinusToNIntegral(dim - 2, theta);
}

// ---------------------------------------------------------

void BCIntegrate::SetMode(std::vector<double> mode)
{
   if((int)mode.size() == fNvar) {
      fBestFitParameters.clear();
      for (int i = 0; i < fNvar; i++)
         fBestFitParameters.push_back(mode[i]);
   }
}

// ---------------------------------------------------------

void BCIntegrate::FCNLikelihood(int & /*npar*/, double * /*grad*/, double &fval, double * par, int /*flag*/)
{
   // copy parameters
   std::vector<double> parameters;

   int n = global_this->GetNvar();

   for (int i = 0; i < n; i++)
      parameters.push_back(par[i]);

   fval = - global_this->LogEval(parameters);

   // delete parameters
   parameters.clear();
}

// ---------------------------------------------------------

//void BCIntegrate::FindModeMCMC(int flag_run)
void BCIntegrate::FindModeMCMC()
{
   // call PreRun
//   if (flag_run == 0)
   if (!fMCMCFlagPreRun)
      MCMCMetropolisPreRun();

   // find global maximum
   //   double probmax = (MCMCGetMaximumLogProb()).at(0);

   double probmax = 0;

   if ( int(fBestFitParameters.size()) == fNvar) {
      probmax = Eval( fBestFitParameters );
//      fBestFitParameters = MCMCGetMaximumPoint(0);
   }

   // loop over all chains and find the maximum point
   for (int i = 0; i < fMCMCNChains; ++i) {
      double prob = exp( (MCMCGetMaximumLogProb()).at(i));

      // copy the point into the vector
      if ( (prob >= probmax && !fFlagIgnorePrevOptimization) || fFlagIgnorePrevOptimization) {
         probmax = prob;

         fBestFitParameters.clear();
         fBestFitParameterErrors.clear();
         fBestFitParameters = MCMCGetMaximumPoint(i);

         for (int j = 0; j < fNvar; ++j)
            fBestFitParameterErrors.push_back(-1.); // error undefined

         fOptimizationMethodMode = BCIntegrate::kOptMetropolis;
      }
   }
}

// ---------------------------------------------------------
void BCIntegrate::SetCubaIntegrationMethod(BCIntegrate::BCCubaMethod type)
{
#ifdef HAVE_CUBA_H
   switch(type) {
      case BCIntegrate::kCubaVegas:
      case BCIntegrate::kCubaSuave:
         fCubaIntegrationMethod = type;
         return;
      case BCIntegrate::kCubaDivonne:
      case BCIntegrate::kCubaCuhre:
         BCLog::OutError(TString::Format(
            "BAT does not yet support global setting of Cuba integration method to %s. "
            "To use this method use explicit call to CubaIntegrate() with arguments.",
            DumpCubaIntegrationMethod(type).c_str()));
         return;
      default:
         BCLog::OutError(TString::Format("Integration method of type %d is not defined for Cuba",type));
         return;
   }
#else
   BCLog::OutWarning("!!! This version of BAT is compiled without Cuba.");
   BCLog::OutWarning("    Setting Cuba integration method will have no effect.");
#endif
}

// ---------------------------------------------------------
int BCIntegrate::CubaIntegrand(const int *ndim, const double xx[],
                                              const int * /*ncomp*/, double ff[], void * /*userdata*/)
{
#ifdef HAVE_CUBA_H
   // scale variables
   double jacobian = 1.0;

   std::vector<double> scaled_parameters;

   for (int i = 0; i < *ndim; i++) {
      double range = global_this->fx->at(i)->GetUpperLimit() - global_this->fx->at(i)->GetLowerLimit();

      // multiply range to jacobian
      jacobian *= range;

      // get the scaled parameter value
      scaled_parameters.push_back(global_this->fx->at(i)->GetLowerLimit() + xx[i] * range);
   }

   // call function to integrate
   ff[0] = global_this->Eval(scaled_parameters);

   // multiply jacobian
   ff[0] *= jacobian;

   // multiply fudge factor
   ff[0] *= 1e99;

   // remove parameter vector
   scaled_parameters.clear();
#else
   BCLog::OutError("!!! This version of BAT is compiled without Cuba.");
   BCLog::OutError("    Use other integration methods or install Cuba and recompile BAT.");

   return 1;
#endif

   return 0;
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateCuba(const std::vector<double> &x) {
	return IntegrateCuba(fCubaIntegrationMethod,x);
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateCuba(BCCubaMethod cubatype, const std::vector<double> &x) {
	 
#ifdef HAVE_CUBA_H
	 LogOutputAtStartOfIntegration(kIntCuba, cubatype, x);
	 
   std::vector<double> parameters_double;
   std::vector<double>    parameters_int;

   parameters_double.push_back(fRelativePrecision);
   parameters_double.push_back(fAbsolutePrecision);

   parameters_int.push_back(fCubaVerbosity);
   parameters_int.push_back(fNIterationsMin);
   parameters_int.push_back(fNIterationsMax);

   switch (cubatype) {
      case BCIntegrate::kCubaSuave:
         parameters_int.push_back(fCubaSuaveNNew);
         parameters_double.push_back(fCubaSuaveFlatness);
         break;
      case BCIntegrate::kCubaDivonne:
         break;
      case BCIntegrate::kCubaCuhre:
         break;
      default: // if unknown method run Vegas. Shouldn't ever happen anyway
      case BCIntegrate::kCubaVegas:
         parameters_int.push_back(fCubaVegasNStart);
         parameters_int.push_back(fCubaVegasNIncrease);
   }

   return IntegrateCuba(cubatype, parameters_double, parameters_int,x);
#else
   BCLog::OutError("!!! This version of BAT is compiled without Cuba.");
   BCLog::OutError("    Use other integration methods or install Cuba and recompile BAT.");
   return 0.;
#endif
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateCuba(BCIntegrate::BCCubaMethod method, std::vector<double> parameters_double, std::vector<double> parameters_int, const std::vector<double> &/*x*/)
{
#ifdef HAVE_CUBA_H
   const int NDIM      = int(fx ->size());
   const int NCOMP     = 1;
   const int USERDATA  = 0;
   const int SEED      = 0;
   const int NBATCH    = 1000;
   const int GRIDNO    = -1;
   const char*STATEFILE = "";

   const double EPSREL = parameters_double[0];
   const double EPSABS = parameters_double[1];
   const int VERBOSE   = int(parameters_int[0]);
   const int MINEVAL   = int(parameters_int[1]);
   const int MAXEVAL   = int(parameters_int[2]);

   int fail;
   int nregions;
   double integral[NCOMP];
   double error[NCOMP];
   double prob[NCOMP];

   global_this = this;

   integrand_t an_integrand = &BCIntegrate::CubaIntegrand;

   // reset number of iterations
   fNIterations = 0;

   if (method == 0) {
      // set VEGAS specific parameters
      const int NSTART    = int(parameters_int[3]);
      const int NINCREASE = int(parameters_int[4]);

      // call VEGAS integration method
			Vegas(NDIM, NCOMP, an_integrand, USERDATA,
						EPSREL, EPSABS, VERBOSE, SEED,
						MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
						GRIDNO, STATEFILE,
						&fNIterations, &fail, integral, error, prob);

      // interface for Cuba version 1.5
      /*
      Vegas(NDIM, NCOMP, an_integrand,
            EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
            NSTART, NINCREASE,
            &fNIterations, &fail, integral, error, prob);
      */
   }
   else if (method == 1) {
		 // set SUAVE specific parameters
		 //      const int LAST     = int(parameters_int[3]);
      const int NNEW        = int(parameters_int[3]);
      const double FLATNESS = parameters_double[2];

      // call SUAVE integration method
      Suave(NDIM, NCOMP, an_integrand,
                  USERDATA, EPSREL, EPSABS, VERBOSE, SEED,
                  MINEVAL, MAXEVAL,
                  NNEW, FLATNESS,
						      &nregions, &fNIterations, &fail, integral, error, prob);

      // interface for Cuba version 1.5
      /*
      Suave(NDIM, NCOMP, an_integrand,
         EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
         NNEW, FLATNESS,
         &nregions, &fNIterations, &fail, integral, error, prob);
			*/
   }
   else if (method == 2) {
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
      const int FLAGS    = 0;

			// call DIVONNE integration method
			Divonne(NDIM, NCOMP, an_integrand, USERDATA,
							EPSREL, EPSABS, FLAGS, SEED, MINEVAL, MAXEVAL,
							KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
							NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
							&nregions, &fNIterations, &fail, integral, error, prob);

      // interface for Cuba version 1.5
         /*
            Divonne(NDIM, NCOMP, an_integrand,
            EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
            KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
            NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
            &nregions, &fNIterations, &fail, integral, error, prob);
         */
   }
   else if (method == 3) {
      // set CUHRE specific parameters
      //const int LAST = int(parameters_int[3]);
      const int KEY  = int(parameters_int[4]);
      const int FLAGS    = 0;

      // call CUHRE integration method
			Cuhre(NDIM, NCOMP, an_integrand, USERDATA,
						EPSREL, EPSABS, FLAGS, MINEVAL, MAXEVAL, KEY,
						&nregions, &fNIterations, &fail, integral, error, prob);
			
      // interface for Cuba version 1.5
      /*
      Cuhre(NDIM, NCOMP, an_integrand,
            EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL, KEY,
            &nregions, &fNIterations, &fail, integral, error, prob);
      */
   }
   else {
		 BCLog::OutError(" Integration method not available. ");
		 integral[0] = -1e99;
   }

   if (fail != 0) {
		 BCLog::OutWarning(" Warning, integral did not converge with the given set of parameters. ");
		 BCLog::OutWarning(TString::Format(" neval    = %d", fNIterations));
		 BCLog::OutWarning(TString::Format(" fail     = %d", fail));
		 BCLog::OutWarning(TString::Format(" integral = %e", integral[0]/1e99));
		 BCLog::OutWarning(TString::Format(" error    = %e", error[0]/1e99));
		 BCLog::OutWarning(TString::Format(" prob     = %e", prob[0]));
   }

	 fError = error[0] / 1e99;
   return integral[0] / 1e99;
#else
   BCLog::OutError("!!! This version of BAT is compiled without Cuba.");
   BCLog::OutError("    Use other integration methods or install Cuba and recompile BAT.");
   return 0.;
#endif
}

// ---------------------------------------------------------
void BCIntegrate::MCMCIterationInterface()
{
   // what's within this method will be executed
   // for every iteration of the MCMC

   // fill error band
   MCMCFillErrorBand();

   // do user defined stuff
   MCMCUserIterationInterface();
}

// ---------------------------------------------------------
void BCIntegrate::MCMCFillErrorBand()
{
   if (!fFillErrorBand)
      return;

   // function fitting
   if (fFitFunctionIndexX < 0)
      return;

   // loop over all possible x values ...
   if (fErrorBandContinuous) {
      double x = 0;
      for (int ix = 0; ix < fErrorBandNbinsX; ix++) {
         // calculate x
         x = fErrorBandXY->GetXaxis()->GetBinCenter(ix + 1);

         // calculate y
         std::vector<double> xvec;
         xvec.push_back(x);

         // loop over all chains
         for (int ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
            // calculate y
            double y = FitFunction(xvec, MCMCGetx(ichain));

            // fill histogram
            fErrorBandXY->Fill(x, y);
         }

         xvec.clear();
      }
   }
   // ... or evaluate at the data point x-values
   else {
      int ndatapoints = int(fErrorBandX.size());
      double x = 0;

      for (int ix = 0; ix < ndatapoints; ++ix) {
         // calculate x
         x = fErrorBandX.at(ix);

         // calculate y
         std::vector<double> xvec;
         xvec.push_back(x);

         // loop over all chains
         for (int ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
            // calculate y
            double y = FitFunction(xvec, MCMCGetx(ichain));

            // fill histogram
            fErrorBandXY->Fill(x, y);
         }

         xvec.clear();
      }
   }
}

// ---------------------------------------------------------
void BCIntegrate::SAInitialize()
{
   fSAx.clear();
   fSAx.assign(fNvar, 0.0);
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpIntegrationMethod(BCIntegrate::BCIntegrationMethod type)
{
   switch(type) {
      case BCIntegrate::kIntMonteCarlo:
         return "Sampled Mean Monte Carlo";
      case BCIntegrate::kIntImportance:
         return "Importance Sampling";
      case BCIntegrate::kIntMetropolis:
         return "Metropolis";
      case BCIntegrate::kIntCuba:
         return "Cuba";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpMarginalizationMethod(BCIntegrate::BCMarginalizationMethod type)
{
   switch(type) {
      case BCIntegrate::kMargMonteCarlo:
         return "Monte Carlo Integration";
      case BCIntegrate::kMargMetropolis:
         return "Metropolis MCMC";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpOptimizationMethod(BCIntegrate::BCOptimizationMethod type)
{
   switch(type) {
      case BCIntegrate::kOptSA:
         return "Simulated Annealing";
      case BCIntegrate::kOptMetropolis:
         return "Metropolis MCMC";
      case BCIntegrate::kOptMinuit:
         return "Minuit";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpCubaIntegrationMethod(BCIntegrate::BCCubaMethod type)
{
   switch(type) {
      case BCIntegrate::kCubaVegas:
         return "Vegas";
      case BCIntegrate::kCubaSuave:
         return "Suave";
      case BCIntegrate::kCubaDivonne:
         return "Divonne";
      case BCIntegrate::kCubaCuhre:
         return "Cuhre";
      default:
         return "Undefined";
   }
}

// ---------------------------------------------------------

