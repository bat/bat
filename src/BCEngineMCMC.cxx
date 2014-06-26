/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCEngineMCMC.h"

#include "BCH1D.h"
#include "BCH2D.h"
#include "BCLog.h"
#include "BCMath.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TF1.h>

#include <math.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const char * name)
	: fMCMCFlagWriteChainToFile(false)
	, fMCMCFlagWritePreRunToFile(false)
	, fMCMCFlagWritePreRunObservablesToFile(true)
	, fMCMCScaleFactorLowerLimit(0)
	, fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max())
	, fMCMCEfficiencyMin(0.15)
	, fMCMCEfficiencyMax(0.50)
	, fMCMCFlagInitialPosition(1)
	, fMCMCFlagOrderParameters(true)
	, fMCMCRValueUseStrict(false)
	, fMCMCRValueCriterion(0.1)
	, fMCMCRValueParametersCriterion(0.1)
	, fRandom(new TRandom3())
{
	SetName(name);
	MCMCSetPrecision(BCEngineMCMC::kMedium);
	MCMCSetRandomSeed(0);
}

// ---------------------------------------------------------
void BCEngineMCMC::SetName(const char * name)
{
	fName = name;
	fSafeName = name;
	fSafeName.erase(std::remove_if(fSafeName.begin(),fSafeName.end(),::isspace),fSafeName.end());
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetPrecision(BCEngineMCMC::Precision precision) {
	switch(precision) {

	case BCEngineMCMC::kLow:
		fMCMCNChains                    = 1;
		fMCMCNLag                       = 1;
		fMCMCNIterationsMax             = 10000;
		fMCMCNIterationsRun             = 10000;
		fMCMCNIterationsPreRunMin       = 1500;
		fMCMCNIterationsEfficiencyCheck = 500;
		fMCMCNIterationsUpdate          = 1000;
		fMCMCNIterationsUpdateClear     = 5000;
		break;
		
	case BCEngineMCMC::kQuick:
		fMCMCNChains                    = 2;
		fMCMCNLag                       = 1;
		fMCMCNIterationsMax             = 10000;
		fMCMCNIterationsRun             = 10000;
		fMCMCNIterationsPreRunMin       = 1500;
		fMCMCNIterationsEfficiencyCheck = 500;
		fMCMCNIterationsUpdate          = 1000;
		fMCMCNIterationsUpdateClear     = 5000;
		break;

	case  BCEngineMCMC::kMedium:
		fMCMCNChains                    = 5;
		fMCMCNLag                       = 1;
		fMCMCNIterationsMax             = 100000;
		fMCMCNIterationsRun             = 100000;
		fMCMCNIterationsPreRunMin       = 1500;
		fMCMCNIterationsEfficiencyCheck = 500;
		fMCMCNIterationsUpdate          = 1000;
		fMCMCNIterationsUpdateClear     = 5000;
		break;

	case  BCEngineMCMC::kHigh:
		fMCMCNChains                    = 10;
		fMCMCNLag                       = 10;
		fMCMCNIterationsMax             = 1000000;
		fMCMCNIterationsRun             = 1000000;
		fMCMCNIterationsPreRunMin       = 5000;
		fMCMCNIterationsEfficiencyCheck = 1000;
		fMCMCNIterationsUpdate          = 1000;
		fMCMCNIterationsUpdateClear     = 5000;
		break;

	case  BCEngineMCMC::kVeryHigh:
		fMCMCNChains                    = 10;
		fMCMCNLag                       = 10;
		fMCMCNIterationsMax             = 10000000;
		fMCMCNIterationsRun             = 10000000;
		fMCMCNIterationsPreRunMin       = 10000;
		fMCMCNIterationsEfficiencyCheck = 1000;
		fMCMCNIterationsUpdate          = 1000;
		fMCMCNIterationsUpdateClear     = 5000;
		break;
	}

	// re-initialize
	MCMCInitialize();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetPrecision(const BCEngineMCMC * other) {
	fMCMCNChains                    = other -> MCMCGetNChains();
	fMCMCNLag                       = other -> MCMCGetNLag();
	fMCMCNIterationsMax             = other -> MCMCGetNIterationsMax();
	fMCMCNIterationsRun             = other -> MCMCGetNIterationsRun();
	fMCMCNIterationsPreRunMin       = other -> MCMCGetNIterationsPreRunMin();
	fMCMCNIterationsEfficiencyCheck = other -> MCMCGetNIterationsEfficiencyCheck();
	fMCMCNIterationsUpdate          = other -> MCMCGetNIterationsUpdate();
	fMCMCNIterationsUpdateClear     = other -> MCMCGetNIterationsUpdateClear();
	fMCMCRValueCriterion            = other -> MCMCGetRValueCriterion();
	fMCMCRValueParametersCriterion  = other -> MCMCGetRValueParametersCriterion();
	fMCMCEfficiencyMin              = other -> MCMCGetMinimumEfficiency();
	fMCMCEfficiencyMax              = other -> MCMCGetMaximumEfficiency();

	// re-initialize
	MCMCInitialize();
}

// ---------------------------------------------------------
BCEngineMCMC::~BCEngineMCMC()
{
   // delete random number generator
   delete fRandom;

   // delete 1-d marginalized distributions
   for (unsigned i = 0; i < fH1Marginalized.size(); ++i) {
		 delete fH1Marginalized[i];
	 }
   fH1Marginalized.clear();

   // delete 2-d marginalized distributions
   for (unsigned i = 0; i < fH2Marginalized.size(); ++i) {
		 for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
			 delete fH2Marginalized[i][j];
		 fH2Marginalized[i].clear();
	 }
   fH2Marginalized.clear();
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC & other)
{
	Copy(other);
}

// ---------------------------------------------------------
void BCEngineMCMC::Copy(const BCEngineMCMC & other)
{
   fMCMCPointerToGetProposalPoint            = other.fMCMCPointerToGetProposalPoint;
   fMCMCNChains                              = other.fMCMCNChains;
   fMCMCNLag                                 = other.fMCMCNLag;
   fMCMCNIterations                          = other.fMCMCNIterations;
   fMCMCCurrentIteration                     = other.fMCMCCurrentIteration;
   fMCMCCurrentChain                         = other.fMCMCCurrentChain;
	 fMCMCNIterationsEfficiencyCheck           = other.fMCMCNIterationsEfficiencyCheck;
   fMCMCNIterationsUpdate                    = other.fMCMCNIterationsUpdate;
   fMCMCNIterationsUpdateClear               = other.fMCMCNIterationsUpdateClear;
   fMCMCNIterationsConvergenceGlobal         = other.fMCMCNIterationsConvergenceGlobal;
   fMCMCNIterationsMax                       = other.fMCMCNIterationsMax;
   fMCMCNIterationsRun                       = other.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin                 = other.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue                          = other.fMCMCNTrialsTrue;
   fMCMCNTrials                              = other.fMCMCNTrials;
   fMCMCFlagWriteChainToFile                 = other.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile                = other.fMCMCFlagWritePreRunToFile;
	 fMCMCFlagWritePreRunObservablesToFile = other.fMCMCFlagWritePreRunObservablesToFile;
   fMCMCTrialFunctionScaleFactor             = other.fMCMCTrialFunctionScaleFactor;
   fMCMCTrialFunctionScaleFactorStart        = other.fMCMCTrialFunctionScaleFactorStart;
   fMCMCFlagPreRun                           = other.fMCMCFlagPreRun;
   fMCMCFlagRun                              = other.fMCMCFlagRun;
   fMCMCInitialPosition                      = other.fMCMCInitialPosition;
   fMCMCEfficiencies                         = other.fMCMCEfficiencies;
   fMCMCEfficiencyMin                        = other.fMCMCEfficiencyMin;
   fMCMCEfficiencyMax                        = other.fMCMCEfficiencyMax;
	 fMCMCScaleFactorLowerLimit                = other.fMCMCScaleFactorLowerLimit;
	 fMCMCScaleFactorUpperLimit                = other.fMCMCScaleFactorUpperLimit;
   fMCMCFlagInitialPosition                  = other.fMCMCFlagInitialPosition;
   fMCMCFlagOrderParameters                  = other.fMCMCFlagOrderParameters;
   fMCMCPhase                                = other.fMCMCPhase;
   fMCMCx                                    = other.fMCMCx;
   fMCMCObservables                      = other.fMCMCObservables;
   fMCMCxMax                                 = other.fMCMCxMax;
   fMCMCxMean                                = other.fMCMCxMean;
   fMCMCxVar                                 = other.fMCMCxVar;
   fMCMCprob                                 = other.fMCMCprob;
   fMCMCprobMax                              = other.fMCMCprobMax;
   fMCMCprobMean                             = other.fMCMCprobMean;
   fMCMCprobVar                              = other.fMCMCprobVar;
   fMCMCRValueUseStrict                      = other.fMCMCRValueUseStrict;
   fMCMCRValueCriterion                      = other.fMCMCRValueCriterion ;
   fMCMCRValueParametersCriterion            = other.fMCMCRValueParametersCriterion;
   fMCMCRValue                               = other.fMCMCRValue;
   fMCMCRValueParameters                     = other.fMCMCRValueParameters;
	 fRandom                                   = (other.fRandom) ? new TRandom3(*other.fRandom) : NULL;
   fMCMCThreadLocalStorage                   = other.fMCMCThreadLocalStorage;

	 fH1Marginalized = std::vector<TH1D*>(other.fH1Marginalized.size(),0);
   for (unsigned i = 0; i < other.fH1Marginalized.size(); ++i)
		 if (other.fH1Marginalized.at(i))
			 fH1Marginalized[i] = new TH1D(*(other.fH1Marginalized[i]));

	 fH2Marginalized = std::vector<std::vector<TH2D*> > (other.fH2Marginalized.size(), std::vector<TH2D*>(other.fH2Marginalized[0].size(),0));
   for (unsigned i = 0; i < other.fH2Marginalized.size(); ++i)
		 for (unsigned j = 0; j < other.fH2Marginalized[i].size(); ++j)
			 if (other.fH2Marginalized[i][j])
				 fH2Marginalized[i][j] = new TH2D(*(other.fH2Marginalized[i][j]));

	 fMCMCTrees.assign( other.fMCMCTrees.size() , 0 );

   fMarginalModes         = other.fMarginalModes;
   fMCMCBestFitParameters = other.fMCMCBestFitParameters;
   fMCMCLogMaximum        = other.fMCMCLogMaximum;
}

// ---------------------------------------------------------
BCEngineMCMC & BCEngineMCMC::operator = (const BCEngineMCMC & enginemcmc)
{
   Copy(enginemcmc);
   return *this;
}

// --------------------------------------------------------
unsigned int BCEngineMCMC::GetNFreeParameters()
{
  return (GetNParameters() - GetNFixedParameters());
}

// --------------------------------------------------------
unsigned int BCEngineMCMC::GetNFixedParameters()
{
   int n = 0;
   for (unsigned int i = 0; i < fParameters.Size(); ++i) {
      if (fParameters[i]->Fixed())
         ++n;
   }

   return n;
}

// --------------------------------------------------------
void BCEngineMCMC::SetNbins(unsigned int nbins)
{
	fParameters.SetNBins(nbins);
	fObservables.SetNBins(nbins);
}

// --------------------------------------------------------
TH1D * BCEngineMCMC::GetMarginalizedHistogram(unsigned index) const {
	if ( index >= fH1Marginalized.size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.",index));
		return 0;
	}

	if (fH1Marginalized[index])
		return fH1Marginalized[index];
	
	// else output warning
	if ( index<fParameters.Size() ) // Marginalization of model parameter
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for parameter %s", fParameters[index]->GetName().data()));
	else if ( index-fParameters.Size()<fObservables.Size() ) // Marginalization of user-defined observable
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for user-observable %s", fObservables[index-fParameters.Size()]->GetName().data()));
	return 0;
}

// --------------------------------------------------------
TH2D * BCEngineMCMC::GetMarginalizedHistogram(unsigned i, unsigned j) const {
	if (i == j) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Called with identical indices %u.", i));
		return 0;
	}
	
	// redundancy reduced by only storing H2's for i<j
	if (i > j)
		return GetMarginalizedHistogram(j,i);
	// TODO: transpose histogram;

	if ( i >= fH2Marginalized.size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", i));
		return 0;
	}
	if ( j >= fH2Marginalized[i].size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", j));
		return 0;
	}

	if (fH2Marginalized[i][j])
		return fH2Marginalized[i][j];

	// else output warning
	if (!fH1Marginalized[i]) {
		if ( i<fParameters.Size() ) // Marginalization of model parameter
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for parameter %s", fParameters[i]->GetName().data()));
		else if ( i-fParameters.Size()<fObservables.Size() ) // Marginalization of user-defined observable
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for user-observable %s", fObservables[i-fParameters.Size()]->GetName().data()));
	}
	if (!fH1Marginalized[j]) {
		if ( j<fParameters.Size() ) // Marginalization of model parameter
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for parameter %s", fParameters[j]->GetName().data()));
		else if ( j-fParameters.Size()<fObservables.Size() ) // Marginalization of user-defined observable
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for user-observable %s", fObservables[j-fParameters.Size()]->GetName().data()));
	}
	return 0;
}

// --------------------------------------------------------
BCH1D * BCEngineMCMC::GetMarginalized(unsigned index) {
	TH1D * h = GetMarginalizedHistogram(index);

	if ( !h )
		return 0;
	
	BCH1D * hprob = new BCH1D(h);
	
	// if best fit parameters exists, set global mode
	if (GetBestFitParameters().size() == fParameters.Size()) {
		double mode = (index<fParameters.Size()) ? GetBestFitParameter(index) : fObservables[index-fParameters.Size()]->Evaluate(GetBestFitParameters());
		hprob -> SetGlobalMode(mode);
	}
	
	if ( fMarginalModes.size() != fH1Marginalized.size() )
		fMarginalModes.resize(fH1Marginalized.size(), 0);
	fMarginalModes[index] = hprob->GetMode();
	
	return hprob;
}
	
// --------------------------------------------------------
BCH2D * BCEngineMCMC::GetMarginalized(unsigned i, unsigned j) {
	TH2D * h = GetMarginalizedHistogram(i,j);

	if (!h)
		return 0;

	BCH2D * hprob = new BCH2D(h);

	// if best fit parameters exists, set global mode
	if (GetBestFitParameters().size() == fParameters.Size()) {
		double imode = (i<fParameters.Size()) ? GetBestFitParameter(i) : fObservables[i-fParameters.Size()]->Evaluate(GetBestFitParameters());
		double jmode = (j<fParameters.Size()) ? GetBestFitParameter(j) : fObservables[j-fParameters.Size()]->Evaluate(GetBestFitParameters());
		hprob -> SetGlobalMode(imode,jmode);
	}
		
	return hprob;
}

// ---------------------------------------------------------
BCVariable * BCEngineMCMC::GetVariable(unsigned index) const {
	if (index<fParameters.Size())
		return fParameters[index];
	if (index-fParameters.Size()<fObservables.Size())
		return fObservables[index-fParameters.Size()];
	return 0;
}

// ---------------------------------------------------------
const std::vector<double> & BCEngineMCMC::GetBestFitParametersMarginalized() const
{
   if(fMarginalModes.empty())
      BCLog::OutError("BCIntegrate::GetBestFitParameterMarginalized : MCMC not yet run, returning center of the range.");

   return fMarginalModes;
}


// ---------------------------------------------------------
double BCEngineMCMC::GetBestFitParameter(unsigned index) const
{

	if (index > fParameters.Size() + fObservables.Size())
		return -std::numeric_limits<double>::infinity();

	std::vector<double> pars = GetBestFitParameters();

	if (pars.size() != fParameters.Size()) {
		BCLog::OutError("BCEngineMCMC::GetBestFitParameter : MCMC not yet run, returning center of the range.");
		return fParameters[index]->GetRangeCenter();
	}

	// parameter
	if ( index < fParameters.Size() )
		return pars[index];
	
	// user-define observable
	return fObservables[index-fParameters.Size()]->Evaluate(pars);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetMaximumPoint(unsigned i) const {
	if (i < fMCMCxMax.size() )
		return fMCMCxMax[i];
	
	std::vector<double> x;
	return x;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetNChains(unsigned n)
{
   fMCMCNChains = n;

   MCMCInitialize();						// gre-initialize
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(const std::vector<std::vector<double> > & x0s) {
	fMCMCInitialPosition = x0s;
	MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(const std::vector<double> & x0s) {
   // clear the existing initial position vector
   fMCMCInitialPosition.clear();

	 std::vector<double>::const_iterator it = x0s.begin();
	 while (it+fParameters.Size()<=x0s.end()) {
		 fMCMCInitialPosition.push_back(std::vector<double>(it,it+fParameters.Size()));
		 it += fParameters.Size();
	 }																		
	 MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(bool flag)
{
	fParameters.FillHistograms(flag);
	fObservables.FillHistograms(flag);	
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetMarkovChainTrees(const std::vector<TTree *> & trees)
{
	// clear vector
	fMCMCTrees.clear();

	// copy tree
	for (unsigned i = 0; i < trees.size(); ++i)
		fMCMCTrees.push_back(trees[i]);
}

void BCEngineMCMC::MCMCSetRandomSeed(unsigned seed)
{
	// set main generator
	fRandom->SetSeed(seed);

	// call once so return value of GetSeed() fixed
	fRandom->Rndm();

	SyncThreadStorage();

	// type conversion to avoid compiler warnings
	if (size_t(fMCMCNChains) != fMCMCThreadLocalStorage.size())
		BCLog::OutError(Form("#chains does not match #(thread local storages): %d vs %u",
												 fMCMCNChains, unsigned(fMCMCThreadLocalStorage.size())));

	// set all single chain generators
	for (unsigned i = 0; i < fMCMCNChains ; ++i){
		// call once so return value of GetSeed() fixed
		fMCMCThreadLocalStorage[i].rng->SetSeed(fRandom->GetSeed() + i);
		fMCMCThreadLocalStorage[i].rng->Rndm();
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInitializeMarkovChainTrees()
{
	// clear vector
	fMCMCTrees.clear();

	// create new trees
	for (unsigned i = 0; i < fMCMCNChains; ++i) {
		fMCMCTrees.push_back(new TTree(TString::Format("MarkovChainTree_%i", i), "MarkovChainTree"));
		fMCMCTrees[i]->Branch("Iteration",       &fMCMCNIterations[i],  "iteration/i");
		// todo check example and parallel_TEST how to automatically determine #parameters
		//      fMCMCTrees[i]->Branch("NParameters",     fParameters.Size(),   "parameters/I");
		fMCMCTrees[i]->Branch("LogProbability",  &fMCMCprob[i],         "log(probability)/D");
		fMCMCTrees[i]->Branch("Phase",           &fMCMCPhase,           "phase/I");

		for (unsigned j = 0; j < fParameters.Size(); ++j) {
			fMCMCTrees[i]->Branch(TString::Format("Parameter%i", j), &fMCMCx[i][j], TString::Format("parameter %i/D", j));
			fMCMCTrees[i]->SetAlias(fParameters[j]->GetName().data(),TString::Format("Parameter%i", j));
		}
		for (unsigned j = 0; j < fObservables.Size(); ++j) {
			fMCMCTrees[i]->Branch(TString::Format("Observable%i", j), &fMCMCObservables[i][j], TString::Format("observable %i/D", j));
			fMCMCTrees[i]->SetAlias(fObservables[j]->GetName().data(),TString::Format("Observable%i", j));
		}		
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunction(unsigned ichain, std::vector<double> &x)
{
	// call MCMCTrialFunctionSingle() for all parameters by default
	for (unsigned i = 0; i < fParameters.Size(); ++i)
		x[i] = MCMCTrialFunctionSingle(ichain, i);
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCTrialFunctionSingle(unsigned ichain, unsigned iparameter) {
	// no check of range for performance reasons

	// use uniform distribution
	//   return = fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter] * 2.0 * (0.5 - fRandom->Rndm());

	// Breit-Wigner width adjustable width
	return fMCMCThreadLocalStorage[ichain].rng->BreitWigner(0.0, fMCMCTrialFunctionScaleFactor[ichain][iparameter]);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(unsigned ichain) const {
	if ( ichain < fMCMCTrialFunctionScaleFactor.size() )
		return fMCMCTrialFunctionScaleFactor[ichain];

	std::vector<double> x;
	return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(unsigned ichain, unsigned ipar) {
	// check if ichain is in range
	if (ichain >= fMCMCTrialFunctionScaleFactor.size() or ipar >= fMCMCTrialFunctionScaleFactor[ichain].size() )
		return 0;
	return fMCMCTrialFunctionScaleFactor[ichain][ipar];
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetx(unsigned ichain) {
	if ( ichain < fMCMCx.size() )
		return fMCMCx[ichain];

	std::vector<double> x;
	return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetx(unsigned ichain, unsigned ipar) const
{
	if ( ichain >= fMCMCx.size() or ipar >= fMCMCx[ipar].size() )
		return 0;
	return fMCMCx[ichain][ipar];
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetLogProbx(unsigned ichain)
{
	// check if ichain is in range
	if (ichain >= fMCMCNChains)
		return -std::numeric_limits<double>::max();

	// return log of the probability at the current point in the ith chain
	return fMCMCprob.at(ichain);
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned chain, std::vector<double> &x)
{
	// get unscaled random point. this point might not be in the correct volume.
	MCMCTrialFunction(chain, x);

	// get a proposal point from the trial function and scale it
	for (unsigned i = 0; i < fParameters.Size(); ++i)
		if (fParameters[i]->Fixed())
			x[i] = fParameters[i] -> GetFixedValue();
		else
			x[i] = fMCMCx[chain][i] + x[i]*(fParameters[i]->GetUpperLimit()-fParameters[i]->GetLowerLimit());

	// check if the point is in the correct volume.
	for (unsigned i = 0; i < fParameters.Size(); ++i)
		if (!fParameters[i]->IsValid(x[i]))
			return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned ichain, unsigned ipar, std::vector<double> &x)
{
  // copy the old point into the new
	x = fMCMCx[ichain];
	
  // check if parameter is fixed
  if (fParameters[ipar]->Fixed()) {
    x[ipar] = fParameters[ipar]->GetFixedValue();
    return true; // assume that value is inside allowed region
  }

  // get unscaled random point in the dimension of the chosen
  // parameter. this point might not be in the correct volume.
  double proposal = MCMCTrialFunctionSingle(ichain, ipar);

  // modify the parameter under study
  x[ipar] += proposal * (fParameters[ipar]->GetUpperLimit() - fParameters[ipar]->GetLowerLimit());

  // check if the point is in the correct volume.
	return fParameters[ipar] -> IsValid(x[ipar]);
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(unsigned chain, unsigned parameter) {
	fMCMCCurrentChain = chain;
	
	// increase counter
	fMCMCNIterations[chain]++;

	// get proposal point
	if ( MCMCGetProposalPointMetropolis(chain, parameter, fMCMCThreadLocalStorage[chain].xLocal) ) {
		// calculate probabilities of the old and new points
		double p0 = (std::isfinite(fMCMCprob[chain])) ? fMCMCprob[chain] : -std::numeric_limits<double>::max();
		double p1 = LogEval(fMCMCThreadLocalStorage[chain].xLocal);

		if (std::isfinite(p1)) {
			// if the new point is more probable, keep it; or else throw dice
			if ( p1 >= p0 or log(fMCMCThreadLocalStorage[chain].rng->Rndm()) < (p1 - p0) ) {
					// increase counter
				fMCMCNTrialsTrue[chain][parameter]++;
				// copy the point
				fMCMCx[chain][parameter] = fMCMCThreadLocalStorage[chain].xLocal[parameter];
				// save the probability of the point
				fMCMCprob[chain] = p1;
				 
				// execute user code
				MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, true);
				return true;
			}
		} else {						// new log(likelihood) was not a finite number
			BCLog::OutDebug(Form("Log(likelihood) evaluated to nan or inf in chain %i while varying parameter %s to %.3e",chain,fParameters[parameter]->GetName().data(),fMCMCThreadLocalStorage[chain].xLocal[parameter]));
			// print parameter point
		}
	}
	
	// execute user code
	MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, false);
	return false;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(unsigned chain) {
	fMCMCCurrentChain = chain;

	// increase counter
	fMCMCNIterations[chain]++;

	// get proposal point
	if ( MCMCGetProposalPointMetropolis(chain, fMCMCThreadLocalStorage[chain].xLocal) ) {
		// calculate probabilities of the old and new points
		double p0 = (std::isfinite(fMCMCprob[chain])) ? fMCMCprob[chain] : -std::numeric_limits<double>::max();
		double p1 = LogEval(fMCMCThreadLocalStorage[chain].xLocal);

		if (std::isfinite(p1)) {
			// if the new point is more probable, keep it; or else throw dice
			if ( p1 >= p0 or log(fMCMCThreadLocalStorage[chain].rng->Rndm()) < (p1-p0) ) {
				// increase counters
				for (unsigned i = 0; i < fParameters.Size(); ++i)
					fMCMCNTrialsTrue[chain][i]++;
				// copy the point
				fMCMCx[chain] = fMCMCThreadLocalStorage[chain].xLocal;
				// save the probability of the point
				fMCMCprob[chain] = p1;

				// execute user code
				MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, true);
				return true;
			}
		} else { // new log(likelihood) was not a finite number
			BCLog::OutDebug("Log(likelihood) evaluated to nan or inf at");
			// TODO print parameter point
		}
	}

	// execute user code for every point
	MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, false);
	return false;
}

//--------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis() {
	bool return_value = true;

	if ( fMCMCFlagOrderParameters ) { // run over pars one at a time

		for (unsigned ipar = 0; ipar < fParameters.Size(); ++ipar)

			if ( !fParameters[ipar]->Fixed() ) {

				//loop over chains
				unsigned chunk = 1; (void) chunk;
				unsigned ichain;    (void) ichain;
#pragma omp parallel for shared(chunk) private(ichain) schedule(static, chunk)
				for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain)
					return_value *= MCMCGetNewPointMetropolis(ichain,ipar);
				
				fMCMCCurrentChain = -1;
				MCMCInChainCheckMaximum();

			}
	}

	else {											// run over all pars at once

		//loop over chains
		unsigned chunk = 1; (void) chunk;
		unsigned ichain;    (void) ichain;
#pragma omp parallel for shared(chunk) private(ichain) schedule(static, chunk)
		for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain)
			return_value *= MCMCGetNewPointMetropolis(ichain);
		
		fMCMCCurrentChain = -1;
		MCMCInChainCheckMaximum();

	}

	++fMCMCCurrentIteration;
	return return_value;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainCheckMaximum() {
	// loop over all chains
	for (unsigned i = 0; i < fMCMCNChains; ++i)
		// check if new maximum is found or chain is at the beginning
		if (fMCMCprob[i] > fMCMCprobMax[i] || fMCMCNIterations[i] == 1) {
			// copy maximum value
			fMCMCprobMax[i] = fMCMCprob[i];
			// copy mode of chain
			fMCMCxMax[i] = fMCMCx[i];
		}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainUpdateStatistics() {
	for (unsigned int c = 0; c<fMCMCx.size(); ++c)
		for (unsigned int p = 0; p<fMCMCx[c].size(); ++p) {
			fMCMCxMean[c][p] += (fMCMCx[c][p] - fMCMCxMean[c][p]) / fMCMCNTrials;
			
			if ( fMCMCNTrials > 1 ) // calculate variance of each chain for this part
				fMCMCxVar[c][p] = (1 - 1./fMCMCNTrials)*fMCMCxVar[c][p] + (fMCMCx[c][p]-fMCMCxMean[c][p])*(fMCMCx[c][p]-fMCMCxMean[c][p]) / (fMCMCNTrials - 1);
		}

	// loop over chains
	for (unsigned i = 0; i < fMCMCNChains; ++i) {
		// calculate mean value of each chain for this part
		fMCMCprobMean[i] += (fMCMCprob[i] - fMCMCprobMean[i]) / fMCMCNTrials;
		
		if (fMCMCNTrials > 1) // calculate variance of each chain for this part
			fMCMCprobVar[i] = (1 - 1./fMCMCNTrials)*fMCMCprobVar[i]	+ (fMCMCprob[i]-fMCMCprobMean[i])*(fMCMCprob[i]-fMCMCprobMean[i]) / (fMCMCNTrials - 1);
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainFillHistograms() {
	// loop over chains
	for (unsigned i = 0; i < fMCMCNChains; ++i) {

		////////////////////////////////////////
		// fill each 1-dimensional histogram (if supposed to be filled)

		// Parameters
		for (unsigned j = 0; j < fParameters.Size(); ++j)
			if (TH1 * h = fH1Marginalized[j])
				h -> Fill(fMCMCx[i][j]);
		// User-defined Observables
		for (unsigned j = 0; j < fObservables.Size(); ++j)
			if (TH1 * h = fH1Marginalized[j+fParameters.Size()])
				h -> Fill(fMCMCObservables[i][j]);


		////////////////////////////////////////
		// fill each 2-dimensional histogram (if supposed to be filled)

		for (unsigned j = 0; j < fParameters.Size(); ++j) {
			// Parameter vs Parameter
			for (unsigned k = j+1; k < fParameters.Size(); ++k)
				if (TH2D * h = fH2Marginalized[j][k])
					h -> Fill(fMCMCx[i][j],fMCMCx[i][k]);
			// User-defined Observable vs Parameter
			for (unsigned k = 0; k < fObservables.Size(); ++k)
				if (TH2D * h = fH2Marginalized[j][k+fParameters.Size()])
					h -> Fill(fMCMCx[i][j],fMCMCObservables[i][k]);
		}
		// User-defined Observable vs User-defined Observable
		for (unsigned j = 0; j < fObservables.Size(); ++j)
			for (unsigned k = j+1; k < fObservables.Size(); ++k)
				if (TH2D * h = fH2Marginalized[j+fParameters.Size()][k+fParameters.Size()])
					h -> Fill(fMCMCObservables[i][j],fMCMCObservables[i][k]);
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainTestConvergenceAllChains() {
	if (fMCMCNChains <= 1 or fMCMCNTrials <= 1)
		return;

	// define flag for convergence
	bool flag_convergence = true;
	
	// extract means and variances
	std::vector<double> means(fMCMCNChains);
	std::vector<double> variances(fMCMCNChains);

	// loop over parameters
	for (unsigned iparameter = 0; iparameter < fParameters.Size(); ++iparameter) {
		if ( fParameters[iparameter] -> Fixed() )
			continue;

		// loop over chains
		for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain) {
			means[ichain] = fMCMCxMean[ichain][iparameter];
			variances[ichain] = fMCMCxVar[ichain][iparameter];
		}

		fMCMCRValueParameters[iparameter] = BCMath::Rvalue(means, variances, fMCMCNTrials, fMCMCRValueUseStrict);

		// set flag to false if convergence criterion is not fulfilled for the parameter
		flag_convergence *= ((fMCMCRValueParameters[iparameter]-1) < fMCMCRValueParametersCriterion);

	}

	fMCMCRValue = BCMath::Rvalue(fMCMCprobMean, fMCMCprobVar, fMCMCNTrials, fMCMCRValueUseStrict);

	// set flag to false if convergence criterion is not fulfilled for the log-likelihood
	flag_convergence *= ((fMCMCRValue - 1) < fMCMCRValueCriterion);
			
	// remember number of iterations needed to converge
	if (fMCMCNIterationsConvergenceGlobal == -1 && flag_convergence == true)
		fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0] / GetNFreeParameters();
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainWriteChains()
{
	// loop over all chains
	for (unsigned i = 0; i < fMCMCNChains; ++i)
		fMCMCTrees[i]->Fill();
}

// --------------------------------------------------------
double BCEngineMCMC::LogEval(const std::vector<double> & /*parameters*/)
{
	// test function for now
	// this will be overloaded by the user
	return 0.0;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCMetropolisPreRun() {
	// print on screen
	BCLog::OutSummary("Pre-run Metropolis MCMC...");
	
	// initialize Markov chain
	MCMCInitialize();
	MCMCInitializeMarkovChains();
	
	// reset run statistics
	MCMCResetRunStatistics();

	// perform run
	BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsMax));

	// set phase and cycle number
	fMCMCPhase = 1;

	//////////////////////////////////////////////////
	// Adjust scales until all parameters are in correct efficiency range in all chains

	fMCMCCurrentIteration = 0;
	fMCMCNTrials = 0;
	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(fParameters.Size(),0));
	fMCMCEfficiencies.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
	bool allEfficient = false;
	bool inefficientScalesAdjustable = true;
	while (fMCMCCurrentIteration < (int)fMCMCNIterationsMax and !allEfficient and inefficientScalesAdjustable) {

		MCMCGetNewPointMetropolis();
		++fMCMCNTrials;

		if (fMCMCFlagWritePreRunToFile) {
			if (fMCMCFlagWritePreRunObservablesToFile)
				CalculateObservables();
			MCMCInChainWriteChains();
		}

		if ( fMCMCNTrials != fMCMCNIterationsEfficiencyCheck)
			continue;

		allEfficient = true;
		inefficientScalesAdjustable = false;
		for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain) {
			for (unsigned ipar = 0; ipar < fParameters.Size(); ++ipar) {

				if (fParameters[ipar]->Fixed())
					continue;

				fMCMCEfficiencies[ichain][ipar] = 1. * fMCMCNTrialsTrue[ichain][ipar] / fMCMCNIterationsEfficiencyCheck;

				if ( fMCMCEfficiencies[ichain][ipar] < fMCMCEfficiencyMin ) {
					// if efficiency too low ...
					
					if (allEfficient)			// print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					
					allEfficient = false;

					fMCMCTrialFunctionScaleFactor[ichain][ipar] /= (fMCMCEfficiencies[ichain][ipar] < 0.5*fMCMCEfficiencyMin) ? 4 : 2;
					
					if ( fMCMCTrialFunctionScaleFactor[ichain][ipar] > fMCMCScaleFactorLowerLimit ) {
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale decreased to %6.2f %%", fParameters.MaxNameLength(), fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					}	else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorLowerLimit;
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale now at lower limit (%6.2f %%)",	fParameters.MaxNameLength(), fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCScaleFactorLowerLimit));
					}
					
				} else if (fMCMCEfficiencies[ichain][ipar] > fMCMCEfficiencyMax ) {
					// if efficiency too high ...

					if (allEfficient)		 // print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					
					allEfficient = false;

					fMCMCTrialFunctionScaleFactor[ichain][ipar] *= 2;

					if ( fMCMCTrialFunctionScaleFactor[ichain][ipar] < fMCMCScaleFactorUpperLimit ) {
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale increased to %6.2f %%", fParameters.MaxNameLength(), fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					} else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorUpperLimit;
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale now at upper limit (%6.2f %%)", fParameters.MaxNameLength(), fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, fMCMCScaleFactorUpperLimit));														 
					}
				}
			}
		}
		fMCMCNTrials = 0;
		fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(fParameters.Size(),0));
	}
	if (allEfficient)
		BCLog::OutDetail(Form("     * Efficiency status: Efficiencies within predefined range after %i iterations.",fMCMCCurrentIteration));
	else if (!inefficientScalesAdjustable)
		BCLog::OutWarning(Form("     * Efficiency status: Some efficiencies outside predefined range, but scales are at limits after %i iterations.",fMCMCCurrentIteration));
	else
		BCLog::OutDetail(Form("     * Efficiency status: Some efficiencies outside predefined range, but maximum number of iterations (%i) reached.",fMCMCNIterationsMax));


	// continue measuring efficiency
	unsigned NTrialsForEff = fMCMCNTrials;

	if (fMCMCNChains > 1) {
		//////////////////////////////////////////////////
		// Run until all chains have converged

		unsigned nIterationsCheckConvergence = fMCMCNIterationsUpdate;
		if ( nIterationsCheckConvergence > fMCMCNIterationsUpdateClear )
			nIterationsCheckConvergence = fMCMCNIterationsUpdateClear;

		fMCMCNTrials = fMCMCNIterationsUpdateClear;
	
		fMCMCNIterationsConvergenceGlobal = -1;
		
		if (fMCMCCurrentIteration >= (int)fMCMCNIterationsMax) {
      BCLog::OutWarning(" Convergence never checked !");
      BCLog::OutWarning("   Increase maximum number of iterations in the pre-run /MCMCSetNIterationsMax()/");
		}

		while ( fMCMCCurrentIteration < (int)fMCMCNIterationsMax and fMCMCNIterationsConvergenceGlobal < 0 ) {
		
			if ( fMCMCNTrials == fMCMCNIterationsUpdateClear ) {
				fMCMCxMean = fMCMCx;
				fMCMCxVar.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
				fMCMCprobMean = fMCMCprob;
				fMCMCprobVar.assign(fMCMCNChains,0);
				fMCMCNTrials = 0;
			}

			MCMCGetNewPointMetropolis();
			++fMCMCNTrials;
			++NTrialsForEff;

			if (fMCMCFlagWritePreRunToFile) {
				if (fMCMCFlagWritePreRunObservablesToFile)
					CalculateObservables();
				MCMCInChainWriteChains();
			}

			MCMCInChainUpdateStatistics();

			if ( fMCMCNTrials % nIterationsCheckConvergence != 0 and fMCMCCurrentIteration < (int)fMCMCNIterationsMax )
				continue;

			MCMCInChainTestConvergenceAllChains();

			if ( fMCMCNIterationsConvergenceGlobal > 0 )
				continue;

			BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));
			BCLog::OutDetail(Form("       - %-*s : R-Value",fParameters.MaxNameLength(),"Parameter"));
				
			for (unsigned ipar = 0; ipar < fParameters.Size(); ++ipar) {

				if ( fParameters[ipar]->Fixed() )
					continue;

				if( fMCMCRValueParameters[ipar]-1 < fMCMCRValueParametersCriterion )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f",fParameters.MaxNameLength(),fParameters[ipar]->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else if ( fMCMCRValueParameters.at(ipar) != std::numeric_limits<double>::max() )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f <--",fParameters.MaxNameLength(),fParameters[ipar]->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else
					BCLog::OutDetail(TString::Format("         %-*s :  MAX_DOUBLE <--",fParameters.MaxNameLength(),fParameters[ipar]->GetName().data()));
			}

			if( fMCMCRValue-1 < fMCMCRValueCriterion )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f", fMCMCRValue));
			else if ( fMCMCRValue != std::numeric_limits<double>::max() )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f <--", fMCMCRValue));
			else
				BCLog::OutDetail("       - Log-Likelihood :  MAX_DOUBLE <--");
		}

		if ( fMCMCNIterationsConvergenceGlobal > 0 )
			if (allEfficient)
				BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations, and all scales are adjusted.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
			else
				BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations, but could not adjust all scales.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
		else
			if (allEfficient)
				BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, but all scales are adjusted.", fMCMCNChains, fMCMCNIterationsMax));
			else
				BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, and could not adjust all scales.", fMCMCNChains, fMCMCNIterationsMax));
	}

		
	//////////////////////////////////////////////////
	// Run until pre-run iteration min criteria met
	if ( fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin ) {

		unsigned N = fMCMCNIterationsPreRunMin - fMCMCCurrentIteration;
		unsigned nwrite = UpdateFrequency(N);

		BCLog::OutDetail(Form(" Current iteration (%d) is below minimum for pre-run. Running %d more iterations.",fMCMCCurrentIteration,N));

		unsigned n = 0;
		while ( fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin ) {
			MCMCGetNewPointMetropolis();
			++NTrialsForEff;
			++n;

			if (fMCMCFlagWritePreRunToFile) {
				if (fMCMCFlagWritePreRunObservablesToFile)
					CalculateObservables();
				MCMCInChainWriteChains();
			}

			if ( n % nwrite == 0)
				BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, 100.*n/N));
		}
		BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations in pre-run.", fMCMCCurrentIteration));
	}

	// print scale factors and efficiencies
	std::vector<double> scalefactors (fParameters.Size(),0);
	std::vector<double> efficiencies (fParameters.Size(),0);

	BCLog::OutDetail(Form(" --> Average scale factors and efficiencies (measured in %d iterations):",(NTrialsForEff==0) ? fMCMCNIterationsEfficiencyCheck : NTrialsForEff));
	BCLog::OutDetail(Form("       - %-*s : Scale factor    Efficiency",fParameters.MaxNameLength(),"Parameter"));
	for (unsigned i = 0; i < fParameters.Size(); ++i) {
		if (fParameters[i]->Fixed())
			continue;
		for (unsigned j = 0; j < fMCMCNChains; ++j) {
			efficiencies[i] += ( NTrialsForEff==0 ) ? fMCMCEfficiencies[j][i] / fMCMCNChains : 1.*fMCMCNTrialsTrue[j][i]/NTrialsForEff/fMCMCNChains;
			scalefactors[i] += fMCMCTrialFunctionScaleFactor[j][i] / fMCMCNChains;
		}
		BCLog::OutDetail(Form("         %-*s :     %6.02f %%        %4.1f %%",fParameters.MaxNameLength(),fParameters[i]->GetName().data(), 100.*scalefactors[i], 100.*efficiencies[i]));
	}

	// reset current iteration
	fMCMCCurrentIteration = -1;

	// reset current chain
	fMCMCCurrentChain = -1;

	// no error
	return 1;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCMetropolis()
{
  // check the number of free parameters
  if (GetNFreeParameters() <= 0) {
    BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Number of free parameters <= 0. Do not run Metropolis.");
    return 0;
  }

	// check if prerun should be performed
	if (fMCMCFlagPreRun)
		MCMCMetropolisPreRun();
	else
		BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Not running prerun. This can cause trouble if the data have changed.");

	// print to screen
	BCLog::OutSummary( "Run Metropolis MCMC ...");

	// reset run statistics
	MCMCResetRunStatistics();

	// set phase and cycle number
	fMCMCPhase = 2;

	BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

	unsigned nwrite = UpdateFrequency(fMCMCNIterationsRun);

	// start the run
	fMCMCCurrentIteration = 0;
	while ( fMCMCCurrentIteration < (int)fMCMCNIterationsRun ) {

		MCMCGetNewPointMetropolis();

		if ( fMCMCCurrentIteration % nwrite == 0 )
			BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration)/(double)fMCMCNIterationsRun*100.));
		
		if (fMCMCCurrentIteration % fMCMCNLag != 0) // apply lag
			continue;
			
		MCMCIterationInterface();		// user action (overloadable)
		
		CalculateObservables();

		// fill histograms
		if ( !fH1Marginalized.empty() or !fH2Marginalized.empty() )
			MCMCInChainFillHistograms();
		
		// write chain to file
		if ( fMCMCFlagWriteChainToFile )
			MCMCInChainWriteChains();
		
	} // end run

	BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));


	// find global maximum and calculate efficiencies
	unsigned probmaxindex = 0;

	fMCMCEfficiencies.assign(fMCMCNTrialsTrue.size(), std::vector<double>(GetNParameters(),0));
	std::vector<double> efficiencies (GetNParameters(),0);

	// loop over all chains
	for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain) {
		// find maximum point
		if (fMCMCprobMax[ichain] > fMCMCprobMax[probmaxindex])
			probmaxindex = ichain;
		// calculate efficiencies
		for (unsigned ipar = 0; ipar<GetNParameters(); ++ipar) {
			fMCMCEfficiencies[ichain][ipar] = 1.*fMCMCNTrialsTrue[ichain][ipar]/fMCMCNIterationsRun;
			efficiencies[ipar] += fMCMCEfficiencies[ichain][ipar]/fMCMCNChains;
		}
	}

	// print efficiencies
	BCLog::OutDetail(Form(" --> Average efficiencies (measured in %d iterations):",fMCMCNIterationsRun));
	BCLog::OutDetail(Form("       - %-*s : Efficiency",fParameters.MaxNameLength(),"Parameter"));
	for (unsigned i = 0; i < fParameters.Size(); ++i) {
		if (fParameters[i]->Fixed())
			continue;
		BCLog::OutDetail(Form("         %-*s :     %4.1f %%",fParameters.MaxNameLength(),fParameters[i]->GetName().data(), 100.*efficiencies[i]));
	}

	
	// save if improved the log posterior
	if (fMCMCBestFitParameters.empty() || fMCMCprobMax[probmaxindex] > fMCMCLogMaximum) {
		fMCMCLogMaximum = fMCMCprobMax[probmaxindex];
		fMCMCBestFitParameters = fMCMCxMean[probmaxindex];
	}

	BCLog::OutDetail(" --> Global mode from MCMC:");
	BCLog::OutDebug(Form(" --> Posterior value: %g", fMCMCprobMax[probmaxindex]));
	PrintParameters(fMCMCBestFitParameters,BCLog::OutDetail);

	// reset counter
	fMCMCCurrentIteration = -1;

	// reset current chain
	fMCMCCurrentChain = -1;

	// set flags
	fMCMCFlagRun = true;

	return 1;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCResetRunStatistics()
{
	fMCMCNTrials     = 0;

	fMCMCNIterations.assign(fMCMCNChains,0);
	fMCMCprobMean.assign(fMCMCNChains,0);
	fMCMCprobVar.assign(fMCMCNChains,0);
	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(fParameters.Size(),0));

	// reset marginalized distributions
	for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
		if (fH1Marginalized[i])
			fH1Marginalized[i]->Reset();

	for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
		for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
			if (fH2Marginalized[i][j])
				fH2Marginalized[i][j]->Reset();

	fMCMCRValue = 100;
}

// --------------------------------------------------------
int BCEngineMCMC::AddParameter(const char * name, double min, double max, const char * latexname)
{
	// todo memory leak:
	//   ==8243== 64 (44 direct, 20 indirect) bytes in 1 blocks are definitely lost in loss record 20,913 of 29,122
	//   ==8243==    at 0x402C9B4: operator new(unsigned int) (in /usr/lib/valgrind/vgpreload_memcheck-x86-linux.so)
	//   ==8243==    by 0x409DED1: BCEngineMCMC::AddParameter(char const*, double, double) (BCEngineMCMC.cxx:1480)
	return AddParameter(new BCParameter(name, min, max, latexname));
}

// --------------------------------------------------------
int BCEngineMCMC::AddParameter(BCParameter * par)
{
	return fParameters.Add(par);
}

// --------------------------------------------------------
int BCEngineMCMC::AddObservable(const char * name, double min, double max, ObservableFunction fn, const char * latexname)
{
	return AddObservable(new BCObservable(name, min, max, fn, latexname));
}

// --------------------------------------------------------
int BCEngineMCMC::AddObservable(BCObservable * par)
{
	return fObservables.Add(par);
}

// --------------------------------------------------------
void BCEngineMCMC::CalculateObservables() {
	for (unsigned i = 0; i < fMCMCNChains; ++i )
		for (unsigned j = 0; j < fObservables.Size(); ++j)
			fMCMCObservables[i][j] = fObservables[j] -> Evaluate(fMCMCx[i]);
}


// --------------------------------------------------------
void BCEngineMCMC::MCMCInitializeMarkovChains()
{
	// evaluate function at the starting point
	for (unsigned j = 0; j < fMCMCNChains; ++j) {
		double l = LogEval(fMCMCx[j]);
		fMCMCprob[j] = (std::isfinite(l)) ? l : -std::numeric_limits<double>::max();
	}
}

// --------------------------------------------------------
void BCEngineMCMC::ResetResults()
{
	// reset variables
	fMCMCCurrentIteration = -1;
	fMCMCCurrentChain = -1;
	fMCMCNIterations.clear();
	fMCMCNTrialsTrue.clear();
	fMCMCNTrials = 0;
	fMCMCTrialFunctionScaleFactor.clear();
	fMCMCprobMean.clear();
	fMCMCprobVar.clear();
	fMCMCxMean.clear();
	fMCMCxVar.clear();
	fMCMCx.clear();
	fMCMCObservables.clear();
	fMCMCprob.clear();
	fMCMCxMax.clear();
	fMCMCprobMax.clear();
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCRValueParameters.clear();
	fMCMCRValue = std::numeric_limits<double>::max();

	for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
		if (fH1Marginalized[i])
			delete fH1Marginalized[i];

	for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
		for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
			if (fH2Marginalized[i][j])
				delete fH2Marginalized[i][j];

	// clear plots
	fH1Marginalized.clear();
	fH2Marginalized.clear();

	// reset flags
	fMCMCFlagPreRun = true;
	fMCMCFlagRun = false;

	fMCMCBestFitParameters.clear();
	fMCMCLogMaximum = -std::numeric_limits<double>::max();
	fMarginalModes.clear();
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCInitialize()
{
	// resource allocation must be done only by one thread
	// reset values
	ResetResults();

	// free memory for vectors
	fMCMCNIterations.assign(fMCMCNChains, 0);
	fMCMCprobMean.assign(fMCMCNChains, 0);
	fMCMCprobVar.assign(fMCMCNChains, 0);
	fMCMCprob.assign(fMCMCNChains, -1.0);
	fMCMCprobMax.assign(fMCMCNChains, -1.0);

	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(fParameters.Size(), 0));
	fMCMCNTrials = 0;
	fMCMCxMax.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
	fMCMCxMean.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
	fMCMCxVar.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));

	fMCMCRValueParameters.assign(fParameters.Size(), 0);

	SyncThreadStorage();

   if (fMCMCTrialFunctionScaleFactorStart.size() == 0 or fMCMCTrialFunctionScaleFactorStart.size()!=fParameters.Size())
		 fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,std::vector<double>(fParameters.Size(), 1.0));
   else
		 fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,fMCMCTrialFunctionScaleFactorStart);

   // set initial position
   if (fMCMCFlagInitialPosition == 2) { // user defined points
		 // check the length of the array of initial positions
		 if (fMCMCInitialPosition.size() < fMCMCNChains) {
			 BCLog::OutError("BCEngine::MCMCInitialize : Too few initial positions provided.");
			 fMCMCFlagInitialPosition = 1;
		 } else {
			 // check the boundaries
			 for (unsigned j = 0; j < fMCMCNChains; ++j) {
				 if (fMCMCInitialPosition[j].size() == fParameters.Size()) {
					 for (unsigned i = 0; i < fParameters.Size(); ++i)
						 if (!fParameters[i]->IsValid(fMCMCInitialPosition[j][i])) {
							 BCLog::OutError(TString::Format("BCEngine::MCMCInitialize : Initial position of parameter %u (\"%s\") is out of boundaries in chain %u.",i,fParameters[i]->GetName().data(),j).Data());
							 fMCMCFlagInitialPosition = 1;
						 }
				 } else {
					 BCLog::OutError(Form("BCEngine::MCMCInitialize : Too few coordinates provdided in initial position for chain %i",j));
					 fMCMCFlagInitialPosition = 1;
				 }
			 }
		 }

		 if (fMCMCFlagInitialPosition == 1)
			 BCLog::OutError("BCEngine::MCMCInitialize : Using random initial positions instead of user-provided ones.");
		 else {
			 fMCMCx = fMCMCInitialPosition;
			 for (unsigned i = 0; i < fParameters.Size(); ++i)
				 if (fParameters[i]->Fixed())
					 for (unsigned j = 0; j < fMCMCNChains; ++j)
						 fMCMCx[j][i] = fParameters[i]->GetFixedValue();
		 }
	 }
	 
	 if (fMCMCFlagInitialPosition == 0) { // center of the interval
		 fMCMCx.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
		 for (unsigned j = 0; j < fMCMCNChains; ++j)
			 for (unsigned i = 0; i < fParameters.Size(); ++i)
				 fMCMCx[j][i] = (fParameters[i]->Fixed()) ? fParameters[i]->GetFixedValue() : 0.5 * (fParameters[i]->GetLowerLimit() + fParameters[i]->GetUpperLimit());

	 } else { // random number (default)
		 fMCMCx.assign(fMCMCNChains,std::vector<double>(fParameters.Size(),0));
		 for (unsigned j = 0; j < fMCMCNChains; ++j)
			 for (unsigned i = 0; i < fParameters.Size(); ++i)
				 fMCMCx[j][i] = (fParameters[i]->Fixed()) ? fParameters[i]->GetFixedValue() : fParameters[i]->GetLowerLimit() + fMCMCThreadLocalStorage[j].rng->Rndm() * (fParameters[i]->GetUpperLimit()-fParameters[i]->GetLowerLimit());
	 }
	 
	 // initialize user-defined observables
	 fMCMCObservables.assign(fMCMCNChains,std::vector<double>(fObservables.Size(),0));

   // define 1-dimensional histograms for marginalization
   bool fillAny = false;
	 fH1Marginalized.assign(fParameters.Size()+fObservables.Size(),NULL);

   for(unsigned i = 0; i < fParameters.Size(); ++i)
		 if (fParameters[i]->FillHistograms() && !fParameters[i]->Fixed()) {
			 fH1Marginalized[i] = fParameters[i] -> CreateH1(TString::Format("h1_%s_parameter_%i", GetSafeName().data() ,i));
			 fH1Marginalized[i] -> SetStats(kFALSE);
			 fillAny = true;
		 }
   for(unsigned i = 0; i < fObservables.Size(); ++i)
		 if (fObservables[i]->FillHistograms()) {
			 fH1Marginalized[i+fParameters.Size()] = fObservables[i] -> CreateH1(TString::Format("h1_%s_observable_%i", GetSafeName().data() ,i));
			 fH1Marginalized[i+fParameters.Size()] -> SetStats(kFALSE);
			 fillAny = true;
		 }
	 
   // if filling no histograms, set H1 vector to zero size, implies no 2D histograms either
   if (!fillAny)
		 fH1Marginalized.clear();
   else {
		 fH2Marginalized.assign(fParameters.Size()+fObservables.Size(),std::vector<TH2D*>(fParameters.Size()+fObservables.Size(),0));
		 // define 2-dimensional histograms for marginalization
		 for(unsigned i = 0; i < fParameters.Size(); ++i)
			 if (fParameters[i]->FillHistograms() and !fParameters[i]->Fixed()) {
				 // paramater vs parameter
				 for (unsigned j = i + 1; j < fParameters.Size(); ++j)
					 if (fParameters[j]->FillHistograms() and !fParameters[j]->Fixed()) {
						 fH2Marginalized[i][j] = fParameters[i] -> CreateH2(Form("h2_%s_parameters_%i_vs_%i", GetSafeName().data(), i, j), fParameters[j]);
						 fH2Marginalized[i][j] -> SetStats(kFALSE);
					 }
				 // user-defined observable vs parameter
				 for (unsigned j = 0; j < fObservables.Size(); ++j)
					 if (fObservables[j]->FillHistograms()) {
						 fH2Marginalized[i][j+fParameters.Size()] = fParameters[i] -> CreateH2(Form("h2_%s_par_%i_vs_obs_%i", GetSafeName().data(), i, j), fObservables[j]);
						 fH2Marginalized[i][j+fParameters.Size()] -> SetStats(kFALSE);
					 }
			 }
		 // user-defined observable vs user-defined observable
		 for(unsigned i = 0; i < fObservables.Size(); ++i)
			 if (fObservables[i]->FillHistograms())
				 for (unsigned j = i + 1; j < fObservables.Size(); ++j)
					 if (fObservables[j]->FillHistograms()) {
						 fH2Marginalized[i+fParameters.Size()][j+fParameters.Size()] = fObservables[i] -> CreateH2(Form("h2_%s_observables_%i_vs_%i", GetSafeName().data(), i, j), fObservables[j]);
						 fH2Marginalized[i+fParameters.Size()][j+fParameters.Size()] -> SetStats(kFALSE);
					 }

   }
   return 1;
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCIterationInterface()
{
   // do user defined stuff
   MCMCUserIterationInterface();
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintSummary()
{
   // model summary
   BCLog::OutSummary(Form("Model : %s", fName.data()));
   BCLog::OutSummary(Form("Number of parameters : %u", GetNParameters()));
   BCLog::OutSummary("Parameters:");

   // parameter summary
   for (unsigned i = 0; i < GetNParameters(); i++)
      fParameters[i]->PrintSummary();

	 if (GetNObservables()>0) {
		 BCLog::OutSummary(Form("Number of observables : %u", GetNObservables()));
		 BCLog::OutSummary("Observables:");

		 // observable summary
		 for (unsigned i = 0; i < GetNObservables(); i++)
			 fObservables[i]->PrintSummary();
	 }

   // best fit parameters
   if ( !GetBestFitParameters().empty()) {
     BCLog::OutSummary(Form("Log of the maximum posterior: %f", GetLogMaximum()));
		 BCLog::OutSummary("Best fit results:");

		 for (unsigned i = 0; i < GetNVariables(); i++) {
			 if (i < GetNParameters() and GetParameter(i)->Fixed() )
				 BCLog::OutSummary(Form(" %s = %.*f (fixed)",  GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(),GetBestFitParameter(i)));
			 else
				 BCLog::OutSummary(Form(" %s = %.*f (global)", GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(),GetBestFitParameter(i)));
			 
			 if ( fMarginalModes.size() == GetNVariables())
				 BCLog::OutSummary(Form(" %s = %.*f (marginalized)", GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(), GetBestFitParametersMarginalized()[i]));
		 }
   }
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintResults(const char * file) {
	// open file
	std::ofstream ofi(file);

	// check if file is open
	if (!ofi.is_open()) {
		BCLog::OutError(Form("Couldn't open file %s.", file));
		return;
	}

	PrintSummaryToStream(ofi);
	PrintBestFitToStream(ofi);
	PrintMarginalizationToStream(ofi);
	
	// close file
	ofi.close();

}

// ---------------------------------------------------------
void BCEngineMCMC::PrintSummaryToStream(std::ofstream & ofi) {
	ofi << std::endl
			<< " -----------------------------------------------------" << std::endl
			<< " Summary" << std::endl
			<< " -----------------------------------------------------" << std::endl
			<< std::endl;
	
	ofi << " Model summary" << std::endl
			<< " =============" << std::endl
			<< " Model: " << GetName() << std::endl
			<< " Number of parameters: " << GetNParameters() << std::endl
			<< " List of Parameters and ranges:" << std::endl;

	// Parameters & Observables
	for (unsigned i = 0; i < GetNVariables(); ++i) {
		if (i==GetNParameters())
			ofi << " List of Observables and ranges:" << std::endl;
		ofi << " (" << i << ") " << GetVariable(i)->OneLineSummary() << std::endl;
	}
	ofi << std::endl;
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintBestFitToStream(std::ofstream & ofi) {
	if (GetBestFitParameters().empty()) {
		ofi << "No best fit information available." << std::endl << std::endl;
		return;
	}
		
	unsigned nletters = GetMaximumParameterNameLength();

	ofi << " Best Fit Results" << std::endl
			<< " ===========================" << std::endl
			<< " Log of the maximum posterior: " << GetLogMaximum() << std::endl
			<< " List of parameters and global mode:" << std::endl;

	for (unsigned i = 0; i < GetNVariables(); ++i) {
		ofi << TString::Format(" (%d) %10s \"%*s\" : %.*f", i, (i<GetNParameters()) ? "Parameter" : "Observable",
													 nletters, GetVariable(i)->GetName().data(),
													 GetVariable(i)->GetPrecision(),GetBestFitParameter(i));
		if (i<GetNParameters() and GetParameter(i)->Fixed())
			ofi << " (fixed)";
		ofi << std::endl;
	}
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintMarginalizationToStream(std::ofstream & ofi) {
	if (fMCMCFlagRun)
		ofi << " Results of the marginalization" << std::endl
				<< " ==============================" << std::endl;

	// give warning if MCMC did not converge
	if (!MCMCGetNIterationsConvergenceGlobal()<=0 && fMCMCFlagRun)
		ofi << " WARNING: the Markov Chain did not converge!" << std::endl
				<< " Be cautious using the following results!" << std::endl
				<< std::endl;
	
	ofi << " List of parameters and properties of the marginalized" << std::endl
			<< " distributions:" << std::endl;

	for (unsigned i = 0; i < GetNVariables(); ++i) {
		if ( ! GetVariable(i)->FillHistograms())
			continue;
		
		unsigned prec = GetVariable(i) -> GetPrecision();

		ofi << "  (" << i << ") " << ((i<GetNParameters()) ? "Parameter" : "Observable")
				<< " \"" << GetVariable(i)->GetName() << "\":";
		
		if (!MarginalizedHistogramExists(i)) {
			if (i<GetNParameters() and GetParameter(i)->Fixed())
				ofi << TString::Format(" fixed at %.*g.",prec,GetParameter(i)->GetFixedValue()) << std::endl;
			else
				ofi << " histogram does not exist." << std::endl;
			continue;
		}
		ofi << std::endl;
		
		// get marginalized histogram
		BCH1D * bch1d = GetMarginalized(i);
		
		ofi << TString::Format("      Mean +- sqrt(V):                %.*g +- %.*g\n", prec,bch1d->GetMean(), prec,bch1d->GetRMS())
				<< TString::Format("      Median +- central 68%% interval: %.*g +  %.*g - %.*g\n",
													 prec,bch1d->GetMedian(), prec,bch1d->GetQuantile(0.84)-bch1d->GetMedian(), prec,bch1d->GetMedian()-bch1d->GetQuantile(0.16))
				<< TString::Format("      (Marginalized) mode:            %.*g\n",  prec, bch1d->GetMode())
				<< TString::Format("       5%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.05))
				<< TString::Format("      10%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.10))
				<< TString::Format("      16%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.16))
				<< TString::Format("      84%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.85))
				<< TString::Format("      90%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.90))
				<< TString::Format("      95%% quantile:                   %.*g\n", prec, bch1d->GetQuantile(0.95));
		
		std::vector<std::vector<double> > v = bch1d -> GetSmallestIntervals(0.68);
		ofi << TString::Format("      Smallest interval%s containing at least 68%% and local mode%s:",(v.size()>1 ? "s":""),(v.size()>1 ? "s":"")) << std::endl;
		for (unsigned j = 0; j < v.size(); ++j)
			ofi << TString::Format("       (%.*g, %.*g) (local mode at %.*g with rel. height %.*g; rel. area %.*g)\n",
														 prec,v[j][0], prec,v[j][1], prec,v[j][3], prec,v[j][2], prec,v[j][4]);
		ofi << std::endl;
	}
	
	if (fMCMCFlagRun) {

		ofi << " Status of the MCMC" << std::endl
				<< " ==================" << std::endl
				<< " Convergence reached:                    ";
		
		if (MCMCGetNIterationsConvergenceGlobal()>0)
			ofi << "yes" << std::endl
					<< " Number of iterations until convergence: "
					<< MCMCGetNIterationsConvergenceGlobal() << std::endl;
		else
			ofi << "no" << std::endl;
		
		ofi << " Number of chains:                       " << MCMCGetNChains() << std::endl
				<< " Number of iterations per chain:         " << MCMCGetNIterationsRun() << std::endl
				<< " Average run efficiencies:" << std::endl;
		
		unsigned nletters = GetMaximumParameterNameLength(false);

		for (unsigned i = 0; i < GetNParameters(); ++i)
			if (GetParameter(i)->Fixed())
				ofi << TString::Format("  (%d) Parameter \"%*s\" : (fixed)", i, nletters, GetParameter(i)->GetName().data()) << std::endl;
			else {
				double eff = 0;
				for (unsigned j = 0; j < fMCMCEfficiencies.size(); ++j)
					eff += fMCMCEfficiencies[j][i] / fMCMCEfficiencies.size();
				ofi << TString::Format("  (%d) Parameter \"%*s\" : %5.2f %%", i, nletters, GetParameter(i)->GetName().data(),eff*100) << std::endl;
			}
	}

	ofi << " -----------------------------------------------------" << std::endl
			<< " Notation:" << std::endl
			<< " Mean        : mean value of the marg. pdf" << std::endl
			<< " Median      : median of the marg. pdf" << std::endl
			<< " Marg. mode  : most probable value of the marg. pdf" << std::endl
			<< " V           : Variance of the marg. pdf" << std::endl
			<< " Quantiles   : most commonly used quantiles" <<std::endl
			<< " -----------------------------------------------------" << std::endl
			<< std::endl;
}



// ---------------------------------------------------------
void BCEngineMCMC::PrintParameters(std::vector<double> const & P, void (*output)(const char *) ) {
	if ( P.size() != fParameters.Size() )
		return;

	for (unsigned i = 0; i < fParameters.Size(); ++i)
		output(TString::Format("          %-*s :   % .*g", fParameters.MaxNameLength(), fParameters[i]->GetName().data(),fParameters[i]->GetPrecision(),P[i]));
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintAllMarginalized1D(const char * filebase) {
	if (fH1Marginalized.size() == 0) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : No marginalized distributions stored.");
		return 0;
	}

	int nh = 0;
	for (unsigned i = 0; i < fH1Marginalized.size(); ++i) {
		std::string name = (i<fParameters.Size()) ? fParameters[i]->GetName() : fObservables[i-fParameters.Size()]->GetName();
		if (BCH1D * h = GetMarginalized(i)) {
			h->Print(Form("%s_1D_%s.pdf", filebase, name.data()));
			nh++;
		}
	}

	return nh;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintAllMarginalized2D(const char * filebase) {
	if (fH2Marginalized.size() == 0) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : Marginalized distributions not stored.");
		return 0;
	}
	
	int nh = 0;
	for (unsigned i = 0; i<fH2Marginalized.size(); ++i) {
		std::string iname = (i<fParameters.Size()) ? fParameters[i]->GetName() : fObservables[i-fParameters.Size()]->GetName();
		for (unsigned j = 1; j<fH2Marginalized[i].size(); ++i) {
			std::string jname = (j<fParameters.Size()) ? fParameters[j]->GetName() : fObservables[j-fParameters.Size()]->GetName();
			
			if (BCH2D * h = GetMarginalized(i,j)) {
				h -> Print(Form("%s_2D_%s_%s",filebase,iname.data(),jname.data()));
				nh++;
			}
		}
	}

	return nh;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintAllMarginalized(const char * file, std::string options1d, std::string options2d, unsigned int hdiv, unsigned int vdiv) {
	if (fH1Marginalized.empty() and fH2Marginalized.empty()) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : Marginalized distributions not stored.");
		return 0;
	}

	 // Find nonempty H1's
   std::vector<BCH1D *> h1;
   for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
		 if (fH1Marginalized[i]) {
			 h1.push_back(GetMarginalized(i));
      // check if histogram exists
			 if (!h1.back())
				 h1.pop_back();
		 }

	 // Find nonempty H2's
   std::vector<BCH2D *> h2;
	 // fill h2 vector in order (par vs par; obs vs par; obs vs obs)
	 for (unsigned i = 0; i<fParameters.Size() && i<fH2Marginalized.size(); ++i) {
		 // parameter vs parameter
		 for (unsigned j = i+1; j<fParameters.Size() && j<fH2Marginalized[i].size(); ++j)
			 if (fH2Marginalized[i][j]) {
				 h2.push_back(GetMarginalized(i,j));
				 // check if histogram exists
				 if (!h2.back())
					 h2.pop_back();
			 }
		 // user-defined observable vs parameter
		 for (unsigned j = 0; j<fObservables.Size() && j+fParameters.Size()<fH2Marginalized[i].size(); ++j)
			 if (fH2Marginalized[i][j+fParameters.Size()]) {
				 h2.push_back(GetMarginalized(i,j+fParameters.Size()));
				 // check if histogram exists
				 if (!h2.back())
					 h2.pop_back();
			 }
	 }
	 // user-defined observable vs user-define observable
	 for (unsigned i = 0; i<fObservables.Size() && i+fParameters.Size()<fH2Marginalized.size(); ++i)
		 for (unsigned j = i+1; j<fObservables.Size() && j+fParameters.Size()<fH2Marginalized[i].size(); ++j)
			 if (fH2Marginalized[i+fParameters.Size()][j+fParameters.Size()]) {
				 h2.push_back(GetMarginalized(i+fParameters.Size(),j+fParameters.Size()));
				 // check if histogram exists
				 if (!h2.back())
					 h2.pop_back();
			 }

	 if (h1.empty() and h2.empty()) {
		 BCLog::OutWarning("BCEngineMCMC::PrintAllMarginalized : No marginalizations to print");
		 return 0;
	 }

	 if (options1d.empty())
		 options1d = "BTsiB3CS1D0pdf0Lmeanmode";
	 if (options2d.empty())
		 options2d = "BTfB3CS1meangmode";

   std::string filename(file);

   // check if file extension does not exist or is not pdf or ps
   if ( (filename.find_last_of(".") == std::string::npos) or
         ((filename.substr(filename.find_last_of(".")+1) != "pdf") and
               (filename.substr(filename.find_last_of(".")+1) != "ps"))) {
      // make it a PDF file
      filename += ".pdf";
   }

   int c_width  = gStyle->GetCanvasDefW(); // default canvas width
   int c_height = gStyle->GetCanvasDefH(); // default canvas height

   if (hdiv > vdiv) {
      if (hdiv > 3) {
         c_width = 1000;
         c_height = 700;
      }
      else {
         c_width = 800;
         c_height = 600;
      }
   }
   else if (hdiv < vdiv) {
      if (hdiv > 3) {
         c_height = 1000;
         c_width = 700;
      }
      else {
         c_height = 800;
         c_width = 600;
      }
   }

   const unsigned nplots = h1.size() + h2.size();

   // give out warning if too many plots
   BCLog::OutSummary(Form("Printing all marginalized distributions (%d x 1D + %d x 2D = %d) into file %s",
                          (int)h1.size(), (int)h2.size(), nplots, filename.c_str()));
   if (nplots > 100)
      BCLog::OutDetail("This can take a while...");

   // setup the canvas and file
   TCanvas c("c", "canvas", c_width, c_height);
   c.Divide(hdiv, vdiv);
	 c.Print(std::string( filename + "[").c_str());

	 // Draw BCH1D's
   for (unsigned i = 0; i < h1.size(); ++i) {
		 // if current page is full, switch to new page
		 if (i != 0 && i % (hdiv * vdiv) == 0) {
			 c.Print(filename.c_str());
			 c.Clear("D");
		 }

		 // go to next pad
		 c.cd(i % (hdiv * vdiv) + 1);

		 h1[i] -> Draw(options1d);

		 if (i+1 % 100 == 0)
			 BCLog::OutDetail(Form(" --> %d plots done", i+1));
	 }
	 if (h1.size()>0) {
		 c.Print(filename.c_str());
		 c.Clear("D");
	 }

	 // Draw BCH2D's
	 for (unsigned i=0; i<h2.size(); ++i) {
		 // if current page is full, switch to new page
		 if (i != 0 && i % (hdiv * vdiv) == 0) {
			 c.Print(filename.c_str());
			 c.Clear("D");
		 }
				 
		 // go to next pad
		 c.cd(i % (hdiv * vdiv) + 1);
				 
		 h2[i] -> Draw(options2d);
				 
		 if ((i + h1.size() + 1) % 100 == 0)
			 BCLog::OutDetail(Form(" --> %d plots done", i + (int)h1.size() + 1));
	 }
	 if (h2.size()>0) {
		 c.Print(filename.c_str());
		 c.Clear("D");
	 }

   c.Print(std::string( filename + "]").c_str());

   if (nplots > 100 && nplots % 100 != 0)
		 BCLog::OutDetail(Form(" --> %d plots done", nplots));

   // return total number of drawn histograms
   return nplots;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintParameterPlot(const char * filename, int npar, double interval_content, std::vector<double> quantiles) {

	TCanvas * c_par = new TCanvas("c_parplot_init");
	c_par -> Print(Form("%s[",filename));

	int return_val = 1;
	if (npar<=0) { // all parameters on one page, all user-defined observables on the next
		return_val *= PrintParameterPlot(0,GetNParameters(), filename, interval_content, quantiles);
		return_val *= PrintParameterPlot(GetNParameters(),GetNObservables(), filename, interval_content, quantiles);
	}

	else { // npar per page

		// parameters first
		for (unsigned i = 0; i<GetNParameters(); i += npar)
			return_val *= PrintParameterPlot(i,std::min<int>(npar,GetNParameters()-i), filename, interval_content, quantiles);

		// then user-defined observables
		for (unsigned i = GetNParameters(); i<GetNVariables(); i += npar)
			return_val *= PrintParameterPlot(i,std::min<int>(npar,GetNVariables()-i), filename, interval_content, quantiles);
	}

	c_par -> Print(Form("%s]",filename));
	return return_val;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintParameterPlot(unsigned i0, unsigned npar, const char * filename, double interval_content, std::vector<double> quantiles) {

	// if npar==0, print all remaining observables
	unsigned i1 = (npar>0 && npar<=GetNVariables()) ? i0+npar : GetNVariables();

	if (i1 <= i0) {
		BCLog::OutError(Form("BCSummaryTool::PrintParameterPlot : invalid parameter range [%d, %d)",i0,i1));
		return 0;
	}

	// default quantiles
	// TO DO: allow for option to use defualt
	// quantiles, and option for no quantiles
	if (quantiles.empty()) {
		quantiles.push_back(05.00e-2);
		quantiles.push_back(10.00e-2);
		quantiles.push_back(15.87e-2);
		quantiles.push_back(50.00e-2);
		quantiles.push_back(84.13e-2);
		quantiles.push_back(90.00e-2);
		quantiles.push_back(95.00e-2);
	}

	/////////////////////////
	// Gather information

	std::vector<double> x_quantiles;
	std::vector<double> quantile_vals;
	std::vector<double> x_i;
	std::vector<double> x_i_bf;
	std::vector<double> mean;
	std::vector<double> rms;
	std::vector<double> global_mode;
	std::vector<double> local_mode;
	std::vector<double> interval_lo;
	std::vector<double> interval_hi;
	
	for (unsigned i = i0; i < i1; ++i) {
		
		// Global Mode
		x_i_bf.push_back(i);
		global_mode.push_back(GetVariable(i)->PositionInRange(GetBestFitParameter(i)));

		if (!fH1Marginalized[i])
			continue;

		BCH1D * bch1d_temp = GetMarginalized(i);
		if (!bch1d_temp)
			continue;

		x_i.push_back(i);

		// quantiles
		x_quantiles.insert(x_quantiles.end(),quantiles.size(),i);
		for (unsigned j = 0; j < quantiles.size(); ++j)
			quantile_vals.push_back(GetVariable(i)->PositionInRange(bch1d_temp->GetQuantile(quantiles[j])));

		// mean
		mean.push_back(GetVariable(i)->PositionInRange(bch1d_temp->GetMean()));
		rms.push_back(bch1d_temp->GetRMS()/GetVariable(i)->GetRangeWidth());
		 
		// Local Mode
		local_mode.push_back(GetVariable(i)->PositionInRange(bch1d_temp->GetMode()));

		// smallest interval
		std::vector<double> intervals = bch1d_temp->GetSmallestIntervals(interval_content).front();
		if (intervals.size()>3) {
			interval_lo.push_back(fabs(intervals[3]-intervals[0])/GetVariable(i)->GetRangeWidth());
			interval_hi.push_back(fabs(intervals[3]-intervals[1])/GetVariable(i)->GetRangeWidth());
		} else {
			interval_lo.push_back(0);
			interval_hi.push_back(0);
		}
	}
	
	if (x_i.empty() and x_i_bf.empty())
		return 0;

	/////////////////////////
	// Draw it all

	TCanvas * c_par = new TCanvas(TString::Format("c_parplot_%d_%d",i0,i1));
	c_par -> cd();

	// Create, label, and draw axes
	TH2D * hist_axes = new TH2D(TString::Format("h2_axes_parplot_%s_%d_%d",GetSafeName().data(),i0,i1), ";;Scaled range [a.u.]",
															i1-i0, i0-0.5, i1-0.5, 10, -0.1, 1.1);
	hist_axes -> SetStats(kFALSE);
	hist_axes -> GetXaxis() -> SetLabelOffset(0.015);
	hist_axes -> GetXaxis() -> SetLabelSize(std::max<double>(0.01,0.05*std::min<double>(1,4./hist_axes->GetNbinsX())));
	hist_axes -> GetXaxis() -> SetTickLength(0.0);
	// set bin labels
	for (int i=0; i<hist_axes->GetNbinsX(); ++i)
		hist_axes -> GetXaxis() -> SetBinLabel(i+1, GetVariable(i0+i)->GetLatexName().c_str());
	hist_axes -> Draw("");

	// Draw lines
	TLine * line = new TLine();
	line -> SetLineColor(kBlack);
	line -> SetLineStyle(1);
	line -> SetLineWidth(2);
	line -> DrawLine(hist_axes->GetXaxis()->GetXmin(), 0.0, hist_axes->GetXaxis()->GetXmax(), 0.0);
	line -> DrawLine(hist_axes->GetXaxis()->GetXmin(), 1.0, hist_axes->GetXaxis()->GetXmax(), 1.0);

	// Mark parameter ranges
	TLatex * latex = new TLatex();
	latex -> SetTextSize(0.02);
	latex -> SetTextAlign(11);
	latex -> DrawLatex(hist_axes->GetXaxis()->GetXmax(),  1.03, "  Par. max.");
	latex -> SetTextAlign(13);
	latex -> DrawLatex(hist_axes->GetXaxis()->GetXmax(), -0.03, "  Par. min.");
	latex -> SetTextAlign(21);
	for (unsigned i = i0; i < i1; ++i) {
		latex -> SetTextAlign(21);
		latex->DrawLatex((double)i,  1.03, Form("%+.*g", GetVariable(i)->GetPrecision(),GetVariable(i)->GetUpperLimit()));
		latex -> SetTextAlign(23);
		latex->DrawLatex((double)i, -0.03, Form("%+.*g", GetVariable(i)->GetPrecision(),GetVariable(i)->GetLowerLimit()));
	}

	// create legend
	TLegend * legend = new TLegend(0.1, 0.91, 0.9, 0.99);
	legend -> SetBorderSize(0);
	legend -> SetFillColor(0);
	legend -> SetNColumns(2);

	if (!x_i.empty()) {

		// Smallest Interval
		std::vector<double> x_i_err(x_i.size(),0.5);
		TGraphAsymmErrors * graph_intervals = new TGraphAsymmErrors(x_i.size(), x_i.data(), local_mode.data(), x_i_err.data(), x_i_err.data(), interval_lo.data(), interval_hi.data());
		graph_intervals->SetFillColor(kYellow);
		graph_intervals->SetLineStyle(2);
		graph_intervals->SetLineColor(kRed);
		graph_intervals->SetMarkerSize(0);
		graph_intervals->DrawClone("SAME2"); // draw area
		//set y-error zero, to draw line at local mode
		for (int i = 0; i < graph_intervals->GetN(); ++i)
			graph_intervals->SetPointError(i, 0.5, 0.5, 0.0, 0.0);
		graph_intervals->Draw("SAMEZ"); // draw local mode

		// Quantiles graph
		if (!quantiles.empty()) {
			std::vector<double> quantiles_err(x_quantiles.size(),0.5);
			TGraphErrors * graph_quantiles = new TGraphErrors(x_quantiles.size(), x_quantiles.data(), quantile_vals.data(), quantiles_err.data(), 0);
			graph_quantiles->SetMarkerSize(0);
			graph_quantiles->SetLineColor(38);
			graph_quantiles->SetLineStyle(2);
			graph_quantiles->Draw("SAMEZ");
			std::string quantiles_text = "Quantiles (";
			for (unsigned i=0; i<quantiles.size()-1; ++i)
				quantiles_text += Form("%.0f%%, ",quantiles[i]*100);
			quantiles_text += (quantiles.size()>0) ? Form("%.0f%%)",quantiles.back()*100) : "none)";
			legend -> AddEntry(graph_quantiles, quantiles_text.c_str(), "L");
		}

		// Means & RMSs
		TGraphErrors * graph_mean = new TGraphErrors(x_i.size(), x_i.data(), mean.data(), 0, rms.data());
		graph_mean->SetMarkerColor(kBlack);
		graph_mean->SetMarkerStyle(20);
		graph_mean->Draw("SAMEP");

		legend -> AddEntry(graph_mean, "Mean and RMS", "LEP");
		legend -> AddEntry(graph_intervals, TString::Format("Smallest %.0f%% interval and local mode",100.*interval_content), "FL");
	}

	// Global Modes
	if (!x_i_bf.empty()) {
		TGraph * graph_mode = new TGraph(x_i_bf.size(), x_i_bf.data(), global_mode.data());
		graph_mode->SetMarkerColor(kRed);
		graph_mode->SetMarkerStyle(20);
		graph_mode->Draw("SAMEP");
		legend->AddEntry(graph_mode, "Global mode", "P");
	}

	legend->Draw("SAME");
	gPad->RedrawAxis();
	c_par->Print(filename);

	// no error
	return 1;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintCorrelationMatrix(const char * filename) {
   // create histogram
	TH2D * hist_corr = new TH2D(TString::Format("hist_correlation_matrix_%s",GetSafeName().data()),";;",GetNVariables(), -0.5, GetNVariables()-0.5, GetNVariables(), -0.5, GetNVariables()-0.5);
   hist_corr -> SetStats(false);
   hist_corr -> GetXaxis() -> SetTickLength(0.0);
   hist_corr -> GetYaxis() -> SetTickLength(0.0);
	 hist_corr -> GetXaxis() -> SetLabelSize(0);
	 hist_corr -> GetYaxis() -> SetLabelSize(0);

   for (unsigned i = 0; i < GetNVariables(); ++i) {
	 	 hist_corr->SetBinContent(i+1, GetNVariables()-i, 1);
	 	 // Fill Plot
	 	 if (i>=fH2Marginalized.size())
	 		 continue;
	 	 for (unsigned j = i+1; j < GetNVariables(); ++j)
	 		 if (j<fH2Marginalized[i].size() and fH2Marginalized[i][j]) {
	 			 hist_corr -> SetBinContent(i+1, GetNVariables()-j, fH2Marginalized[i][j]->GetCorrelationFactor());
	 			 hist_corr -> SetBinContent(j+1, GetNVariables()-i, fH2Marginalized[i][j]->GetCorrelationFactor());
	 		 }
	 }

   // print to file
   TCanvas * c_corr = new TCanvas("c_corr_matrix");
   c_corr -> cd();

	 double text_size = std::max<double>(0.005,0.02*std::min<double>(1.,5./GetNVariables()));

	 TLatex * xlabel = new TLatex();
	 xlabel -> SetTextFont(62);
	 xlabel -> SetTextSize(text_size);
	 xlabel -> SetTextAlign(22);

	 TLatex * ylabel = new TLatex();
	 ylabel -> SetTextFont(62);
	 ylabel -> SetTextSize(text_size);
	 ylabel -> SetTextAlign(22);
	 ylabel -> SetTextAngle(90);

   hist_corr -> Draw("AXIS");

	 // Draw labels and square colors for correlations
	 TBox * bcorr = new TBox();
	 for (int i = 1; i <= hist_corr->GetNbinsX(); ++i) {
		 // labels
		 xlabel -> DrawLatex(hist_corr->GetXaxis()->GetBinCenter(i),
												 hist_corr->GetYaxis()->GetXmax()+20e-2,
												 GetVariable(i-1)->GetLatexName().c_str());

		 ylabel -> DrawLatex(hist_corr->GetXaxis()->GetXmin()-20e-2,
												 hist_corr->GetYaxis()->GetBinCenter(GetNVariables()-i+1),
												 GetVariable(i-1)->GetLatexName().c_str());
		 
		 //squares
	 	 for (int j = 1; j <= hist_corr->GetNbinsY(); ++j) {
	 		 if (i==hist_corr->GetNbinsY()-j+1)
	 			 bcorr ->SetFillColor(kWhite);
	 		 else if (hist_corr->GetBinContent(i,j) > 0.5)
	 			 bcorr -> SetFillColor(kGreen);
	 		 else if (hist_corr->GetBinContent(i,j) >= -0.5)
	 			 bcorr -> SetFillColor(kYellow);
	 		 else 
	 			 bcorr -> SetFillColor(kRed);
	 		 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(i), hist_corr->GetYaxis()->GetBinLowEdge(j),
	 											hist_corr->GetXaxis()->GetBinUpEdge(i),  hist_corr->GetYaxis()->GetBinUpEdge(j));
	 	 }
	 }

	 // write numbers in
   hist_corr -> Draw("text same");

	 // Blank out empty squares
	 bcorr -> SetFillColor(kWhite);
   for (unsigned i = 0; i < GetNVariables(); ++i)
	 	 for (unsigned j = i+1; j < GetNVariables(); ++j)
	 		 if (i>=fH2Marginalized.size() or j>=fH2Marginalized[i].size() or !fH2Marginalized[i][j]) {
	 			 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(i+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-j),
	 												hist_corr->GetXaxis()->GetBinUpEdge (i+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-j));
	 			 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(j+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-i),
	 												hist_corr->GetXaxis()->GetBinUpEdge (j+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-i));
	 		 }
	 
   // redraw top and right lines
	 TLine * lA = new TLine();
   lA -> DrawLine(hist_corr->GetXaxis()->GetXmin(),hist_corr->GetYaxis()->GetXmax(),hist_corr->GetXaxis()->GetXmax(),hist_corr->GetYaxis()->GetXmax());
   lA -> DrawLine(hist_corr->GetXaxis()->GetXmax(),hist_corr->GetYaxis()->GetXmin(),hist_corr->GetXaxis()->GetXmax(),hist_corr->GetYaxis()->GetXmax());
	 // draw line between parameters and user-defined observables
	 if (GetNObservables()>0) {
	 	 lA -> DrawLine(hist_corr->GetXaxis()->GetXmin()-0.40,hist_corr->GetYaxis()->GetBinLowEdge(hist_corr->GetNbinsY()-GetNParameters()+1),
	 									hist_corr->GetXaxis()->GetBinUpEdge(GetNParameters()),hist_corr->GetYaxis()->GetBinLowEdge(hist_corr->GetNbinsY()-GetNParameters()+1));
	 	 lA -> DrawLine(hist_corr->GetXaxis()->GetBinUpEdge(GetNParameters()),hist_corr->GetYaxis()->GetBinLowEdge(hist_corr->GetNbinsY()-GetNParameters()+1),
	 									hist_corr->GetXaxis()->GetBinUpEdge(GetNParameters()),hist_corr->GetYaxis()->GetXmax()+0.45);
	 }

   gPad->RedrawAxis();
   c_corr->Print(filename);

	 delete lA;
   delete hist_corr;
   delete c_corr;

   // no error
   return 1;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintCorrelationPlot(const char * filename) {

	// Array of indices for which any maginalizations were stored
	std::vector<unsigned> I;
	for (unsigned i = 0; i < fH1Marginalized.size(); ++i) {
		if (fH1Marginalized[i])
			I.push_back(i);
		else 
			for (unsigned j = i+1; j < fH2Marginalized[i].size(); ++j)
				if (fH2Marginalized[i][j]) {
					I.push_back(i);
					break;
				}
	}
	
	if (I.size() == 0)
		return 0;
	
	TCanvas * c = new TCanvas("c_correlation_plot");
	c->cd();
	
	double margin = 0.1;
	double padsize = (1 - 2*margin) / I.size();

	// array with pads holding the histograms
	std::vector<std::vector<TPad*> > pad (I.size(), std::vector<TPad*>(I.size(),0));
	
	// position of pads
	double xlow, xup, ylow, yup;
	double margintop    = 0.01;
	double marginbottom = margintop;
	double marginleft   = 4*margintop;
	double marginright  = marginleft;
	
	TLatex * ylabel = new TLatex();
	ylabel -> SetTextFont(62);
	ylabel -> SetTextSize(1e-1/I.size());
	ylabel -> SetTextAlign(22);			// TODO: set to 32, if latex names too long
	ylabel -> SetNDC();
	ylabel -> SetTextAngle(90);			// TODO: set to 80, if latex names too long
	
	TLatex * xlabel = new TLatex();
	xlabel -> SetTextFont(62);
	xlabel -> SetTextSize(1e-1/I.size());
	xlabel -> SetTextAlign(22);			// TODO: set to 12, if latex names too long
	xlabel -> SetNDC();
	xlabel -> SetTextAngle(0);			// TODO: set to 350, if latex names too long

	// Box + Text for empty squares:
	TBox * box_na = new TBox();
	box_na -> SetLineWidth(1);
	box_na -> SetLineColor(kGray+1);
	box_na -> SetFillColor(kWhite);
	TText * text_na = new TText();
	text_na -> SetTextFont(42);
	text_na -> SetTextAlign(22);
	text_na -> SetTextSize(8e-1/I.size());
	text_na -> SetTextColor(kGray+1);

	// drawing all histograms
	for (unsigned i = 0; i < I.size(); ++i) {
		xlow = i*padsize + margin;
		xup = xlow + padsize;

		for (unsigned j = i; j < I.size(); ++j) {
			yup = 1. - j*padsize - margin;
			ylow = yup - padsize;

			// preparing the pad
			pad[i][j] =  new TPad(TString::Format("pad_correlation_plots_%d_%d",i,j), "", xlow, ylow, xup, yup);
			pad[i][j] -> SetMargin(marginleft,marginright,marginbottom,margintop);
			pad[i][j] -> SetFillColor(kWhite);
			pad[i][j] -> Draw();
			pad[i][j] -> cd();

			// get the histogram
			BCH1D * bh1 = 0;
			BCH2D * bh2 = 0;
			TH1   * hh  = 0;

			if (i==j) {
				bh1 = GetMarginalized(I[i]);
				hh = (bh1) ? bh1->GetHistogram() : 0;
			} else {
				bh2 = GetMarginalized(I[i],I[j]);
				hh = (bh2) ? bh2->GetHistogram() : 0;
			}
			
			if (!bh1 and !bh2) { // if the histogram is not available, draw "N/A"
				pad[i][j] -> SetFillColor(kGray);
				box_na -> DrawBox(marginleft,marginbottom,1.-marginright,1.-margintop);
				text_na -> DrawText(.5,.5,"N/A");

			}	else {									// draw histogram
				if (bh1)
					bh1 -> Draw("BTsiB3CS1D0");
				else
					bh2 -> Draw("BTfB3CS1nL");

				hh -> GetXaxis() -> SetLabelSize(0);
				hh -> GetYaxis() -> SetLabelSize(0);
				hh -> GetXaxis() -> SetTitleSize(0);
				hh -> GetYaxis() -> SetTitleSize(0);
			}

			c->cd();
				
			if(i==0)								// y-axis labels
				ylabel -> DrawLatex(margin*(1-8*ylabel->GetTextSize()), yup-padsize/2., GetVariable(I[j])->GetLatexName().c_str());
			if(j==I.size()-1)				// x-axis labels
				xlabel -> DrawLatex(xlow+padsize/2., margin*(1-8*xlabel->GetTextSize()), GetVariable(I[i])->GetLatexName().c_str());
		}
	}
	
	// gPad->RedrawAxis();
	c->Print(filename);
	
	return 1;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintParameterLatex(const char * filename) {
	// open file
	std::ofstream ofi(filename);
	ofi.precision(3);

	// check if file is open
	if(!ofi.is_open()) {
		std::cerr << "Couldn't open file " << filename <<std::endl;
		return 0;
	}

	const char * blank = "---";
	unsigned texwidth = 15;

	// print table
	ofi	<< "\\documentclass[11pt, a4paper]{article}\n\n"
			<< "\\usepackage[landscape]{geometry}\n\n"
			<< "\\begin{document}\n\n"
			<< "  \\begin{table}[ht!]\n\n"
			<< "    \\begin{center}\n\n"
			<< "      \\begin{tabular}{llllllll}\n\n"
			<< "        Parameter &\n"
			<< TString::Format("        %*s & %*s & %*s & %*s & %*s & %*s & %*s\\\\\n",
												 texwidth, "Mean",
												 texwidth, "RMS",
												 texwidth, "Gl. Mode",
												 texwidth, "Mode",
												 texwidth, "Median",
												 texwidth, "16\\% quant.",
												 texwidth, "84\\% quant.")
			<< "        \\hline\n" << std::endl;
	
	for (unsigned i = 0; i < GetNVariables(); ++i) {

		if (i==GetNParameters())		// first user-defined observable
			ofi << "        &\\\\\n"
					<< "        User-defined Observable &\n"
					<< TString::Format("        %*s & %*s & %*s & %*s & %*s & %*s & %*s\\\\\n",
														 texwidth, "Mean",
														 texwidth, "RMS",
														 texwidth, "Gl. Mode",
														 texwidth, "Mode",
														 texwidth, "Median",
														 texwidth, "16\\% quant.",
														 texwidth, "84\\% quant.")
					<< "        \\hline\n" << std::endl;
		
		// formate LaTeX name for LaTeX ('#' -> '\')
		std::string latexName(GetVariable(i)->GetLatexName());
		std::replace(latexName.begin(),latexName.end(),'#','\\');

		ofi << "        \\ensuremath{" << latexName << "} &" << std::endl;

		unsigned prec = GetVariable(i) -> GetPrecision();

		// fixed parameter
		if (i<GetNParameters() and GetParameter(i)->Fixed())
			ofi << TString::Format("        \\multicolumn{7}{c}{--- fixed to %*.*g ---}\\\\\n",texwidth,prec,GetParameter(i)->GetFixedValue()) << std::endl;
		
		// not fixed
		else {
			BCH1D * bch1d = (i<fH1Marginalized.size() and fH1Marginalized[i]) ? GetMarginalized(i) : 0;

			// marginalization exists
			if (bch1d)
				ofi << TString::Format("        %*.*g & %*.*g & %*.*g & %*.*g & %*.*g & %*.*g & %*.*g\\\\\n",
															 texwidth, prec, bch1d->GetMean(),
															 texwidth, prec, bch1d->GetRMS(),
															 texwidth, prec, GetBestFitParameter(i),
															 texwidth, prec, bch1d -> GetMode(),
															 texwidth, prec, bch1d -> GetMedian(),
															 texwidth, prec, bch1d -> GetQuantile(0.16),
															 texwidth, prec, bch1d -> GetQuantile(0.84))
						<< std::endl;
			
			// marginalization does not exist
			else
				ofi << TString::Format("        %*s & %*s & %*.*g & %*s & %*s & %*s & %*s\\\\\n",
															 texwidth, blank,
															 texwidth, blank,
															 texwidth, prec, GetBestFitParameter(i),
															 texwidth, blank,
															 texwidth, blank,
															 texwidth, blank,
															 texwidth, blank)
						<< std::endl;
		}
	}

	ofi	<< "        \\hline\n" << std::endl
			<< "      \\end{tabular}\n" << std::endl
			<< "      \\caption{Summary of the parameter estimates.}\n" << std::endl
			<< "    \\end{center}\n" << std::endl
			<< "  \\end{table}" << std::endl
			<< std::endl
			<< "\\end{document}" << std::endl;
	
	// close file
	ofi.close();
	
	// no error
	return 1;
}



// ---------------------------------------------------------
unsigned BCEngineMCMC::UpdateFrequency(unsigned N) {
	int n = N/10;
	if (n < 100)
		return 100;
	if (n < 10000)
		return 1000;
	return 10000;
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCThreadLocalStorage::MCMCThreadLocalStorage(const unsigned & dim) :
   xLocal(dim, 0.0),
   rng(new TRandom3(0))
{
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCThreadLocalStorage::MCMCThreadLocalStorage(const MCMCThreadLocalStorage & other)    :
   xLocal(other.xLocal),
   rng(new TRandom3(*other.rng))
{
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCThreadLocalStorage & BCEngineMCMC::MCMCThreadLocalStorage::operator = (const MCMCThreadLocalStorage & other)
{
   xLocal = other.xLocal;
   if (rng)
   {
      // call = operator
      *rng = *other.rng;
   }
   else
      rng = new TRandom3(*other.rng);

   return *this;
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCThreadLocalStorage::~MCMCThreadLocalStorage()
{
   delete rng;
}

// ---------------------------------------------------------
void BCEngineMCMC::SyncThreadStorage() {
	if ( fMCMCNChains > fMCMCThreadLocalStorage.size() )
		fRandom -> Rndm();					// fix return value of GetSeed()

	// add storage until equal to number of chains
	while (fMCMCThreadLocalStorage.size() < fMCMCNChains) {
		fMCMCThreadLocalStorage.push_back(MCMCThreadLocalStorage(fParameters.Size()));
		// each chains gets a different seed. fRandom always returns same seed after the fixing done above
		fMCMCThreadLocalStorage.back().rng->SetSeed(fRandom->GetSeed() + fMCMCThreadLocalStorage.size());
	}

	// remove storage until equal to number of chain
	while (fMCMCThreadLocalStorage.size() > fMCMCNChains)
		fMCMCThreadLocalStorage.pop_back();

	// update parameter size for each chain
	for (unsigned i = 0 ; i < fMCMCThreadLocalStorage.size(); ++i)
		fMCMCThreadLocalStorage[i].xLocal.assign(fParameters.Size(), 0.0);
}
