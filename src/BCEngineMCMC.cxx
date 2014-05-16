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
#include "BCParameter.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TTree.h>

#include <math.h>
#include <limits>
#include <cmath>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC()
{
   // set default parameters for the mcmc
   MCMCSetValuesDefault();

   // initialize random number generator
   fRandom = new TRandom3();
   MCMCSetRandomSeed(0);
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesDefault()
{
   fMCMCFlagWriteChainToFile  = false;
   fMCMCFlagWritePreRunToFile = false;
   fMCMCFlagPreRun            = true;
   fMCMCFlagRun               = false;
   fMCMCEfficiencyMin         = 0.15;
   fMCMCEfficiencyMax         = 0.50;
	 fMCMCScaleFactorLowerLimit = 0;
	 fMCMCScaleFactorUpperLimit = std::numeric_limits<double>::max();
   fMCMCFlagInitialPosition   = 1;
   fMCMCNLag                  = 1;
   fMCMCCurrentIteration      = -1;
   fMCMCCurrentChain          = -1;
   fMCMCLogMaximum            = -std::numeric_limits<double>::max();

   MCMCSetValuesDetail();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesQuick()
{
   fMCMCNChains              = 2;
   fMCMCNIterationsMax       = 10000;
   fMCMCNIterationsRun       = 10000;
   fMCMCNIterationsPreRunMin = 500;
   fMCMCFlagInitialPosition  = 1;
   fMCMCRValueUseStrict      = false;
   fMCMCRValueCriterion      = 0.1;
   fMCMCRValueParametersCriterion = 0.1;
   fMCMCNIterationsConvergenceGlobal = -1;
   fMCMCFlagConvergenceGlobal = false;
   fMCMCRValue               = 100;
	 fMCMCNIterationsEfficiencyCheck = 500;
   fMCMCNIterationsUpdate    = 1000;
   fMCMCNIterationsUpdateClear = 5000;
   fMCMCNIterationsUpdateMax = 10000;
   fMCMCFlagOrderParameters  = true;
   fMCMCCurrentIteration     = -1;
   fMCMCCurrentChain         = -1;
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesDetail()
{
   fMCMCNChains              = 5;
   fMCMCNIterationsMax       = 1000000;
   fMCMCNIterationsRun       = 100000;
   fMCMCNIterationsPreRunMin = 500;
   fMCMCFlagInitialPosition  = 1;
   fMCMCRValueUseStrict = false;
   fMCMCRValueCriterion      = 0.1;
   fMCMCRValueParametersCriterion = 0.1;
   fMCMCNIterationsConvergenceGlobal = -1;
   fMCMCFlagConvergenceGlobal = false;
   fMCMCRValue               = 100;
	 fMCMCNIterationsEfficiencyCheck = 500;
   fMCMCNIterationsUpdate    = 1000;
   fMCMCNIterationsUpdateClear = 5000;
   fMCMCNIterationsUpdateMax = 10000;
   fMCMCFlagOrderParameters  = true;
   fMCMCCurrentIteration     = -1;
   fMCMCCurrentChain         = -1;
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetPrecision(BCEngineMCMC::Precision precision)
{
   switch(precision) {
   case BCEngineMCMC::kLow:
      fMCMCNChains              = 1;
      fMCMCNLag                 = 1;
      fMCMCNIterationsMax       = 10000;
      fMCMCNIterationsRun       = 10000;
      fMCMCNIterationsPreRunMin = 100;
			fMCMCNIterationsEfficiencyCheck = 500;
      fMCMCNIterationsUpdate    = 1000;
      fMCMCNIterationsUpdateClear= 5000;
      fMCMCNIterationsUpdateMax = 10000;
      fMCMCRValueCriterion      = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue               = 100;
      break;
   case  BCEngineMCMC::kMedium:
      fMCMCNChains              = 5;
      fMCMCNLag                 = 1;
      fMCMCNIterationsMax       = 100000;
      fMCMCNIterationsRun       = 100000;
      fMCMCNIterationsPreRunMin = 100;
			fMCMCNIterationsEfficiencyCheck = 1000;
      fMCMCNIterationsUpdate    = 1000;
      fMCMCNIterationsUpdateClear = 5000;
      fMCMCNIterationsUpdateMax = 10000;
      fMCMCRValueCriterion      = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue               = 100;
      break;
   case  BCEngineMCMC::kHigh:
      fMCMCNChains                = 10;
      fMCMCNLag                   = 10;
      fMCMCNIterationsMax         = 1000000;
      fMCMCNIterationsRun         = 1000000;
      fMCMCNIterationsPreRunMin   = 100;
      fMCMCNIterationsUpdate      = 1000;
      fMCMCNIterationsUpdateClear = 5000;
      fMCMCNIterationsUpdateMax   = 10000;
      fMCMCRValueCriterion        = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue                 = 100;
      break;
   case  BCEngineMCMC::kVeryHigh:
      fMCMCNChains              = 10;
      fMCMCNLag                 = 10;
      fMCMCNIterationsMax       = 10000000;
      fMCMCNIterationsRun       = 10000000;
      fMCMCNIterationsPreRunMin = 100;
			fMCMCNIterationsEfficiencyCheck = 1000;
      fMCMCNIterationsUpdate    = 1000;
      fMCMCNIterationsUpdateClear = 5000;
      fMCMCNIterationsUpdateMax = 10000;
      fMCMCRValueCriterion      = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue               = 100;
      break;
   }

   // re-initialize
   MCMCInitialize();
}

// ---------------------------------------------------------
BCEngineMCMC::~BCEngineMCMC()
{
   // delete random number generator
   delete fRandom;

   // delete 1-d marginalized distributions
   for (unsigned i = 0; i < fMCMCH1Marginalized.size(); ++i)
      delete fMCMCH1Marginalized[i];
   fMCMCH1Marginalized.clear();

   // delete 2-d marginalized distributions
   for (unsigned i = 0; i < fMCMCH2Marginalized.size(); ++i)
      delete fMCMCH2Marginalized[i];
   fMCMCH2Marginalized.clear();
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC & other)
{
	Copy(other);
}

// ---------------------------------------------------------
void BCEngineMCMC::Copy(const BCEngineMCMC & other)
{
   fMCMCPointerToGetProposalPoint = other.fMCMCPointerToGetProposalPoint;
   fMCMCNChains               = other.fMCMCNChains;
   fMCMCNLag                  = other.fMCMCNLag;
   fMCMCNIterations           = other.fMCMCNIterations;
   fMCMCCurrentIteration      = other.fMCMCCurrentIteration;
   fMCMCCurrentChain          = other.fMCMCCurrentChain;
	 fMCMCNIterationsEfficiencyCheck = other.fMCMCNIterationsEfficiencyCheck;
   fMCMCNIterationsUpdate     = other.fMCMCNIterationsUpdate;
   fMCMCNIterationsUpdateClear= other.fMCMCNIterationsUpdateClear;
   fMCMCNIterationsUpdateMax  = other.fMCMCNIterationsUpdateMax;
   fMCMCNIterationsConvergenceGlobal = other.fMCMCNIterationsConvergenceGlobal;
   fMCMCFlagConvergenceGlobal = other.fMCMCFlagConvergenceGlobal;
   fMCMCNIterationsMax        = other.fMCMCNIterationsMax;
   fMCMCNIterationsRun        = other.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin  = other.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue           = other.fMCMCNTrialsTrue;
   fMCMCNTrials               = other.fMCMCNTrials;
   fMCMCFlagWriteChainToFile  = other.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile = other.fMCMCFlagWritePreRunToFile;
   fMCMCTrialFunctionScaleFactor = other.fMCMCTrialFunctionScaleFactor;
   fMCMCTrialFunctionScaleFactorStart = other.fMCMCTrialFunctionScaleFactorStart;
   fMCMCFlagPreRun            = other.fMCMCFlagPreRun;
   fMCMCFlagRun               = other.fMCMCFlagRun;
   fMCMCInitialPosition       = other.fMCMCInitialPosition;
   fMCMCEfficiencies          = other.fMCMCEfficiencies;
   fMCMCEfficiencyMin         = other.fMCMCEfficiencyMin;
   fMCMCEfficiencyMax         = other.fMCMCEfficiencyMax;
	 fMCMCScaleFactorLowerLimit = other.fMCMCScaleFactorLowerLimit;
	 fMCMCScaleFactorUpperLimit = other.fMCMCScaleFactorUpperLimit;
   fMCMCFlagInitialPosition   = other.fMCMCFlagInitialPosition;
   fMCMCFlagOrderParameters   = other.fMCMCFlagOrderParameters;
   fMCMCPhase                 = other.fMCMCPhase;
   fMCMCx                     = other.fMCMCx;
   fMCMCxMax                  = other.fMCMCxMax;
   fMCMCxMean                 = other.fMCMCxMean;
   fMCMCxVar                  = other.fMCMCxVar;
   fMCMCprob                  = other.fMCMCprob;
   fMCMCprobMax               = other.fMCMCprobMax;
   fMCMCprobMean              = other.fMCMCprobMean;
   fMCMCprobVar               = other.fMCMCprobVar;
   fMCMCRValueUseStrict       = other.fMCMCRValueUseStrict;
   fMCMCRValueCriterion       = other.fMCMCRValueCriterion ;
   fMCMCRValueParametersCriterion = other.fMCMCRValueParametersCriterion;
   fMCMCRValue                = other.fMCMCRValue;
   fMCMCRValueParameters      = other.fMCMCRValueParameters;
   if (other.fRandom)
   {
      fRandom = new TRandom3(*other.fRandom);
   }
   else
   {
      fRandom = NULL;
   }

   fMCMCThreadLocalStorage    = other.fMCMCThreadLocalStorage;

   for (unsigned i = 0; i < other.fMCMCH1Marginalized.size(); ++i) {
      if (other.fMCMCH1Marginalized.at(i))
         fMCMCH1Marginalized.push_back(new TH1D(*(other.fMCMCH1Marginalized.at(i))));
      else
         fMCMCH1Marginalized.push_back(0);
   }

   for (unsigned i = 0; i < other.fMCMCH2Marginalized.size(); ++i) {
      if (other.fMCMCH2Marginalized.at(i))
         fMCMCH2Marginalized.push_back(new TH2D(*(other.fMCMCH2Marginalized.at(i))));
      else
         fMCMCH2Marginalized.push_back(0);
   }

   for (unsigned i = 0; i < other.fMCMCTrees.size(); ++i) {
      fMCMCTrees.push_back(0);
   }

   fMarginalModes = other.fMarginalModes;
   fMCMCBestFitParameters = other.fMCMCBestFitParameters;
   fMCMCLogMaximum = other.fMCMCLogMaximum;
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
   for (unsigned i = 0 ; i < fParameters.Size() ; ++i)
      fParameters[i]->SetNbins(nbins);
}

// --------------------------------------------------------
BCH1D * BCEngineMCMC::MCMCGetH1Marginalized(unsigned index)
{
   if ( !fParameters.ValidIndex(index)) {
      BCLog::OutError(Form("BCEngineMCMC::MCMCGetH1Marginalized. Index %u out of range.", index));
      return 0;
   }

   // use when BCIntegrate acts MCMC marginalization and only some marginals have been computed
   if (!fMCMCH1Marginalized[index]) {
      BCLog::OutWarning(Form("BCEngineMCMC::MCMCGetH1Marginalized: marginal distribution not computed/stored for par. %d", index));
      return 0;
   }

   // set histogram
   BCH1D * hprob = new BCH1D();
   hprob->SetHistogram(fMCMCH1Marginalized[index]);

   if (fMarginalModes.empty())
      fMarginalModes.assign(fParameters.Size(), 0.0);
   fMarginalModes[index] = hprob->GetMode();

   return hprob;
}

// --------------------------------------------------------
BCH2D * BCEngineMCMC::MCMCGetH2Marginalized(unsigned i, unsigned j)
{
   if ( !fParameters.ValidIndex(i)) {
      BCLog::OutError(Form("BCEngineMCMC::MCMCGetH2Marginalized. Index %u out of range.", i));
      return 0;
   }
   if ( !fParameters.ValidIndex(j)) {
      BCLog::OutError(Form("BCEngineMCMC::MCMCGetH2Marginalized. Index %u out of range.", j));
      return 0;
   }
   if (i == j) {
      BCLog::OutError(Form("BCEngineMCMC::MCMCGetH2Marginalized. Called with identical indices %u.", i));
      return 0;
   }

   // swap indices
   if (i > j) {
      unsigned indexTemp = i;
      i = j;
      j = indexTemp;
   }

   // memory layout for n parameters and indices i, j:
   // first (n-1) elements for first parameter vs. all other parameters
   // then (n-2) elements for second parameter and all others etc
   // so the first combination for which i is the lower index is at n*i - i(i+1)/2
   // and the offset is given  by (j-i-1)
   TH2D * h =  fMCMCH2Marginalized.at(GetNParameters() * i - (i * i + 3 * i) / 2 + j - 1);
   if ( !h)
      return 0;

   BCH2D * hprob = new BCH2D();
   hprob->SetHistogram(h);

   return hprob;
}

// ---------------------------------------------------------
const std::vector<double> & BCEngineMCMC::GetBestFitParametersMarginalized() const
{
   if(fMarginalModes.empty())
      BCLog::OutError("BCIntegrate::GetBestFitParameterMarginalized : MCMC not yet run, returning center of the range.");

   return fMarginalModes;
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

   // re-initialize
   MCMCInitialize();
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
	for (unsigned i = 0; i < fParameters.Size(); ++i)
		fParameters[i]->FillHistograms(flag);
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
	if (!fRandom)
		fRandom = new TRandom3();

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

		for (unsigned j = 0; j < fParameters.Size(); ++j)
			fMCMCTrees[i]->Branch(TString::Format("Parameter%i", j), &fMCMCx[i][j], TString::Format("parameter %i/D", j));
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

		// fill each 1-dimensional histogram (if supposed to be filled)
		for (unsigned j = 0; j < fParameters.Size(); ++j)
			if (TH1 * h = fMCMCH1Marginalized[j])
				h -> Fill(fMCMCx[i][j]);
		
		// fill each 2-dimensional histogram (if supposed to be filled)
		unsigned counter = 0;

		for (unsigned j = 0; j < fParameters.Size(); ++j)
			for (unsigned k = j+1; k < fParameters.Size(); ++k) {
				if (TH2D * h = fMCMCH2Marginalized[counter])
					h -> Fill(fMCMCx[i][j],fMCMCx[i][k]);
				counter ++;
			}

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
	
	int nletters = 20;
	
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

		if (fMCMCFlagWritePreRunToFile)
			MCMCInChainWriteChains();

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
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale decreased to %6.2f %%", nletters, fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					}	else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorLowerLimit;
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale now at lower limit (%6.2f %%)",	nletters, fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCScaleFactorLowerLimit));
					}
					
				} else if (fMCMCEfficiencies[ichain][ipar] > fMCMCEfficiencyMax ) {
					// if efficiency too high ...

					if (allEfficient)		 // print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					
					allEfficient = false;

					fMCMCTrialFunctionScaleFactor[ichain][ipar] *= 2;

					if ( fMCMCTrialFunctionScaleFactor[ichain][ipar] < fMCMCScaleFactorUpperLimit ) {
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale increased to %6.2f %%", nletters, fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					} else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorUpperLimit;
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale now at upper limit (%6.2f %%)", nletters, fParameters[ipar]->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, fMCMCScaleFactorUpperLimit));														 
					}
				}
			}
		}
		fMCMCNTrials = 0;
		fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(fParameters.Size(),0));
	}
	if (allEfficient)
		BCLog::OutDetail("     * Efficiency status: Efficiencies within predefined range.");
	else if (!inefficientScalesAdjustable)
		BCLog::OutWarning("     * Efficiency status: Some efficiencies outside predefined range, but scales are at limits.");
	else
		BCLog::OutDetail("     * Efficiency status: Some efficiencies outside predefined range, but maximum number of iterations reached.");


	// continue measuring efficiency
	unsigned NTrialsForEff = fMCMCNTrials;

	if (fMCMCNChains > 1) {
		//////////////////////////////////////////////////
		// Run until all chains have converged

		unsigned nIterationsCheckConvergence = fMCMCNIterationsUpdate;
		if ( nIterationsCheckConvergence > fMCMCNIterationsUpdateMax )
			nIterationsCheckConvergence = fMCMCNIterationsUpdateMax;
		if ( nIterationsCheckConvergence > fMCMCNIterationsUpdateClear )
			nIterationsCheckConvergence = fMCMCNIterationsUpdateClear;

		fMCMCNTrials = fMCMCNIterationsUpdateClear;
	
		fMCMCNIterationsConvergenceGlobal = -1;
		
		if (fMCMCCurrentIteration >= (int)fMCMCNIterationsMax) {
      BCLog::OutWarning(" Convergence never checked !");
      BCLog::OutWarning("   Increase maximum number of iterations in the pre-run /MCMCSetNIterationsMax()/");
      BCLog::OutWarning("   or decrease maximum number of iterations for update  /MCMCSetNIterationsUpdateMax()/");
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

			if (fMCMCFlagWritePreRunToFile)
				MCMCInChainWriteChains();

			MCMCInChainUpdateStatistics();

			if ( fMCMCNTrials % nIterationsCheckConvergence != 0 and fMCMCCurrentIteration < (int)fMCMCNIterationsMax )
				continue;

			MCMCInChainTestConvergenceAllChains();

			if ( fMCMCNIterationsConvergenceGlobal > 0 )
				continue;

			BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));
			BCLog::OutDetail(Form("       - %-*s : R-Value",nletters,"Parameter"));
				
			for (unsigned ipar = 0; ipar < fParameters.Size(); ++ipar) {

				if ( fParameters[ipar]->Fixed() )
					continue;

				if( fMCMCRValueParameters[ipar]-1 < fMCMCRValueParametersCriterion )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f",nletters,fParameters[ipar]->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else if ( fMCMCRValueParameters.at(ipar) != std::numeric_limits<double>::max() )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f <--",nletters,fParameters[ipar]->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else
					BCLog::OutDetail(TString::Format("         %-*s :  MAX_DOUBLE <--",nletters,fParameters[ipar]->GetName().data()));
			}

			if( fMCMCRValue-1 < fMCMCRValueCriterion )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f", fMCMCRValue));
			else if ( fMCMCRValue != std::numeric_limits<double>::max() )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f <--", fMCMCRValue));
			else
				BCLog::OutDetail("       - Log-Likelihood :  MAX_DOUBLE <--");
		}

		fMCMCFlagConvergenceGlobal = (fMCMCNIterationsConvergenceGlobal > 0);

		if ( fMCMCFlagConvergenceGlobal )
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
	else													// only one chain
		if (allEfficient)
			BCLog::OutSummary(" --> All scales adjusted.");
		else
			BCLog::OutSummary(" --> Could not adjust all scales.");

		
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

			if (fMCMCFlagWritePreRunToFile)
				MCMCInChainWriteChains();

			if ( n % nwrite == 0)
				BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, 100.*n/N));
		}
		BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations in pre-run.", fMCMCCurrentIteration));
	}

	// print scale factors and efficiencies
	std::vector<double> scalefactors (fParameters.Size(),0);
	std::vector<double> efficiencies (fParameters.Size(),0);

	BCLog::OutDetail(Form(" --> Average scale factors and efficiencies (measured in %d iterations):",(NTrialsForEff==0) ? fMCMCNIterationsEfficiencyCheck : NTrialsForEff));
	BCLog::OutDetail(Form("       - %-*s : Scale factor    Efficiency",nletters,"Parameter"));
	for (unsigned i = 0; i < fParameters.Size(); ++i) {
		if (fParameters[i]->Fixed())
			continue;
		for (unsigned j = 0; j < fMCMCNChains; ++j) {
			efficiencies[i] += ( NTrialsForEff==0 ) ? fMCMCEfficiencies[j][i] / fMCMCNChains : 1.*fMCMCNTrialsTrue[j][i]/NTrialsForEff/fMCMCNChains;
			scalefactors[i] += fMCMCTrialFunctionScaleFactor[j][i] / fMCMCNChains;
		}
		BCLog::OutDetail(Form("         %-*s :     %6.02f %%        %4.1f %%",nletters,fParameters[i]->GetName().data(), 100.*scalefactors[i], 100.*efficiencies[i]));
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
		
		// fill histograms
		if ( !fMCMCH1Marginalized.empty() or !fMCMCH2Marginalized.empty() )
			MCMCInChainFillHistograms();
		
		// write chain to file
		if ( fMCMCFlagWriteChainToFile )
			MCMCInChainWriteChains();
		
	} // end run

	BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

	// find global maximum
	unsigned probmaxindex = 0;

	// loop over all chains and find the maximum point
	for (unsigned i = 1; i < fMCMCNChains; ++i)
		if (fMCMCprobMax[i] > fMCMCprobMax[probmaxindex])
			probmaxindex = i;
	
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
	for (unsigned i = 0; i < fMCMCH1Marginalized.size(); ++i)
		if (fMCMCH1Marginalized[i])
			fMCMCH1Marginalized[i]->Reset();

	for (unsigned i = 0; i < fMCMCH2Marginalized.size(); ++i)
		if (fMCMCH2Marginalized[i])
			fMCMCH2Marginalized[i]->Reset();

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
	fMCMCNIterations.clear();
	fMCMCNTrialsTrue.clear();
	fMCMCNTrials = 0;
	fMCMCTrialFunctionScaleFactor.clear();
	fMCMCprobMean.clear();
	fMCMCprobVar.clear();
	fMCMCxMean.clear();
	fMCMCxVar.clear();
	fMCMCx.clear();
	fMCMCprob.clear();
	fMCMCxMax.clear();
	fMCMCprobMax.clear();
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCRValueParameters.clear();

	for (unsigned i = 0; i < fMCMCH1Marginalized.size(); ++i)
		if (fMCMCH1Marginalized[i])
			delete fMCMCH1Marginalized[i];

	for (unsigned i = 0; i < fMCMCH2Marginalized.size(); ++i)
		if (fMCMCH2Marginalized[i])
			delete fMCMCH2Marginalized[i];

	// clear plots
	fMCMCH1Marginalized.clear();
	fMCMCH2Marginalized.clear();

	// reset flags
	fMCMCFlagPreRun = true;
	fMCMCFlagRun = false;
	fMCMCFlagConvergenceGlobal = false;

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
	 
   // define 1-dimensional histograms for marginalization
   bool fillAny = false;
   for(unsigned i = 0; i < fParameters.Size(); ++i) {
		 const BCParameter * p = fParameters[i];
		 TH1D * h1 = NULL;
		 if (p->FillHistograms() && ! p->Fixed()) {
			 h1 = new TH1D(TString::Format("h1_%d_parameter_%i", BCLog::GetHIndex() ,i),
										 TString::Format(";%s;", p->GetLatexName().c_str()),
										 p->GetNbins(), p->GetLowerLimit(), p->GetUpperLimit());
			 h1->SetStats(kFALSE);
			 
			 fillAny = true;
		 }
		 fMCMCH1Marginalized.push_back(h1);
   }

   // if filling no histograms, set H1 vector to zero size, implies no 2D histograms either
   if (!fillAny)
		 fMCMCH1Marginalized.clear();
   else {
		 // define 2-dimensional histograms for marginalization
		 for(unsigned i = 0; i < fParameters.Size(); ++i) {
			 BCParameter * p1 =  fParameters[i];
			 for (unsigned j = i + 1; j < fParameters.Size(); ++j) {
				 TH2D * h2 = 0;
				 BCParameter * p2 =  fParameters[j];
				 if (p2->FillHistograms() && p1->FillHistograms() && ! p2->Fixed() && ! p1->Fixed()) {
					 h2 = new TH2D(Form("h2_%d_parameters_%i_vs_%i", BCLog::GetHIndex(), i, j), "",
												 p1->GetNbins(), p1->GetLowerLimit(), p1->GetUpperLimit(),
												 p2->GetNbins(), p2->GetLowerLimit(), p2->GetUpperLimit());
					 
					 // decorate histogram
					 h2->SetXTitle(Form("%s", p1->GetLatexName().data()));
					 h2->SetYTitle(Form("%s", p2->GetLatexName().data()));
					 h2->SetStats(kFALSE);
				 }
				 fMCMCH2Marginalized.push_back(h2);
			 }
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
int BCEngineMCMC::SetMarginalized(unsigned index, TH1D * h)
{
   if(fMCMCH1Marginalized.size() <= index)
      return 0;

   if(h==0)
      return 0;

   if(fMCMCH1Marginalized.size() == index)
      fMCMCH1Marginalized.push_back(h);
   else
      fMCMCH1Marginalized[index]=h;

   return index;
}

// ---------------------------------------------------------
int BCEngineMCMC::SetMarginalized(unsigned index1, unsigned index2, TH2D * h)
{
   unsigned counter = 0;
   unsigned index = 0;

   // search for correct combination
   for(unsigned i = 0; i < fParameters.Size(); i++)
      for (unsigned j = 0; j < i; ++j)
      {
         if(j == index1 && i == index2)
            index = counter;
         counter++;
      }

   if(fMCMCH2Marginalized.size()<=index)
      return 0;

   if(h==0)
      return 0;

   if(fMCMCH2Marginalized.size()==index)
      fMCMCH2Marginalized.push_back(h);
   else
      fMCMCH2Marginalized[index]=h;

   return index;
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
void BCEngineMCMC::SyncThreadStorage()
{
   // always need as many local storage as #chains
   const int n = fMCMCThreadLocalStorage.size() - fMCMCNChains;
   if (n < 0)
   {
      // fix return value of GetSeed()
      fRandom->Rndm();

      for (int i = 0; i < -n; ++i){
         // append one new storage
         fMCMCThreadLocalStorage.push_back(MCMCThreadLocalStorage(fParameters.Size()));
         // each chains gets a different seed
         // We assume that fRandom always returns same seed, presumably as it has generated at least one random number
         fMCMCThreadLocalStorage.back().rng->SetSeed(fRandom->GetSeed() + fMCMCThreadLocalStorage.size());
      }
   }
   else if (n > 0)
   {
      for (int i = 0; i < n; ++i){
         fMCMCThreadLocalStorage.pop_back();
      }
   }

   // update parameter size for each chain
   for (unsigned i = 0 ; i < fMCMCThreadLocalStorage.size(); ++i)
      fMCMCThreadLocalStorage[i].xLocal.assign(fParameters.Size(), 0.0);
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintParameters(std::vector<double> const & P, void (*output)(const char *) ) {
	if ( P.size() != fParameters.Size() )
		return;

	unsigned int nletters = 0;
	for (unsigned i = 0; i < fParameters.Size(); ++i)
		if (fParameters[i]->GetName().size() > nletters)
			nletters = fParameters[i]->GetName().size();

	for (unsigned i = 0; i < fParameters.Size(); ++i)
		output(TString::Format("          %-*s :   %.6g", nletters, fParameters[i]->GetName().data(),P[i]));
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
