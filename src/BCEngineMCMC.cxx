/*
 * Copyright (C) 2008-2012, Daniel Kollar, Kevin Kroeninger, and Daniel Greenwald.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
   fMCMCFlagWriteChainToFile = false;
   fMCMCFlagWritePreRunToFile = false;
   fMCMCFlagPreRun           = true;
   fMCMCFlagRun              = false;
   fMCMCEfficiencyMin        = 0.15;
   fMCMCEfficiencyMax        = 0.50;
   fMCMCFlagInitialPosition  = 1;
   fMCMCNLag                 = 1;
   fMCMCCurrentIteration     = -1;
   fMCMCCurrentChain         = -1;

   fFillErrorBand = false;
   fFitFunctionIndexX = -1;
   fFitFunctionIndexY = -1;
   fErrorBandContinuous = true;
   fErrorBandXY = NULL;
   fErrorBandNbinsX = 100;
   fErrorBandNbinsY = 500;

   fOptimizationMethod = BCEngineMCMC::kOptMinuit;
   fOptimizationMethodMode = BCEngineMCMC::NOptMethods;
   fLogMaximum = -std::numeric_limits<double>::max();

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
   fMCMCRValueUseStrict = false;
   fMCMCRValueCriterion      = 0.1;
   fMCMCRValueParametersCriterion = 0.1;
   fMCMCNIterationsConvergenceGlobal = -1;
   fMCMCFlagConvergenceGlobal = false;
   fMCMCRValue               = 100;
   fMCMCNIterationsUpdate    = 1000;
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
   fMCMCNIterationsUpdate    = 1000;
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
      fMCMCNIterationsUpdate    = 1000;
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
      fMCMCNIterationsUpdate    = 1000;
      fMCMCNIterationsUpdateMax = 10000;
      fMCMCRValueCriterion      = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue               = 100;
      break;
   case  BCEngineMCMC::kHigh:
      fMCMCNChains              = 10;
      fMCMCNLag                 = 10;
      fMCMCNIterationsMax       = 1000000;
      fMCMCNIterationsRun       = 1000000;
      fMCMCNIterationsPreRunMin = 100;
      fMCMCNIterationsUpdate    = 1000;
      fMCMCNIterationsUpdateMax = 10000;
      fMCMCRValueCriterion      = 0.1;
      fMCMCRValueParametersCriterion = 0.1;
      fMCMCRValue               = 100;
      break;
   case  BCEngineMCMC::kVeryHigh:
      fMCMCNChains              = 10;
      fMCMCNLag                 = 10;
      fMCMCNIterationsMax       = 10000000;
      fMCMCNIterationsRun       = 10000000;
      fMCMCNIterationsPreRunMin = 100;
      fMCMCNIterationsUpdate    = 1000;
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
   fMCMCNIterationsUpdate     = other.fMCMCNIterationsUpdate;
   fMCMCNIterationsUpdateMax  = other.fMCMCNIterationsUpdateMax;
   fMCMCNIterationsConvergenceGlobal = other.fMCMCNIterationsConvergenceGlobal;
   fMCMCFlagConvergenceGlobal = other.fMCMCFlagConvergenceGlobal;
   fMCMCNIterationsMax        = other.fMCMCNIterationsMax;
   fMCMCNIterationsRun        = other.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin  = other.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue           = other.fMCMCNTrialsTrue;
   fMCMCNTrialsFalse          = other.fMCMCNTrialsFalse;
   fMCMCFlagWriteChainToFile  = other.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile = other.fMCMCFlagWritePreRunToFile;
   fMCMCTrialFunctionScaleFactor = other.fMCMCTrialFunctionScaleFactor;
   fMCMCTrialFunctionScaleFactorStart = other.fMCMCTrialFunctionScaleFactorStart;
   fMCMCFlagPreRun            = other.fMCMCFlagPreRun;
   fMCMCFlagRun               = other.fMCMCFlagRun;
   fMCMCInitialPosition       = other.fMCMCInitialPosition;
   fMCMCEfficiencyMin         = other.fMCMCEfficiencyMin;
   fMCMCEfficiencyMax         = other.fMCMCEfficiencyMax;
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
   fOptimizationMethod = other.fOptimizationMethod;
   fOptimizationMethodMode = other.fOptimizationMethodMode;
   fBestFitParameters = other.fBestFitParameters;
   fLogMaximum = other.fLogMaximum;
}

// ---------------------------------------------------------
BCEngineMCMC & BCEngineMCMC::operator = (const BCEngineMCMC & enginemcmc)
{
   Copy(enginemcmc);
   return *this;
}

// --------------------------------------------------------
unsigned int BCEngineMCMC::GetNFixedParameters()
{
   int n = 0;
   for (unsigned int i = 0; i < fParameters.Size(); ++i) {
      if (!fParameters[0]->Fixed())
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

   hprob->SetGlobalMode(fBestFitParameters.at(index));

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

   double gmode[] = { fBestFitParameters.at(i), fBestFitParameters.at(j) };
   hprob->SetGlobalMode(gmode);

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
std::vector<double> BCEngineMCMC::MCMCGetMaximumPoint(unsigned i) const
{
   // create a new vector with the length of fMCMCNParameters
   std::vector<double> x;

   // check if i is in range
   if (i >= fMCMCNChains)
      return x;

   // copy the point in the ith chain into the temporary vector
   for (unsigned j = 0; j < fParameters.Size(); ++j)
      x.push_back(fMCMCxMax.at(i * fParameters.Size() + j));

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
void BCEngineMCMC::MCMCSetInitialPositions(const std::vector<double> & x0s)
{
   // clear the existing initial position vector
   fMCMCInitialPosition.clear();

   // copy the initial positions
   unsigned n = x0s.size();

   for (unsigned i = 0; i < n; ++i)
      fMCMCInitialPosition.push_back(x0s.at(i));

   // use these initial positions for the Markov chain
   MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(std::vector< std::vector<double> > x0s)
{
   // create new vector
   std::vector<double> y0s;

   // loop over vector elements
   for (unsigned i = 0; i < x0s.size(); ++i)
      for (unsigned j = 0; j < x0s.at(i).size(); ++j)
         y0s.push_back((x0s.at(i)).at(j));

   MCMCSetInitialPositions(y0s);
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
         fMCMCTrees[i]->Branch(TString::Format("Parameter%i", j),
               &fMCMCx[i * fParameters.Size() + j],
               TString::Format("parameter %i/D", j));
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
double BCEngineMCMC::MCMCTrialFunctionSingle(unsigned ichain, unsigned iparameter)
{
   // no check of range for performance reasons

   // use uniform distribution
    //   return = fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter] * 2.0 * (0.5 - fRandom->Rndm());

   // Breit-Wigner width adjustable width
   return fMCMCThreadLocalStorage[ichain].rng->BreitWigner(0.0,
                                                               fMCMCTrialFunctionScaleFactor[ichain * fParameters.Size() + iparameter]);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(unsigned ichain) const
{
   // create a new vector with the length of fParameters
   std::vector<double> x;

   // check if ichain is in range
   if (ichain >= fMCMCNChains)
      return x;

   // copy the scale factors into the temporary vector
   for (unsigned j = 0; j < fParameters.Size(); ++j)
      x.push_back(fMCMCTrialFunctionScaleFactor.at(ichain * fParameters.Size() + j));

   return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(unsigned ichain, unsigned ipar)
{
   // check if ichain is in range
   if (ichain >= fMCMCNChains)
      return 0;

   // check if ipar is in range
   if (ipar >= fParameters.Size())
      return 0;

   // return component of ipar point in the ichain chain
   return fMCMCTrialFunctionScaleFactor.at(ichain *  fMCMCNChains + ipar);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetx(unsigned ichain)
{
   // create a new vector with the length of fParameters.Size()
   std::vector<double> x;

   // check if ichain is in range
   if (ichain >= fMCMCNChains)
      return x;

   // copy the point in the ichain chain into the temporary vector
   for (unsigned j = 0; j < fParameters.Size(); ++j)
      x.push_back(fMCMCx.at(ichain * fParameters.Size() + j));

   return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetx(unsigned ichain, unsigned ipar) const
{
   // check if ichain is in range
   if (ichain >= fMCMCNChains)
      return 0;

   // check if ipar is in range
   if (ipar >= fParameters.Size())
      return 0;

   // return component of jth point in the ith chain
   return fMCMCx.at(ichain *  fParameters.Size() + ipar);
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetLogProbx(unsigned ichain)
{
   // check if ichain is in range
   if (ichain >= fMCMCNChains)
      return -1;

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
       x[i] = fMCMCx[chain * fParameters.Size() + i] + x[i] * fParameters[i]->GetRangeWidth();

   // check if the point is in the correct volume.
   for (unsigned i = 0; i < fParameters.Size(); ++i)
       if (fParameters[i]->IsValid(x[i]))
         return false;

   return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned ichain, unsigned ipar, std::vector<double> &x)
{
   // get unscaled random point in the dimension of the chosen
   // parameter. this point might not be in the correct volume.
   double proposal = MCMCTrialFunctionSingle(ichain, ipar);

   // copy the old point into the new
   for (unsigned i = 0; i < fParameters.Size(); ++i)
      x[i] = fMCMCx[ichain * fParameters.Size() + i];

   // modify the parameter under study
   x[ipar] += proposal * fParameters[ipar]->GetRangeWidth();

   // check if the point is in the correct volume.
   if (fParameters[ipar]->IsValid(x[ipar]))
      return true;

   return false;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(unsigned chain, unsigned parameter)
{
   // calculate index
   unsigned index = chain * fParameters.Size();

   fMCMCCurrentChain = chain;

   // increase counter
   fMCMCNIterations[chain]++;

   // get proposal point
   if (!MCMCGetProposalPointMetropolis(chain, parameter, fMCMCThreadLocalStorage[chain].xLocal))
   {
      // increase counter
      fMCMCNTrialsFalse[chain * fParameters.Size() + parameter]++;

      // execute user code for every point
      MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, false);

      return false;
   }

   // calculate probabilities of the old and new points
   double p0 = fMCMCprob[chain];
   double p1 = LogEval(fMCMCThreadLocalStorage[chain].xLocal);

   // flag for accept
   bool accept = false;

   // if the new point is more probable, keep it ...
   if (p1 >= p0)
      accept = true;
   // ... or else throw dice.
   else
   {
      double r = log(fMCMCThreadLocalStorage[chain].rng->Rndm());

      if(r < p1 - p0)
         accept = true;
   }

   // fill the new point
   if(accept)
   {
      // increase counter
      fMCMCNTrialsTrue[chain * fParameters.Size() + parameter]++;

      // copy the point
      for(unsigned i = 0; i < fParameters.Size(); ++i)
      {
         // save the point
         fMCMCx[index + i] = fMCMCThreadLocalStorage[chain].xLocal[i];

         // save the probability of the point
         fMCMCprob[chain] = p1;
      }
   }
   else
   {
      // increase counter
      fMCMCNTrialsFalse[chain * fParameters.Size() + parameter]++;
   }

   // execute user code for every point
   MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, accept);

   return accept;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(unsigned chain)
{
   // calculate index
   unsigned index = chain * fParameters.Size();

   fMCMCCurrentChain = chain;

   // increase counter
   fMCMCNIterations[chain]++;

   // get proposal point
   if (!MCMCGetProposalPointMetropolis(chain, fMCMCThreadLocalStorage[chain].xLocal))
   {
      // increase counter
      for (unsigned i = 0; i < fParameters.Size(); ++i)
         fMCMCNTrialsFalse[chain * fParameters.Size() + i]++;

      // execute user code for every point
      MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, false);

      return false;
   }

   // calculate probabilities of the old and new points
   double p0 = fMCMCprob[chain];
   double p1 = LogEval(fMCMCThreadLocalStorage[chain].xLocal);

   // flag for accept
   bool accept = false;

   // if the new point is more probable, keep it ...
   if (p1 >= p0)
      accept = true;

   // ... or else throw dice.
   else
   {
      double r = log(fMCMCThreadLocalStorage[chain].rng->Rndm());

      if(r < p1 - p0)
         accept = true;
   }

   // fill the new point
   if(accept)
   {
      // increase counter
      for (unsigned i = 0; i < fParameters.Size(); ++i)
         fMCMCNTrialsTrue[chain * fParameters.Size() + i]++;

      // copy the point
      for(unsigned i = 0; i < fParameters.Size(); ++i)
      {
         // save the point
         fMCMCx[index + i] = fMCMCThreadLocalStorage[chain].xLocal[i];

         // save the probability of the point
         fMCMCprob[chain] = p1;
      }
   }
   else
   {
      // increase counter
      for (unsigned i = 0; i < fParameters.Size(); ++i)
         fMCMCNTrialsFalse[chain * fParameters.Size() + i]++;
   }

   // execute user code for every point
   MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, accept);

   return accept;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainCheckMaximum()
{
   // loop over all chains
   for (unsigned i = 0; i < fMCMCNChains; ++i)
   {
      // check if new maximum is found or chain is at the beginning
      if (fMCMCprob[i] > fMCMCprobMax[i] || fMCMCNIterations[i] == 1)
      {
         // copy maximum value
         fMCMCprobMax[i] = fMCMCprob[i];

         // copy mode of chain
         for (unsigned j = 0; j < fParameters.Size(); ++j)
            fMCMCxMax[i * fParameters.Size() + j] = fMCMCx[i * fParameters.Size() + j];
      }
   }
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainUpdateStatistics()
{
   // calculate number of entries in this part of the chain
   unsigned npoints = fMCMCNTrialsTrue[0] + fMCMCNTrialsFalse[0];

   // length of vectors
   unsigned nentries = fParameters.Size() * fMCMCNChains;

   // loop over all parameters of all chains
   for (unsigned i = 0; i < nentries; ++i) {
      // calculate mean value of each parameter in the chain for this part
      fMCMCxMean[i] += (fMCMCx[i] - fMCMCxMean[i]) / double(npoints);

      // calculate variance of each chain for this part
      if (npoints > 1)
         fMCMCxVar[i] = (1.0 - 1./double(npoints)) * fMCMCxVar[i]
            + (fMCMCx[i] - fMCMCxMean[i]) * (fMCMCx[i] - fMCMCxMean[i]) / double(npoints - 1);
   }

   // loop over chains
   for (unsigned i = 0; i < fMCMCNChains; ++i) {
      // calculate mean value of each chain for this part
      fMCMCprobMean[i] += (fMCMCprob[i] - fMCMCprobMean[i]) / double(npoints);

      // calculate variance of each chain for this part
      if (npoints > 1)
         fMCMCprobVar[i] = (1.0 - 1/double(npoints)) * fMCMCprobVar[i]
            + (fMCMCprob[i] - fMCMCprobMean[i]) * (fMCMCprob[i] - fMCMCprobMean[i]) / double(npoints - 1);
   }

}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainFillHistograms()
{
   // loop over chains
   for (unsigned i = 0; i < fMCMCNChains; ++i)
   {
      // fill each 1-dimensional histogram (if supposed to be filled)
      for (unsigned j = 0; j < fParameters.Size(); ++j)
         if (TH1 * h = fMCMCH1Marginalized[j])
            h->Fill(fMCMCx[i * fParameters.Size() + j]);

      // fill each 2-dimensional histogram (if supposed to be filled)
      unsigned counter = 0;

      for (unsigned j = 0; j < fParameters.Size(); ++j)
         for (unsigned k = j+1; k < fParameters.Size(); ++k)
         {
           if (TH2D * h = fMCMCH2Marginalized[counter])
             h->Fill(fMCMCx[i*fParameters.Size()+j],fMCMCx[i* fParameters.Size()+k]);
           counter ++;
         }
   }
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainTestConvergenceAllChains()
{
   // calculate number of entries in this part of the chain
   unsigned npoints = fMCMCNTrialsTrue[0] + fMCMCNTrialsFalse[0];

   if (fMCMCNChains > 1 && npoints > 1)
   {
      // define flag for convergence
      bool flag_convergence = true;

      // loop over parameters
      for (unsigned iparameters = 0; iparameters < fParameters.Size(); ++iparameters){
         // extract means and variances
         std::vector<double> means(fMCMCNChains);
         std::vector<double> variances(fMCMCNChains);

         // loop over chains
         for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains) {
            // get parameter index
            unsigned index = ichains * fParameters.Size() + iparameters;
            means[ichains] = fMCMCxMean[index];
            variances[ichains] = fMCMCxVar[index];
         }

         fMCMCRValueParameters[iparameters] = BCMath::Rvalue(means, variances, npoints, fMCMCRValueUseStrict);

         // set flag to false if convergence criterion is not fulfilled for the parameter
         if (! ((fMCMCRValueParameters[iparameters]-1.0) < fMCMCRValueParametersCriterion))
            flag_convergence = false;

         // else: leave convergence flag true for that parameter
      }

      fMCMCRValue = BCMath::Rvalue(fMCMCprobMean, fMCMCprobVar, npoints, fMCMCRValueUseStrict);

      // set flag to false if convergence criterion is not fulfilled for the log-likelihood
      if (!((fMCMCRValue - 1.0) < fMCMCRValueCriterion))
         flag_convergence = false;

      // remember number of iterations needed to converge
      if (fMCMCNIterationsConvergenceGlobal == -1 && flag_convergence == true)
         fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0] / fParameters.Size();
   }
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
int BCEngineMCMC::MCMCMetropolisPreRun()
{
   // print on screen
   BCLog::OutSummary("Pre-run Metropolis MCMC...");

   // initialize Markov chain
   MCMCInitialize();
   MCMCInitializeMarkovChains();

   // helper variable containing number of digits in the number of parameters
   int ndigits = (int)log10(fParameters.Size()) +1;
   if(ndigits<4)
      ndigits=4;

   // reset run statistics
   MCMCResetRunStatistics();
   fMCMCNIterationsConvergenceGlobal = -1;

   // perform run
   BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsMax));

   // don't write to file during pre run
   bool tempflag_writetofile = fMCMCFlagWriteChainToFile;
   fMCMCFlagWriteChainToFile = false;

   // initialize counter variables and flags
   fMCMCCurrentIteration = 1;   // counts the number of iterations
   unsigned counterupdate = 1;        // after how many iterations is an update needed?
   bool convergence = false;     // convergence reached?
   bool flagefficiency = false;  // efficiency reached?

   // array of efficiencies
   std::vector<double> efficiency;
   efficiency.assign(fParameters.Size() * fMCMCNChains, 0.0);

   // how often to check convergence and efficiencies?
   // it's either every fMCMCNParameters*nMCMCNIterationsUpdate (for 5 parameters the default would be 5000)
   // or it's fMCMCNIterationsUpdateMax (10000 by default)
   // whichever of the two is smaller
   unsigned updateLimit = ( fMCMCNIterationsUpdateMax<fMCMCNIterationsUpdate*(fParameters.Size())  && fMCMCNIterationsUpdateMax>0 ) ?
      fMCMCNIterationsUpdateMax : fMCMCNIterationsUpdate*(fParameters.Size());

   // loop over chains
   for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains) {
      // loop over parameters
      for (unsigned iparameter = 0; iparameter < fParameters.Size(); ++iparameter){
         // global index of the parameter (throughout all the chains)
         unsigned index = ichains * fParameters.Size() + iparameter;
         // reset counters
         fMCMCNTrialsTrue[index] = 0;
         fMCMCNTrialsFalse[index] = 0;
         fMCMCxMean[index] = fMCMCx[index];
         fMCMCxVar[index] = 0;
      }
      fMCMCprobMean[ichains] = fMCMCprob[ichains];
      fMCMCprobVar[ichains] = 0;
   }

   // set phase and cycle number
   fMCMCPhase = 1;

   // run chain ...
   // (a) for at least a minimum number of iterations,
   // (b) until a maximum number of iterations is reached,
   // (c) or until convergence is reached and the efficiency is in the
   //     specified region
   while (fMCMCCurrentIteration < int(fMCMCNIterationsPreRunMin) ||
          (fMCMCCurrentIteration < int(fMCMCNIterationsMax) && !(convergence && flagefficiency)))
   {
      //-------------------------------------------
      // reset flags and counters
      //-------------------------------------------

      // set convergence to false by default
      convergence = false;

      // set number of iterations needed to converge to negative
      fMCMCNIterationsConvergenceGlobal = -1;

      //-------------------------------------------
      // get new point in n-dim space
      //-------------------------------------------

      // if the flag is set then run over the parameters one after the other.
      if (fMCMCFlagOrderParameters)
      {
         // loop over parameters
         {
            for (unsigned iparameters = 0; iparameters < fParameters.Size(); ++iparameters)
            {
               // loop over chains

               unsigned chunk = 1; (void) chunk;
               unsigned ichains; (void) ichains;
#pragma omp parallel for shared(chunk) private(ichains)  schedule(static, chunk)
               for (ichains = 0; ichains < fMCMCNChains; ++ichains){
                  MCMCGetNewPointMetropolis(ichains, iparameters);
               }
               // search for global maximum
               MCMCInChainCheckMaximum();
            }
         }
      }

      // if the flag is not set then run over the parameters at the same time.
      else
      {
         // loop over chains
         {
            unsigned chunk = 1; (void) chunk;
            unsigned ichains;  (void) ichains;
#pragma omp parallel for shared(chunk) private(ichains)  schedule(static, chunk)
            for (ichains = 0; ichains < fMCMCNChains; ++ichains){
               MCMCGetNewPointMetropolis(ichains);
            }
         }
         // search for global maximum
         MCMCInChainCheckMaximum();
      }

      //-------------------------------------------
      // print out message to log
      //-------------------------------------------

      // progress printout
      if ( fMCMCCurrentIteration > 0 && fMCMCCurrentIteration % fMCMCNIterationsUpdate == 0 )
         BCLog::OutDetail(Form(" --> Iteration %i", fMCMCNIterations[0]/fParameters.Size()));

      //-------------------------------------------
      // update statistics
      //-------------------------------------------

      if (counterupdate > 1)
         MCMCInChainUpdateStatistics();

      //-------------------------------------------
      // update scale factors and check convergence
      //-------------------------------------------

      // debugKK
      // check if this line makes sense
      if ( counterupdate % updateLimit == 0 && counterupdate > 0 && fMCMCCurrentIteration >= int(fMCMCNIterationsPreRunMin))
      {
         // -----------------------------
         // reset flags and counters
         // -----------------------------

         bool rvalues_ok = true;

         static bool has_converged = false;

         // reset the number of iterations needed for convergence to
         // negative
         fMCMCNIterationsConvergenceGlobal = -1;

         // -----------------------------
         // check convergence status
         // -----------------------------

         // test convergence
         MCMCInChainTestConvergenceAllChains();

         // set convergence flag
         if (fMCMCNIterationsConvergenceGlobal > 0)
            convergence = true;

         // print convergence status:
         if (convergence && fMCMCNChains > 1)
            BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains converged within %i iterations.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
         else if (!convergence && fMCMCNChains > 1)
         {
            BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));

            BCLog::OutDetail("       - R-Values:");
            for (unsigned iparameter = 0; iparameter < fParameters.Size(); ++iparameter)
            {
               if(fabs(fMCMCRValueParameters[iparameter]-1.) < fMCMCRValueParametersCriterion)
                  BCLog::OutDetail(Form( TString::Format("         parameter %%%di :  %%.06f",ndigits), iparameter, fMCMCRValueParameters.at(iparameter)));
               else
               {
                  BCLog::OutDetail(Form( TString::Format("         parameter %%%di :  %%.06f <--",ndigits), iparameter, fMCMCRValueParameters.at(iparameter)));
                  rvalues_ok = false;
               }
            }
            if(fabs(fMCMCRValue-1.) < fMCMCRValueCriterion)
               BCLog::OutDetail(Form("         log-likelihood :  %.06f", fMCMCRValue));
            else
            {
               BCLog::OutDetail(Form("         log-likelihood :  %.06f <--", fMCMCRValue));
               rvalues_ok = false;
            }
         }

         // set convergence flag
         if(!has_converged)
            if(rvalues_ok)
               has_converged = true;

         // -----------------------------
         // check efficiency status
         // -----------------------------

         // set flag
         flagefficiency = true;

         bool flagprintefficiency = true;

         // loop over chains
         for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains)
         {
            // loop over parameters
            for (unsigned iparameter = 0; iparameter < fParameters.Size(); ++iparameter)
            {
               // global index of the parameter (throughout all the chains)
               unsigned index = ichains * fParameters.Size() + iparameter;

               // calculate efficiency
               efficiency[index] = double(fMCMCNTrialsTrue[index]) / double(fMCMCNTrialsTrue[index] + fMCMCNTrialsFalse[index]);

               // adjust scale factors if efficiency is too low
               if (efficiency[index] < fMCMCEfficiencyMin && fMCMCTrialFunctionScaleFactor[index] > .01)
               {
                  if (flagprintefficiency)
                  {
                     BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range."));
                     BCLog::OutDetail(Form("       - Efficiencies:"));
                     flagprintefficiency = false;
                  }

                  double fscale=2.;
                  if(has_converged && fMCMCEfficiencyMin/efficiency[index] > 2.)
                     fscale = 4.;
                  fMCMCTrialFunctionScaleFactor[index] /= fscale;

                  BCLog::OutDetail(Form("         Efficiency of parameter %i dropped below %.2f%% (eps = %.2f%%) in chain %i. Set scale to %.4g",
                                        iparameter, 100. * fMCMCEfficiencyMin, 100. * efficiency[index], ichains, fMCMCTrialFunctionScaleFactor[index]));
               }

               // adjust scale factors if efficiency is too high
               else if (efficiency[index] > fMCMCEfficiencyMax && fMCMCTrialFunctionScaleFactor[index] < 1.0)
               {
                  if (flagprintefficiency)
                  {
                     BCLog::OutDetail(Form("   * Efficiency status: Efficiencies not within pre-defined ranges."));
                     BCLog::OutDetail(Form("     - Efficiencies:"));
                     flagprintefficiency = false;
                  }

                  fMCMCTrialFunctionScaleFactor[index] *= 2.;

                  BCLog::OutDetail(Form("         Efficiency of parameter %i above %.2f%% (eps = %.2f%%) in chain %i. Set scale to %.4g",
                                        iparameter, 100.0 * fMCMCEfficiencyMax, 100.0 * efficiency[index], ichains, fMCMCTrialFunctionScaleFactor[index]));
               }

               // check flag
               if ((efficiency[index] < fMCMCEfficiencyMin && fMCMCTrialFunctionScaleFactor[index] > .01)
                   || (efficiency[index] > fMCMCEfficiencyMax && fMCMCTrialFunctionScaleFactor[index] < 1.))
                  flagefficiency = false;
            } // end of running over all parameters
         } // end of running over all chains

         // print to screen
         if (flagefficiency)
            BCLog::OutDetail(Form("     * Efficiency status: Efficiencies within pre-defined ranges."));

         // -----------------------------
         // reset counters
         // -----------------------------

         counterupdate = 0;

         // loop over chains
         for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains) {
            // loop over parameters
            for (unsigned iparameter = 0; iparameter < fParameters.Size(); ++iparameter){
               // global index of the parameter (throughout all the chains)
               unsigned index = ichains * fParameters.Size() + iparameter;
               // reset counters
               fMCMCNTrialsTrue[index] = 0;
               fMCMCNTrialsFalse[index] = 0;
               fMCMCxMean[index] = fMCMCx[index];
               fMCMCxVar[index] = 0;
            }
            fMCMCprobMean[ichains] = fMCMCprob[ichains];
            fMCMCprobVar[ichains] = 0;
         }
      } // end if update scale factors and check convergence

      //-------------------------------------------
      // write chain to file
      //-------------------------------------------

      // write chain to file
      if (fMCMCFlagWritePreRunToFile)
         MCMCInChainWriteChains();

      //-------------------------------------------
      // increase counters
      //-------------------------------------------
      fMCMCCurrentIteration++;
      counterupdate++;

   } // end of running

   // did we check convergence at least once ?
   if (fMCMCCurrentIteration < int(updateLimit))
   {
      BCLog::OutWarning(" Convergence never checked !");
      BCLog::OutWarning("   Increase maximum number of iterations in the pre-run /MCMCSetNIterationsMax()/");
      BCLog::OutWarning("   or decrease maximum number of iterations for update  /MCMCSetNIterationsUpdateMax()/");
   }

   // ---------------
   // after chain run
   // ---------------

   // define convergence status
   if (fMCMCNIterationsConvergenceGlobal > 0)
      fMCMCFlagConvergenceGlobal = true;
   else
      fMCMCFlagConvergenceGlobal = false;

   // print convergence status
   if (fMCMCFlagConvergenceGlobal && fMCMCNChains > 1 && !flagefficiency)
      BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations but could not adjust scales.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));

   else if (fMCMCFlagConvergenceGlobal && fMCMCNChains > 1 && flagefficiency)
      BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations and all scales are adjusted.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));

   else if (!fMCMCFlagConvergenceGlobal && (fMCMCNChains > 1) && flagefficiency)
      BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations.", fMCMCNChains, fMCMCNIterationsMax));

   else if (!fMCMCFlagConvergenceGlobal && (fMCMCNChains > 1) && !flagefficiency)
      BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations and could not adjust scales.", fMCMCNChains, fMCMCNIterationsMax));

   else if(fMCMCNChains == 1)
      BCLog::OutSummary(" --> No convergence criterion for a single chain defined.");

   else
      BCLog::OutSummary(" --> Only one Markov chain. No global convergence criterion defined.");

   BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCCurrentIteration));


   // print efficiencies
   std::vector<double> efficiencies;

   for (unsigned i = 0; i < fParameters.Size(); ++i)
      efficiencies.push_back(0.);

   BCLog::OutDetail(" --> Average efficiencies:");
   for (unsigned i = 0; i < fParameters.Size(); ++i)
   {
      for (unsigned j = 0; j < fMCMCNChains; ++j)
         efficiencies[i] += efficiency[j * fParameters.Size() + i] / double(fMCMCNChains);

      BCLog::OutDetail(Form(TString::Format(" -->      parameter %%%di :  %%.02f%%%%",ndigits), i, 100. * efficiencies[i]));
   }


   // print scale factors
   std::vector<double> scalefactors;

   for (unsigned i = 0; i < fParameters.Size(); ++i)
      scalefactors.push_back(0.0);

   BCLog::OutDetail(" --> Average scale factors:");
   for (unsigned i = 0; i < fParameters.Size(); ++i)
   {
      for (unsigned j = 0; j < fMCMCNChains; ++j)
         scalefactors[i] += fMCMCTrialFunctionScaleFactor[j * fParameters.Size() + i] / double(fMCMCNChains);

      BCLog::OutDetail(Form( TString::Format(" -->      parameter %%%di :  %%.02f%%%%",ndigits), i, 100. * scalefactors[i]));
   }

   // reset flag
   fMCMCFlagWriteChainToFile = tempflag_writetofile;

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
   // check if prerun should be performed
   if (fMCMCFlagPreRun)
      MCMCMetropolisPreRun();
   else
      BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Not running prerun. This can cause trouble if the data have changed.");

   // print to screen
   BCLog::OutSummary( "Run Metropolis MCMC...");

   // reset run statistics
   MCMCResetRunStatistics();

   // set phase and cycle number
   fMCMCPhase = 2;

   // perform run
   BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

   int nwrite = fMCMCNIterationsRun/10;
   if(nwrite < 100)
      nwrite=100;
   else if(nwrite < 500)
      nwrite=1000;
   else if(nwrite < 10000)
      nwrite=1000;
   else
      nwrite=10000;

   // start the run
   for (fMCMCCurrentIteration = 1; fMCMCCurrentIteration <= int(fMCMCNIterationsRun); ++fMCMCCurrentIteration)
   {
      if ( (fMCMCCurrentIteration)%nwrite == 0 )
         BCLog::OutDetail(Form(" --> iteration number %i (%.2f%%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration)/(double)fMCMCNIterationsRun*100.));

      // if the flag is set then run over the parameters one after the other.
      if (fMCMCFlagOrderParameters)
      {
         // loop over parameters
         for (unsigned iparameters = 0; iparameters < fParameters.Size(); ++iparameters)
         {
            // loop over chains
            {
               unsigned chunk = 1; (void) chunk;
               unsigned ichains; (void) ichains;
#pragma omp parallel for shared(chunk) private(ichains)  schedule(static, chunk)
               for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains)
               {
                  MCMCGetNewPointMetropolis(ichains, iparameters);
               }
            }
            // reset current chain
            fMCMCCurrentChain = -1;

            // update search for maximum
            MCMCInChainCheckMaximum();

         } // end loop over all parameters

         // check if the current iteration is consistent with the lag
         if ( fMCMCCurrentIteration % fMCMCNLag == 0)
         {
            // do anything interface
            MCMCIterationInterface();

            // fill histograms
            if ( ! fMCMCH1Marginalized.empty() or ! fMCMCH2Marginalized.empty())
               MCMCInChainFillHistograms();

            // write chain to file
            if (fMCMCFlagWriteChainToFile)
               MCMCInChainWriteChains();
         }
      }
      // if the flag is not set then run over the parameters at the same time.
      else
      {
         // loop over chains
         {
            unsigned chunk = 1; (void) chunk;
            unsigned ichains; (void) ichains;
#pragma omp parallel for shared(chunk) private(ichains)  schedule(static, chunk)
            for (unsigned ichains = 0; ichains < fMCMCNChains; ++ichains)
            {
               // get new point
               MCMCGetNewPointMetropolis(ichains);
            }
         }
         // reset current chain
         fMCMCCurrentChain = -1;

         // update search for maximum
         MCMCInChainCheckMaximum();

         // check if the current iteration is consistent with the lag
         if (fMCMCCurrentIteration % fMCMCNLag == 0)
         {
            // do anything interface
            MCMCIterationInterface();

            // fill histograms
            if ( ! fMCMCH1Marginalized.empty() or ! fMCMCH2Marginalized.empty())
               MCMCInChainFillHistograms();

            // write chain to file
            if (fMCMCFlagWriteChainToFile)
               MCMCInChainWriteChains();
         }
      }

   } // end run

   // print convergence status
   BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

   // print modes

   // find global maximum
   double probmax = fMCMCprobMax.at(0);
   unsigned probmaxindex = 0;

   // loop over all chains and find the maximum point
   for (unsigned i = 1; i < fMCMCNChains; ++i) {
      if (fMCMCprobMax.at(i) > probmax)
      {
         probmax = fMCMCprobMax.at(i);
         probmaxindex = i;
      }
   }

   // save if improved the log posterior
   if (fBestFitParameters.empty() || probmax > fLogMaximum) {
      fLogMaximum = probmax;
      fBestFitParameters.assign(fParameters.Size(), 0.0);
      for (unsigned i = 0; i < fParameters.Size(); ++i)
         fBestFitParameters[i] = fMCMCxMax[probmaxindex * fParameters.Size() + i];

      fOptimizationMethodMode = BCEngineMCMC::kOptMetropolis;
   }

   BCLog::OutDetail(" --> Global mode from MCMC:");
   BCLog::OutDebug(Form(" --> Posterior value: %g", probmax));
   int ndigits = (int) log10(fParameters.Size());
   for (unsigned i = 0; i < fParameters.Size(); ++i)
      BCLog::OutDetail(Form( TString::Format(" -->      parameter %%%di:   %%.4g", ndigits+1),
            i, fBestFitParameters[i]));

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
   for (unsigned j = 0; j < fMCMCNChains; ++j)
   {
      fMCMCNIterations[j]  = 0;
      fMCMCNTrialsTrue[j]  = 0;
      fMCMCNTrialsFalse[j] = 0;
      fMCMCprobMean[j]         = 0;
      fMCMCprobVar[j]     = 0;

      for (unsigned k = 0; k < fParameters.Size(); ++k)
      {
         fMCMCNTrialsTrue[j * fParameters.Size() + k]  = 0;
         fMCMCNTrialsFalse[j * fParameters.Size() + k] = 0;
      }
   }

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
   std::vector<double> x0;

   for (unsigned j = 0; j < fMCMCNChains; ++j)
   {
      x0.clear();
      for (unsigned i = 0; i < fParameters.Size(); ++i)
         x0.push_back(fMCMCx[j * fParameters.Size() + i]);
      fMCMCprob[j] = LogEval(x0);
   }

   x0.clear();
}

// --------------------------------------------------------
void BCEngineMCMC::ResetResults()
{
   // reset variables
   fMCMCNIterations.clear();
   fMCMCNTrialsTrue.clear();
   fMCMCNTrialsFalse.clear();
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

   fBestFitParameters.clear();
   fLogMaximum = -std::numeric_limits<double>::max();
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

   fMCMCNTrialsTrue.assign(fMCMCNChains * fParameters.Size(), 0);
   fMCMCNTrialsFalse.assign(fMCMCNChains * fParameters.Size(), 0);
   fMCMCxMax.assign(fMCMCNChains * fParameters.Size(), 0.);
   fMCMCxMean.assign(fMCMCNChains * fParameters.Size(), 0);
   fMCMCxVar.assign(fMCMCNChains * fParameters.Size(), 0);

   fMCMCRValueParameters.assign(fParameters.Size(), 0.);

   SyncThreadStorage();

   if (fMCMCTrialFunctionScaleFactorStart.size() == 0)
      fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains * fParameters.Size(), 1.0);
   else
      for (unsigned i = 0; i < fMCMCNChains; ++i)
         for (unsigned j = 0; j < fParameters.Size(); ++j)
            fMCMCTrialFunctionScaleFactor.push_back(fMCMCTrialFunctionScaleFactorStart.at(j));

   // set initial position
   if (fMCMCFlagInitialPosition == 2) // user defined points
   {
      // define flag
      bool flag = true;

      // check the length of the array of initial positions
      if (fMCMCInitialPosition.size() != (fMCMCNChains * fParameters.Size()))
      {
         BCLog::OutError("BCEngine::MCMCInitialize : Length of vector containing initial positions does not have required length.");
         flag = false;
      }

      // check the boundaries
      if (flag)
      {
         for (unsigned j = 0; j < fMCMCNChains; ++j)
            for (unsigned i = 0; i < fParameters.Size(); ++i)
               if (fParameters[i]->IsValid(fMCMCInitialPosition[j * fParameters.Size() + i]))
               {
                  BCLog::OutError("BCEngine::MCMCInitialize : Initial position out of boundaries.");
                  flag = false;
               }
      }

      // check flag
      if (!flag)
         fMCMCFlagInitialPosition = 1;
   }

   if (fMCMCFlagInitialPosition == 0) // center of the interval
      for (unsigned j = 0; j < fMCMCNChains; ++j)
         for (unsigned i = 0; i < fParameters.Size(); ++i)
            fMCMCx.push_back(fParameters[i]->GetLowerLimit() + .5 * fParameters[i]->GetRangeWidth());

   else if (fMCMCFlagInitialPosition == 2) // user defined
   {
      for (unsigned j = 0; j < fMCMCNChains; ++j)
         for (unsigned i = 0; i < fParameters.Size(); ++i)
            fMCMCx.push_back(fMCMCInitialPosition.at(j * fParameters.Size() + i));
   }

   else
   {
      for (unsigned j = 0; j < fMCMCNChains; ++j) // random number (default)
         for (unsigned i = 0; i < fParameters.Size(); ++i)
            fMCMCx.push_back(fParameters[i]->GetLowerLimit() + fMCMCThreadLocalStorage[j].rng->Rndm() * fParameters[i]->GetRangeWidth());
   }

   // copy the point of the first chain
   std::copy(fMCMCx.begin(), fMCMCx.begin() + fParameters.Size(), fMCMCThreadLocalStorage.at(0).xLocal.begin());

   // define 1-dimensional histograms for marginalization
   bool fillAny=false;
   for(unsigned i = 0; i < fParameters.Size(); ++i)
   {
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
   if (!fillAny) {
      fMCMCH1Marginalized.clear();
   }
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
   // what's within this method will be executed
   // for every iteration of the MCMC

   // fill error band
   MCMCFillErrorBand();

   // do user defined stuff
   MCMCUserIterationInterface();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCFillErrorBand()
{
   if (!fFillErrorBand)
      return;

   // function fitting
   if (fFitFunctionIndexX < 0)
      return;

   // loop over all possible x values ...
   if (fErrorBandContinuous) {
      double x = 0;
      for (unsigned ix = 0; ix < fErrorBandNbinsX; ix++) {
         // calculate x
         x = fErrorBandXY->GetXaxis()->GetBinCenter(ix + 1);

         // calculate y
         std::vector<double> xvec;
         xvec.push_back(x);

         // loop over all chains
         for (unsigned ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
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
      unsigned ndatapoints = fErrorBandX.size();
      double x = 0;

      for (unsigned ix = 0; ix < ndatapoints; ++ix) {
         // calculate x
         x = fErrorBandX.at(ix);

         // calculate y
         std::vector<double> xvec;
         xvec.push_back(x);

         // loop over all chains
         for (unsigned ichain = 0; ichain < MCMCGetNChains(); ++ichain) {
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
