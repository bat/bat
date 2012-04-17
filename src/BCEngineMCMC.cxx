/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "BCEngineMCMC.h"

#include "BCLog.h"
#include "BCMath.h"
#include "BCParameter.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TRandom3.h>

#include <math.h>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC()
{
   // set default parameters for the mcmc
   MCMCSetValuesDefault();

   // initialize random number generator
   fRandom = new TRandom3(0);
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(int n)
{
   // set number of chains to n
   fMCMCNChains = n;

   // call default constructor
   BCEngineMCMC();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetValuesDefault()
{
   fMCMCNParameters          = 0;
   fMCMCFlagWriteChainToFile = false;
   fMCMCFlagWritePreRunToFile = false;
   fMCMCFlagPreRun           = true;
   fMCMCFlagRun              = false;
   fMCMCFlagFillHistograms   = true;
   fMCMCEfficiencyMin        = 0.15;
   fMCMCEfficiencyMax        = 0.50;
   fMCMCFlagInitialPosition  = 1;
   fMCMCNLag                 = 1;
   fMCMCCurrentIteration     = -1;
   fMCMCCurrentChain         = -1;

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
   if (fRandom)
      delete fRandom;

   // delete 1-d marginalized distributions
   for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
      if (fMCMCH1Marginalized[i])
         delete fMCMCH1Marginalized[i];
   fMCMCH1Marginalized.clear();

   // delete 2-d marginalized distributions
   for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
      if (fMCMCH2Marginalized[i])
         delete fMCMCH2Marginalized[i];
   fMCMCH2Marginalized.clear();
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC & enginemcmc)
{
   fMCMCPointerToGetProposalPoint = enginemcmc.fMCMCPointerToGetProposalPoint;
   fMCMCNParameters           = enginemcmc.fMCMCNParameters;
   fMCMCBoundaryMin           = enginemcmc.fMCMCBoundaryMin;
   fMCMCBoundaryMax           = enginemcmc.fMCMCBoundaryMax;
   fMCMCFlagsFillHistograms   = enginemcmc.fMCMCFlagsFillHistograms;
   fMCMCNChains               = enginemcmc.fMCMCNChains;
   fMCMCNLag                  = enginemcmc.fMCMCNLag;
   fMCMCNIterations           = enginemcmc.fMCMCNIterations;
   fMCMCCurrentIteration      = enginemcmc.fMCMCCurrentIteration;
   fMCMCCurrentChain          = enginemcmc.fMCMCCurrentChain;
   fMCMCNIterationsUpdate     = enginemcmc.fMCMCNIterationsUpdate;
   fMCMCNIterationsUpdateMax  = enginemcmc.fMCMCNIterationsUpdateMax;
   fMCMCNIterationsConvergenceGlobal = enginemcmc.fMCMCNIterationsConvergenceGlobal;
   fMCMCFlagConvergenceGlobal = enginemcmc.fMCMCFlagConvergenceGlobal;
   fMCMCNIterationsMax        = enginemcmc.fMCMCNIterationsMax;
   fMCMCNIterationsRun        = enginemcmc.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin  = enginemcmc.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue           = enginemcmc.fMCMCNTrialsTrue;
   fMCMCNTrialsFalse          = enginemcmc.fMCMCNTrialsFalse;
   fMCMCFlagWriteChainToFile  = enginemcmc.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile = enginemcmc.fMCMCFlagWritePreRunToFile;
   fMCMCTrialFunctionScaleFactor = enginemcmc.fMCMCTrialFunctionScaleFactor;
   fMCMCTrialFunctionScaleFactorStart = enginemcmc.fMCMCTrialFunctionScaleFactorStart;
   fMCMCFlagPreRun            = enginemcmc.fMCMCFlagPreRun;
   fMCMCFlagRun               = enginemcmc.fMCMCFlagRun;
   fMCMCInitialPosition       = enginemcmc.fMCMCInitialPosition;
   fMCMCEfficiencyMin         = enginemcmc.fMCMCEfficiencyMin;
   fMCMCEfficiencyMax         = enginemcmc.fMCMCEfficiencyMax;
   fMCMCFlagInitialPosition   = enginemcmc.fMCMCFlagInitialPosition;
   fMCMCFlagOrderParameters   = enginemcmc.fMCMCFlagOrderParameters;
   fMCMCFlagFillHistograms    = enginemcmc.fMCMCFlagFillHistograms;
   fMCMCPhase                 = enginemcmc.fMCMCPhase;
   fMCMCCycle                 = enginemcmc.fMCMCCycle;
   fMCMCx                     = enginemcmc.fMCMCx;
   fMCMCxMax                  = enginemcmc.fMCMCxMax;
   fMCMCxMean                 = enginemcmc.fMCMCxMean;
   fMCMCxVar                  = enginemcmc.fMCMCxVar;
   fMCMCxLocal                = enginemcmc.fMCMCxLocal;
   fMCMCprob                  = enginemcmc.fMCMCprob;
   fMCMCprobMax               = enginemcmc.fMCMCprobMax;
   fMCMCprobMean              = enginemcmc.fMCMCprobMean;
   fMCMCprobVar               = enginemcmc.fMCMCprobVar;
   fMCMCRValueUseStrict       = enginemcmc.fMCMCRValueUseStrict;
   fMCMCRValueCriterion       = enginemcmc.fMCMCRValueCriterion ;
   fMCMCRValueParametersCriterion = enginemcmc.fMCMCRValueParametersCriterion;
   fMCMCRValue                = enginemcmc.fMCMCRValue;
   fMCMCRValueParameters      = enginemcmc.fMCMCRValueParameters;
   if (enginemcmc.fRandom)
      fRandom = new TRandom3(*enginemcmc.fRandom);
   else
      fRandom = 0;
   fMCMCH1NBins               = enginemcmc.fMCMCH1NBins;

   for (int i = 0; i < int(enginemcmc.fMCMCH1Marginalized.size()); ++i) {
      if (enginemcmc.fMCMCH1Marginalized.at(i))
         fMCMCH1Marginalized.push_back(new TH1D(*(enginemcmc.fMCMCH1Marginalized.at(i))));
      else
         fMCMCH1Marginalized.push_back(0);
   }

   for (int i = 0; i < int(enginemcmc.fMCMCH2Marginalized.size()); ++i) {
      if (enginemcmc.fMCMCH2Marginalized.at(i))
         fMCMCH2Marginalized.push_back(new TH2D(*(enginemcmc.fMCMCH2Marginalized.at(i))));
      else
         fMCMCH2Marginalized.push_back(0);
   }

   for (int i = 0; i < int(enginemcmc.fMCMCTrees.size()); ++i) {
      // debugKK
      fMCMCTrees.push_back(0);
   }
}

// ---------------------------------------------------------
BCEngineMCMC & BCEngineMCMC::operator = (const BCEngineMCMC & enginemcmc)
{
   fMCMCPointerToGetProposalPoint = enginemcmc.fMCMCPointerToGetProposalPoint;
   fMCMCNParameters           = enginemcmc.fMCMCNParameters;
   fMCMCBoundaryMin           = enginemcmc.fMCMCBoundaryMin;
   fMCMCBoundaryMax           = enginemcmc.fMCMCBoundaryMax;
   fMCMCFlagsFillHistograms   = enginemcmc.fMCMCFlagsFillHistograms;
   fMCMCNChains               = enginemcmc.fMCMCNChains;
   fMCMCNLag                  = enginemcmc.fMCMCNLag;
   fMCMCNIterations           = enginemcmc.fMCMCNIterations;
   fMCMCCurrentIteration      = enginemcmc.fMCMCCurrentIteration;
   fMCMCCurrentChain          = enginemcmc.fMCMCCurrentChain;
   fMCMCNIterationsUpdate     = enginemcmc.fMCMCNIterationsUpdate;
   fMCMCNIterationsUpdateMax  = enginemcmc.fMCMCNIterationsUpdateMax;
   fMCMCNIterationsConvergenceGlobal = enginemcmc.fMCMCNIterationsConvergenceGlobal;
   fMCMCFlagConvergenceGlobal = enginemcmc.fMCMCFlagConvergenceGlobal;
   fMCMCNIterationsMax        = enginemcmc.fMCMCNIterationsMax;
   fMCMCNIterationsRun        = enginemcmc.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin  = enginemcmc.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue           = enginemcmc.fMCMCNTrialsTrue;
   fMCMCNTrialsFalse          = enginemcmc.fMCMCNTrialsFalse;
   fMCMCFlagWriteChainToFile  = enginemcmc.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile = enginemcmc.fMCMCFlagWritePreRunToFile;
   fMCMCTrialFunctionScaleFactor = enginemcmc.fMCMCTrialFunctionScaleFactor;
   fMCMCTrialFunctionScaleFactorStart = enginemcmc.fMCMCTrialFunctionScaleFactorStart;
   fMCMCFlagPreRun            = enginemcmc.fMCMCFlagPreRun;
   fMCMCFlagRun               = enginemcmc.fMCMCFlagRun;
   fMCMCInitialPosition       = enginemcmc.fMCMCInitialPosition;
   fMCMCEfficiencyMin         = enginemcmc.fMCMCEfficiencyMin;
   fMCMCEfficiencyMax         = enginemcmc.fMCMCEfficiencyMax;
   fMCMCFlagInitialPosition   = enginemcmc.fMCMCFlagInitialPosition;
   fMCMCFlagOrderParameters   = enginemcmc.fMCMCFlagOrderParameters;
   fMCMCFlagFillHistograms    = enginemcmc.fMCMCFlagFillHistograms;
   fMCMCPhase                 = enginemcmc.fMCMCPhase;
   fMCMCCycle                 = enginemcmc.fMCMCCycle;
   fMCMCx                     = enginemcmc.fMCMCx;
   fMCMCxMax                  = enginemcmc.fMCMCxMax;
   fMCMCxMean                 = enginemcmc.fMCMCxMean;
   fMCMCxVar                  = enginemcmc.fMCMCxVar;
   fMCMCxLocal                = enginemcmc.fMCMCxLocal;
   fMCMCprob                  = enginemcmc.fMCMCprob;
   fMCMCprobMax               = enginemcmc.fMCMCprobMax;
   fMCMCprobMean              = enginemcmc.fMCMCprobMean;
   fMCMCprobVar               = enginemcmc.fMCMCprobVar;
   fMCMCRValueUseStrict       = enginemcmc.fMCMCRValueUseStrict;
   fMCMCRValueCriterion       = enginemcmc.fMCMCRValueCriterion ;
   fMCMCRValueParametersCriterion = enginemcmc.fMCMCRValueParametersCriterion;
   fMCMCRValue                = enginemcmc.fMCMCRValue;
   fMCMCRValueParameters      = enginemcmc.fMCMCRValueParameters;
   if (enginemcmc.fRandom)
      fRandom = new TRandom3(*enginemcmc.fRandom);
   else
      fRandom = 0;
   fMCMCH1NBins               = enginemcmc.fMCMCH1NBins;

   for (int i = 0; i < int(enginemcmc.fMCMCH1Marginalized.size()); ++i) {
      if (enginemcmc.fMCMCH1Marginalized.at(i))
         fMCMCH1Marginalized.push_back(new TH1D(*(enginemcmc.fMCMCH1Marginalized.at(i))));
      else
         fMCMCH1Marginalized.push_back(0);
   }

   for (int i = 0; i < int(enginemcmc.fMCMCH2Marginalized.size()); ++i) {
      if (enginemcmc.fMCMCH2Marginalized.at(i))
         fMCMCH2Marginalized.push_back(new TH2D(*(enginemcmc.fMCMCH2Marginalized.at(i))));
      else
         fMCMCH2Marginalized.push_back(0);
   }

   for (int i = 0; i < int(enginemcmc.fMCMCTrees.size()); ++i) {
      // debugKK
      fMCMCTrees.push_back(0);
   }

   // return this
   return *this;
}


// --------------------------------------------------------
TH1D * BCEngineMCMC::MCMCGetH1Marginalized(int index)
{
      // check index
   if (index < 0 || index >= int(fMCMCH1Marginalized.size()))
   {
      BCLog::OutWarning("BCEngineMCMC::MCMCGetH1Marginalized. Index out of range.");
      return 0;
   }

   return fMCMCH1Marginalized[index];
}

// --------------------------------------------------------
TH2D * BCEngineMCMC::MCMCGetH2Marginalized(int index1, int index2)
{
   int counter = 0;
   int index = 0;

   // search for correct combination
   for(int i = 0; i < fMCMCNParameters; i++)
      for (int j = 0; j < i; ++j)
      {
         if(j == index1 && i == index2)
            index = counter;
         counter++;
      }

   // check index
   if (index < 0 || index >= int(fMCMCH2Marginalized.size()))
   {
      BCLog::OutWarning("BCEngineMCMC::MCMCGetH2Marginalized. Index out of range.");
      return 0;
   }

   return fMCMCH2Marginalized[index];
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetMaximumPoint(int i)
{
   // create a new vector with the lenght of fMCMCNParameters
   std::vector<double> x;

   // check if i is in range
   if (i < 0 || i >= fMCMCNChains)
      return x;

   // copy the point in the ith chain into the temporary vector
   for (int j = 0; j < fMCMCNParameters; ++j)
   x.push_back(fMCMCxMax.at(i * fMCMCNParameters + j));

   return x;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetNChains(int n)
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
   int n = int(x0s.size());

   for (int i = 0; i < n; ++i)
      fMCMCInitialPosition.push_back(x0s.at(i));

   // use these intial positions for the Markov chain
   MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(std::vector< std::vector<double> > x0s)
{
   // create new vector
   std::vector<double> y0s;

   // loop over vector elements
   for (int i = 0; i < int(x0s.size()); ++i)
      for (int j = 0; j < int((x0s.at(i)).size()); ++j)
         y0s.push_back((x0s.at(i)).at(j));

   MCMCSetInitialPositions(y0s);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(bool flag)
{
   fMCMCFlagFillHistograms = flag;

   for (int i = 0; i < fMCMCNParameters; ++i)
      fMCMCFlagsFillHistograms[i] = flag;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(int index, bool flag)
{
   // check if index is within range
   if (index < 0 || index > fMCMCNParameters)
   {
      BCLog::OutWarning("BCEngineMCMC :MCMCSetFlagFillHistograms. Index out of range.");
      return;
   }

   // set flag
   fMCMCFlagsFillHistograms[index] = flag;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetMarkovChainTrees(std::vector<TTree *> trees)
{
// clear vector
   fMCMCTrees.clear();

   // copy tree
   for (int i = 0; i < int(trees.size()); ++i)
      fMCMCTrees.push_back(trees[i]);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInitializeMarkovChainTrees()
{
// clear vector
   fMCMCTrees.clear();

// create new trees
   for (int i = 0; i < fMCMCNChains; ++i) {
      fMCMCTrees.push_back(new TTree(TString::Format("MarkovChainTree_%i", i), "MarkovChainTree"));
      fMCMCTrees[i]->Branch("Iteration",       &fMCMCNIterations[i],  "iteration/I");
      fMCMCTrees[i]->Branch("NParameters",     &fMCMCNParameters,     "parameters/I");
      fMCMCTrees[i]->Branch("LogProbability",  &fMCMCprob[i],         "log(probability)/D");
      fMCMCTrees[i]->Branch("Phase",           &fMCMCPhase,           "phase/I");
      fMCMCTrees[i]->Branch("Cycle",           &fMCMCCycle,           "cycle/I");

      for (int j = 0; j < fMCMCNParameters; ++j)
         fMCMCTrees[i]->Branch(TString::Format("Parameter%i", j),
               &fMCMCx[i * fMCMCNParameters + j],
               TString::Format("parameter %i/D", j));
   }
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunction(int ichain, std::vector<double> &x)
{
   // call MCMCTrialFunctionSingle() for all parameters by default
   for (int i = 0; i < fMCMCNParameters; ++i)
      x[i] = MCMCTrialFunctionSingle(ichain, i);
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCTrialFunctionSingle(int ichain, int iparameter)
{
   // no check of range for performance reasons

   // use uniform distribution
//   return = fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter] * 2.0 * (0.5 - fRandom->Rndm());

   // Breit-Wigner width adjustable width
   return fRandom->BreitWigner(0.0, fMCMCTrialFunctionScaleFactor[ichain * fMCMCNParameters + iparameter]);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int ichain)
{
   // create a new vector with the length of fMCMCNParameters
   std::vector<double> x;

   // check if ichain is in range
   if (ichain < 0 || ichain >= fMCMCNChains)
      return x;

   // copy the scale factors into the temporary vector
   for (int j = 0; j < fMCMCNParameters; ++j)
      x.push_back(fMCMCTrialFunctionScaleFactor.at(ichain * fMCMCNParameters + j));

   return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetTrialFunctionScaleFactor(int ichain, int ipar)
{
   // check if ichain is in range
   if (ichain < 0 || ichain >= fMCMCNChains)
      return 0;

   // check if ipar is in range
   if (ipar < 0 || ipar >= fMCMCNParameters)
      return 0;

   // return component of ipar point in the ichain chain
   return fMCMCTrialFunctionScaleFactor.at(ichain *  fMCMCNChains + ipar);
}

// --------------------------------------------------------
std::vector<double> BCEngineMCMC::MCMCGetx(int ichain)
{
   // create a new vector with the length of fMCMCNParameters
   std::vector<double> x;

   // check if ichain is in range
   if (ichain < 0 || ichain >= fMCMCNChains)
      return x;

   // copy the point in the ichain chain into the temporary vector
   for (int j = 0; j < fMCMCNParameters; ++j)
      x.push_back(fMCMCx.at(ichain * fMCMCNParameters + j));

   return x;
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetx(int ichain, int ipar)
{
   // check if ichain is in range
   if (ichain < 0 || ichain >= fMCMCNChains)
      return 0;

   // check if ipar is in range
   if (ipar < 0 || ipar >= fMCMCNParameters)
      return 0;

   // return component of jth point in the ith chain
   return fMCMCx.at(ichain *  fMCMCNParameters + ipar);
}

// --------------------------------------------------------
double BCEngineMCMC::MCMCGetLogProbx(int ichain)
{
   // check if ichain is in range
   if (ichain < 0 || ichain >= fMCMCNChains)
      return -1;

   // return log of the probability at the current point in the ith chain
   return fMCMCprob.at(ichain);
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int chain, std::vector<double> &x)
{
   // get unscaled random point. this point might not be in the correct volume.
   MCMCTrialFunction(chain, x);

   // get a proposal point from the trial function and scale it
   for (int i = 0; i < fMCMCNParameters; ++i)
      x[i] = fMCMCx[chain * fMCMCNParameters + i] + x[i] * (fMCMCBoundaryMax.at(i) - fMCMCBoundaryMin.at(i));

   // check if the point is in the correct volume.
   for (int i = 0; i < fMCMCNParameters; ++i)
      if ((x[i] < fMCMCBoundaryMin[i]) || (x[i] > fMCMCBoundaryMax[i]))
         return false;

   return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(int ichain, int ipar, std::vector<double> &x)
{
   // get unscaled random point in the dimension of the chosen
   // parameter. this point might not be in the correct volume.
   double proposal = MCMCTrialFunctionSingle(ichain, ipar);

   // copy the old point into the new
   for (int i = 0; i < fMCMCNParameters; ++i)
      x[i] = fMCMCx[ichain * fMCMCNParameters + i];

   // modify the parameter under study
   x[ipar] += proposal * (fMCMCBoundaryMax[ipar] - fMCMCBoundaryMin[ipar]);

   // check if the point is in the correct volume.
   if ((x[ipar] < fMCMCBoundaryMin[ipar]) || (x[ipar] > fMCMCBoundaryMax[ipar]))
      return false;

   return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain, int parameter)
{
   // calculate index
   int index = chain * fMCMCNParameters;

   fMCMCCurrentChain = chain;

   // increase counter
   fMCMCNIterations[chain]++;

   // get proposal point
   if (!MCMCGetProposalPointMetropolis(chain, parameter, fMCMCxLocal))
   {
      // increase counter
      fMCMCNTrialsFalse[chain * fMCMCNParameters + parameter]++;

      // execute user code for every point
      MCMCCurrentPointInterface(fMCMCxLocal, chain, false);

      return false;
   }

   // calculate probabilities of the old and new points
   double p0 = fMCMCprob[chain];
   double p1 = LogEval(fMCMCxLocal);

   // flag for accept
   bool accept = false;

   // if the new point is more probable, keep it ...
   if (p1 >= p0)
      accept = true;
   // ... or else throw dice.
   else
   {
      double r = log(fRandom->Rndm());

      if(r < p1 - p0)
         accept = true;
   }

   // fill the new point
   if(accept)
   {
      // increase counter
      fMCMCNTrialsTrue[chain * fMCMCNParameters + parameter]++;

      // copy the point
      for(int i = 0; i < fMCMCNParameters; ++i)
      {
         // save the point
         fMCMCx[index + i] = fMCMCxLocal[i];

         // save the probability of the point
         fMCMCprob[chain] = p1;
      }
   }
   else
   {
      // increase counter
      fMCMCNTrialsFalse[chain * fMCMCNParameters + parameter]++;
   }

   // execute user code for every point
   MCMCCurrentPointInterface(fMCMCxLocal, chain, accept);

   return accept;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetNewPointMetropolis(int chain)
{
   // calculate index
   int index = chain * fMCMCNParameters;

   fMCMCCurrentChain = chain;

   // increase counter
   fMCMCNIterations[chain]++;

   // get proposal point
   if (!MCMCGetProposalPointMetropolis(chain, fMCMCxLocal))
   {
      // increase counter
      for (int i = 0; i < fMCMCNParameters; ++i)
         fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++;

      // execute user code for every point
      MCMCCurrentPointInterface(fMCMCxLocal, chain, false);

      return false;
   }

   // calculate probabilities of the old and new points
   double p0 = fMCMCprob[chain];
   double p1 = LogEval(fMCMCxLocal);

   // flag for accept
   bool accept = false;

   // if the new point is more probable, keep it ...
   if (p1 >= p0)
      accept = true;

   // ... or else throw dice.
   else
   {
      double r = log(fRandom->Rndm());

      if(r < p1 - p0)
         accept = true;
   }

   // fill the new point
   if(accept)
   {
      // increase counter
      for (int i = 0; i < fMCMCNParameters; ++i)
         fMCMCNTrialsTrue[chain * fMCMCNParameters + i]++;

      // copy the point
      for(int i = 0; i < fMCMCNParameters; ++i)
      {
         // save the point
         fMCMCx[index + i] = fMCMCxLocal[i];

         // save the probability of the point
         fMCMCprob[chain] = p1;
      }
   }
   else
   {
      // increase counter
      for (int i = 0; i < fMCMCNParameters; ++i)
         fMCMCNTrialsFalse[chain * fMCMCNParameters + i]++;
   }

   // execute user code for every point
   MCMCCurrentPointInterface(fMCMCxLocal, chain, accept);

   return accept;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainCheckMaximum()
{
   // loop over all chains
   for (int i = 0; i < fMCMCNChains; ++i)
   {
      // check if new maximum is found or chain is at the beginning
      if (fMCMCprob[i] > fMCMCprobMax[i] || fMCMCNIterations[i] == 1)
      {
         // copy maximum value
         fMCMCprobMax[i] = fMCMCprob[i];

         // copy mode of chain
         for (int j = 0; j < fMCMCNParameters; ++j)
            fMCMCxMax[i * fMCMCNParameters + j] = fMCMCx[i * fMCMCNParameters + j];
      }
   }
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainUpdateStatistics()
{
   // calculate number of entries in this part of the chain
   int npoints = fMCMCNTrialsTrue[0] + fMCMCNTrialsFalse[0];

   // length of vectors
   int nentries = fMCMCNParameters * fMCMCNChains;

   // loop over all parameters of all chains
   for (int i = 0; i < nentries; ++i) {
      // calculate mean value of each parameter in the chain for this part
      fMCMCxMean[i] += (fMCMCx[i] - fMCMCxMean[i]) / double(npoints);

      // calculate variance of each chain for this part
      if (npoints > 1)
         fMCMCxVar[i] = (1.0 - 1./double(npoints)) * fMCMCxVar[i]
            + (fMCMCx[i] - fMCMCxMean[i]) * (fMCMCx[i] - fMCMCxMean[i]) / double(npoints - 1);
   }

   // loop over chains
   for (int i = 0; i < fMCMCNChains; ++i) {
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
   // check if histograms are supposed to be filled
   if (!fMCMCFlagFillHistograms)
      return;

   // loop over chains
   for (int i = 0; i < fMCMCNChains; ++i)
   {
      // fill each 1-dimensional histogram (if supposed to be filled)
      for (int j = 0; j < fMCMCNParameters; ++j)
         if (fMCMCFlagsFillHistograms.at(j))
            fMCMCH1Marginalized[j]->Fill(fMCMCx[i * fMCMCNParameters + j]);

      // fill each 2-dimensional histogram (if supposed to be filled)
      int counter = 0;

      for (int j = 0; j < fMCMCNParameters; ++j)
         for (int k = 0; k < j; ++k)
         {
            if (fMCMCFlagsFillHistograms.at(j) && fMCMCFlagsFillHistograms.at(k))
               fMCMCH2Marginalized[counter]->Fill(fMCMCx[i*fMCMCNParameters+k],fMCMCx[i* fMCMCNParameters+j]);
            counter ++;
         }
   }
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainTestConvergenceAllChains()
{
   // calculate number of entries in this part of the chain
   int npoints = fMCMCNTrialsTrue[0] + fMCMCNTrialsFalse[0];

   if (fMCMCNChains > 1 && npoints > 1)
   {
      // define flag for convergence
      bool flag_convergence = true;

      // loop over parameters
      for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters){
         // extract means and variances
         std::vector<double> means(fMCMCNChains);
         std::vector<double> variances(fMCMCNChains);

         // loop over chains
         for (int ichains = 0; ichains < fMCMCNChains; ++ichains) {
            // get parameter index
            int index = ichains * fMCMCNParameters + iparameters;
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
         fMCMCNIterationsConvergenceGlobal = fMCMCNIterations[0] / fMCMCNParameters;
   }
}
// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainWriteChains()
{
   // loop over all chains
   for (int i = 0; i < fMCMCNChains; ++i)
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
   int ndigits = (int)log10(fMCMCNParameters) +1;
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
   int counterupdate = 1;        // after how many iterations is an update needed?
   bool convergence = false;     // convergence reached?
   bool flagefficiency = false;  // efficiency reached?

   // array of efficiencies
   std::vector<double> efficiency;
   efficiency.assign(fMCMCNParameters * fMCMCNChains, 0.0);

   // how often to check convergence and efficiencies?
   // it's either every fMCMCNParameters*nMCMCNIterationsUpdate (for 5 parameters the default would be 5000)
   // or it's fMCMCNIterationsUpdateMax (10000 by default)
   // whichever of the two is smaller
   int updateLimit = ( fMCMCNIterationsUpdateMax<fMCMCNIterationsUpdate*(fMCMCNParameters)  && fMCMCNIterationsUpdateMax>0 ) ?
         fMCMCNIterationsUpdateMax : fMCMCNIterationsUpdate*(fMCMCNParameters);

   // loop over chains
   for (int ichains = 0; ichains < fMCMCNChains; ++ichains) {
      // loop over parameters
      for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter){
         // global index of the parameter (throughout all the chains)
         int index = ichains * fMCMCNParameters + iparameter;
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
   fMCMCCycle = 0;

   // run chain ...
   // (a) for at least a minimum number of iterations,
   // (b) until a maximum number of iterations is reached,
   // (c) or until convergence is reached and the efficiency is in the
   //     specified region
   while (fMCMCCurrentIteration < fMCMCNIterationsPreRunMin || (fMCMCCurrentIteration < fMCMCNIterationsMax && !(convergence && flagefficiency)))
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
         for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
         {
            // loop over chains
            for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
               MCMCGetNewPointMetropolis(ichains, iparameters);

            // search for global maximum
            MCMCInChainCheckMaximum();
         }
      }

      // if the flag is not set then run over the parameters at the same time.
      else
      {
         // loop over chains
         for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
            // get new point
            MCMCGetNewPointMetropolis(ichains);

         // search for global maximum
         MCMCInChainCheckMaximum();
      }

      //-------------------------------------------
      // print out message to log
      //-------------------------------------------

      // progress printout
      if ( fMCMCCurrentIteration > 0 && fMCMCCurrentIteration % fMCMCNIterationsUpdate == 0 )
         BCLog::OutDetail(Form(" --> Iteration %i", fMCMCNIterations[0]/fMCMCNParameters));

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
      if ( counterupdate % updateLimit == 0 && counterupdate > 0 && fMCMCCurrentIteration >= fMCMCNIterationsPreRunMin)
//      if ( fMCMCCurrentIteration % fMCMCNIterationsUpdate == 0 && counterupdate > 1 && fMCMCCurrentIteration >= fMCMCNIterationsPreRunMin)
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
            for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
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
         for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
         {
            // loop over parameters
            for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter)
            {
               // global index of the parameter (throughout all the chains)
               int index = ichains * fMCMCNParameters + iparameter;

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

                  BCLog::OutDetail(Form("         Efficiency of parameter %i dropped below %.2lf%% (eps = %.2lf%%) in chain %i. Set scale to %.4g",
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

                  BCLog::OutDetail(Form("         Efficiency of parameter %i above %.2lf%% (eps = %.2lf%%) in chain %i. Set scale to %.4g",
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
         for (int ichains = 0; ichains < fMCMCNChains; ++ichains) {
            // loop over parameters
            for (int iparameter = 0; iparameter < fMCMCNParameters; ++iparameter){
               // global index of the parameter (throughout all the chains)
               int index = ichains * fMCMCNParameters + iparameter;
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

      if (counterupdate == 0 && fMCMCCurrentIteration > 1)
        fMCMCCycle++;
      fMCMCCurrentIteration++;
      counterupdate++;

   } // end of running

   // decrease counter by one since it didn't really run that long
   //   fMCMCCurrentIteration--;
   //   counterupdate--;

   // did we check convergence at least once ?
   if (fMCMCCurrentIteration<updateLimit)
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

   for (int i = 0; i < fMCMCNParameters; ++i)
      efficiencies.push_back(0.);

   BCLog::OutDetail(" --> Average efficiencies:");
   for (int i = 0; i < fMCMCNParameters; ++i)
   {
      for (int j = 0; j < fMCMCNChains; ++j)
         efficiencies[i] += efficiency[j * fMCMCNParameters + i] / double(fMCMCNChains);

      BCLog::OutDetail(Form(TString::Format(" -->      parameter %%%di :  %%.02f%%%%",ndigits), i, 100. * efficiencies[i]));
   }


   // print scale factors
   std::vector<double> scalefactors;

   for (int i = 0; i < fMCMCNParameters; ++i)
      scalefactors.push_back(0.0);

   BCLog::OutDetail(" --> Average scale factors:");
   for (int i = 0; i < fMCMCNParameters; ++i)
   {
      for (int j = 0; j < fMCMCNChains; ++j)
         scalefactors[i] += fMCMCTrialFunctionScaleFactor[j * fMCMCNParameters + i] / double(fMCMCNChains);

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
   fMCMCCycle = 0;

   // perform run
   BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

//   int counterupdate = 0;
//   bool convergence = false;
//   bool flagefficiency = false;

//   std::vector<double> efficiency;

//   for (int i = 0; i < fMCMCNParameters; ++i)
//      for (int j = 0; j < fMCMCNChains; ++j)
//         efficiency.push_back(0.0);

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
   for (fMCMCCurrentIteration = 1; fMCMCCurrentIteration <= fMCMCNIterationsRun; ++fMCMCCurrentIteration)
   {
      if ( (fMCMCCurrentIteration)%nwrite == 0 )
         BCLog::OutDetail(Form(" --> iteration number %i (%.2f%%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration)/(double)fMCMCNIterationsRun*100.));

      // if the flag is set then run over the parameters one after the other.
      if (fMCMCFlagOrderParameters)
      {
         // loop over parameters
         for (int iparameters = 0; iparameters < fMCMCNParameters; ++iparameters)
         {
            // loop over chains
            for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
               MCMCGetNewPointMetropolis(ichains, iparameters);

            // reset current chain
            fMCMCCurrentChain = -1;

            // update search for maximum
            MCMCInChainCheckMaximum();

         } // end loop over all parameters

         // check if the current iteration is consistent with the lag
         if ( fMCMCCurrentIteration % fMCMCNLag == 0)
         {
            // fill histograms
            MCMCInChainFillHistograms();

            // write chain to file
            if (fMCMCFlagWriteChainToFile)
               MCMCInChainWriteChains();

            // do anything interface
            MCMCIterationInterface();
         }
      }
      // if the flag is not set then run over the parameters at the same time.
      else
      {
         // loop over chains
         for (int ichains = 0; ichains < fMCMCNChains; ++ichains)
            // get new point
            MCMCGetNewPointMetropolis(ichains);

         // reset current chain
         fMCMCCurrentChain = -1;

         // update search for maximum
         MCMCInChainCheckMaximum();

         // check if the current iteration is consistent with the lag
         if (fMCMCCurrentIteration % fMCMCNLag == 0)
         {
            // fill histograms
            MCMCInChainFillHistograms();

            // write chain to file
            if (fMCMCFlagWriteChainToFile)
               MCMCInChainWriteChains();

            // do anything interface
            MCMCIterationInterface();
         }
      }

   } // end run

   // print convergence status
   BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

   // print modes

   // find global maximum
   double probmax = fMCMCprobMax.at(0);
   int probmaxindex = 0;

   // loop over all chains and find the maximum point
   for (int i = 1; i < fMCMCNChains; ++i) {
      if (fMCMCprobMax.at(i) > probmax)
      {
         probmax = fMCMCprobMax.at(i);
         probmaxindex = i;
      }
   }

   BCLog::OutDetail(" --> Global mode from MCMC:");
   BCLog::OutDebug(Form(" --> Posterior value: %g", probmax));
   int ndigits = (int) log10(fMCMCNParameters);
   for (int i = 0; i < fMCMCNParameters; ++i)
      BCLog::OutDetail(Form( TString::Format(" -->      parameter %%%di:   %%.4g", ndigits+1),
            i, fMCMCxMax[probmaxindex * fMCMCNParameters + i]));

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
   for (int j = 0; j < fMCMCNChains; ++j)
   {
      fMCMCNIterations[j]  = 0;
      fMCMCNTrialsTrue[j]  = 0;
      fMCMCNTrialsFalse[j] = 0;
      fMCMCprobMean[j]         = 0;
      fMCMCprobVar[j]     = 0;

      for (int k = 0; k < fMCMCNParameters; ++k)
      {
         fMCMCNTrialsTrue[j * fMCMCNParameters + k]  = 0;
         fMCMCNTrialsFalse[j * fMCMCNParameters + k] = 0;
      }
   }

   // reset marginalized distributions
   for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
      if (fMCMCH1Marginalized[i])
         fMCMCH1Marginalized[i]->Reset();

   for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
      if (fMCMCH2Marginalized[i])
         fMCMCH2Marginalized[i]->Reset();

   fMCMCRValue = 100;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCAddParameter(double min, double max)
{
   // add the boundaries to the corresponding vectors
   fMCMCBoundaryMin.push_back(min);
   fMCMCBoundaryMax.push_back(max);

   // set flag for individual parameters
   fMCMCFlagsFillHistograms.push_back(true);

   // increase the number of parameters by one
   fMCMCNParameters++;

   // return the number of parameters
   return fMCMCNParameters;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInitializeMarkovChains()
{
   // evaluate function at the starting point
   std::vector<double> x0;

   for (int j = 0; j < fMCMCNChains; ++j)
   {
      x0.clear();
      for (int i = 0; i < fMCMCNParameters; ++i)
         x0.push_back(fMCMCx[j * fMCMCNParameters + i]);
      fMCMCprob[j] = LogEval(x0);
   }

   x0.clear();
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCResetResults()
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

   for (int i = 0; i < int(fMCMCH1Marginalized.size()); ++i)
      if (fMCMCH1Marginalized[i])
         delete fMCMCH1Marginalized[i];

   for (int i = 0; i < int(fMCMCH2Marginalized.size()); ++i)
      if (fMCMCH2Marginalized[i])
         delete fMCMCH2Marginalized[i];

   // clear plots
   fMCMCH1Marginalized.clear();
   fMCMCH2Marginalized.clear();

   // reset flags
   fMCMCFlagPreRun = true;
   fMCMCFlagRun = false;
   fMCMCFlagConvergenceGlobal = false;

   // no errors
   return 1;
}

// --------------------------------------------------------
int BCEngineMCMC::MCMCInitialize()
{
   // reset values
   MCMCResetResults();

   // free memory for vectors
   fMCMCNIterations.assign(fMCMCNChains, 0);
   fMCMCprobMean.assign(fMCMCNChains, 0);
   fMCMCprobVar.assign(fMCMCNChains, 0);
   fMCMCprob.assign(fMCMCNChains, -1.0);
   fMCMCprobMax.assign(fMCMCNChains, -1.0);

   fMCMCNTrialsTrue.assign(fMCMCNChains * fMCMCNParameters, 0);
   fMCMCNTrialsFalse.assign(fMCMCNChains * fMCMCNParameters, 0);
   fMCMCxMax.assign(fMCMCNChains * fMCMCNParameters, 0.);
   fMCMCxMean.assign(fMCMCNChains * fMCMCNParameters, 0);
   fMCMCxVar.assign(fMCMCNChains * fMCMCNParameters, 0);

   fMCMCRValueParameters.assign(fMCMCNParameters, 0.);

   if (fMCMCTrialFunctionScaleFactorStart.size() == 0)
      fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains * fMCMCNParameters, 1.0);
   else
      for (int i = 0; i < fMCMCNChains; ++i)
         for (int j = 0; j < fMCMCNParameters; ++j)
            fMCMCTrialFunctionScaleFactor.push_back(fMCMCTrialFunctionScaleFactorStart.at(j));

   // set initial position
   if (fMCMCFlagInitialPosition == 2) // user defined points
   {
      // define flag
      bool flag = true;

      // check the length of the array of initial positions
      if (int(fMCMCInitialPosition.size()) != (fMCMCNChains * fMCMCNParameters))
      {
         BCLog::OutError("BCEngine::MCMCInitialize : Length of vector containing initial positions does not have required length.");
         flag = false;
      }

      // check the boundaries
      if (flag)
      {
         for (int j = 0; j < fMCMCNChains; ++j)
            for (int i = 0; i < fMCMCNParameters; ++i)
               if (fMCMCInitialPosition[j * fMCMCNParameters + i] < fMCMCBoundaryMin[i] ||
                     fMCMCInitialPosition[j * fMCMCNParameters + i] > fMCMCBoundaryMax[i])
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
      for (int j = 0; j < fMCMCNChains; ++j)
         for (int i = 0; i < fMCMCNParameters; ++i)
            fMCMCx.push_back(fMCMCBoundaryMin[i] + .5 * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i]));

   else if (fMCMCFlagInitialPosition == 2) // user defined
   {
      for (int j = 0; j < fMCMCNChains; ++j)
         for (int i = 0; i < fMCMCNParameters; ++i)
            fMCMCx.push_back(fMCMCInitialPosition.at(j * fMCMCNParameters + i));
   }

   else
      for (int j = 0; j < fMCMCNChains; ++j) // random number (default)
         for (int i = 0; i < fMCMCNParameters; ++i)
            fMCMCx.push_back(fMCMCBoundaryMin[i] + fRandom->Rndm() * (fMCMCBoundaryMax[i] - fMCMCBoundaryMin[i]));

   // copy the point of the first chain
   fMCMCxLocal.clear();
   for (int i = 0; i < fMCMCNParameters; ++i)
      fMCMCxLocal.push_back(fMCMCx[i]);

   // define 1-dimensional histograms for marginalization
   for(int i = 0; i < fMCMCNParameters; ++i)
   {
      TH1D * h1 = 0;
      if (fMCMCFlagsFillHistograms.at(i))
         h1 = new TH1D(TString::Format("h1_%d_parameter_%i", BCLog::GetHIndex() ,i), "",
               fMCMCH1NBins[i], fMCMCBoundaryMin[i], fMCMCBoundaryMax[i]);
      fMCMCH1Marginalized.push_back(h1);
   }

   for(int i = 0; i < fMCMCNParameters; ++i)
      for (int k = 0; k < i; ++k)
      {
         TH2D * h2 = 0;
         if (fMCMCFlagsFillHistograms.at(i) && fMCMCFlagsFillHistograms.at(k))
            h2 = new TH2D(Form("h2_%d_parameters_%i_vs_%i", BCLog::GetHIndex(), i, k), "",
                  fMCMCH1NBins[k], fMCMCBoundaryMin[k], fMCMCBoundaryMax[k],
                  fMCMCH1NBins[i], fMCMCBoundaryMin[i], fMCMCBoundaryMax[i] );
         fMCMCH2Marginalized.push_back(h2);
      }

   return 1;
}

// ---------------------------------------------------------
int BCEngineMCMC::SetMarginalized(int index, TH1D * h)
{
   if((int)fMCMCH1Marginalized.size()<=index)
      return 0;

   if(h==0)
      return 0;

   if((int)fMCMCH1Marginalized.size()==index)
      fMCMCH1Marginalized.push_back(h);
   else
      fMCMCH1Marginalized[index]=h;

   return index;
}

// ---------------------------------------------------------
int BCEngineMCMC::SetMarginalized(int index1, int index2, TH2D * h)
{
   int counter = 0;
   int index = 0;

   // search for correct combination
   for(int i = 0; i < fMCMCNParameters; i++)
      for (int j = 0; j < i; ++j)
      {
         if(j == index1 && i == index2)
            index = counter;
         counter++;
      }

   if((int)fMCMCH2Marginalized.size()<=index)
      return 0;

   if(h==0)
      return 0;

   if((int)fMCMCH2Marginalized.size()==index)
      fMCMCH2Marginalized.push_back(h);
   else
      fMCMCH2Marginalized[index]=h;

   return index;
}

// ---------------------------------------------------------

