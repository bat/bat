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
#include <TFile.h>
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
#include <TObject.h>
#include <TKey.h>
#include <TList.h>

#include <math.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const char * name)
	: fMCMCFlagWriteChainToFile(false)
	, fMCMCFlagWritePreRunToFile(false)
	, fMCMCOutputFile(0)
	, fMCMCOutputFilename("")
	, fMCMCOutputFileOption("")
	, fMCMCOutputFileAutoclose(false)
	, fMCMCScaleFactorLowerLimit(0)
	, fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max())
	, fMCMCEfficiencyMin(0.15)
	, fMCMCEfficiencyMax(0.50)
	, fMCMCFlagInitialPosition(1)
	, fMCMCFlagOrderParameters(true)
	, fMCMCPhase(BCEngineMCMC::MCMCUnsetPhase)
	, fMCMCRValueUseStrict(false)
	, fMCMCRValueCriterion(0.1)
	, fMCMCRValueParametersCriterion(0.1)
	, fRandom(new TRandom3())
	, fMCMCTree(0)
	, fMCMCTreeLoaded(false)
	, fMCMCTreeReuseObservables(true)
	, fParameterTree(0)
{
	SetName(name);
	MCMCSetPrecision(BCEngineMCMC::kMedium);
	MCMCSetRandomSeed(0);
	
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(std::string filename, std::string name, bool reuseObservables)
	: fMCMCFlagWriteChainToFile(false)
	, fMCMCFlagWritePreRunToFile(false)
	, fMCMCOutputFile(0)
	, fMCMCOutputFilename("")
	, fMCMCOutputFileOption("")
	, fMCMCOutputFileAutoclose(false)
	, fMCMCScaleFactorLowerLimit(0)
	, fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max())
	, fMCMCEfficiencyMin(0.15)
	, fMCMCEfficiencyMax(0.50)
	, fMCMCFlagInitialPosition(1)
	, fMCMCFlagOrderParameters(true)
	, fMCMCPhase(BCEngineMCMC::MCMCUnsetPhase)
	, fMCMCRValueUseStrict(false)
	, fMCMCRValueCriterion(0.1)
	, fMCMCRValueParametersCriterion(0.1)
	, fRandom(new TRandom3())
	, fMCMCTree(0)
	, fMCMCTreeLoaded(false)
	, fMCMCTreeReuseObservables(true)
	, fParameterTree(0)
{
	SetName(name);
	MCMCSetPrecision(BCEngineMCMC::kMedium);
	MCMCSetRandomSeed(0);
	LoadMCMC(filename,"","",reuseObservables);
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
		fMCMCNChains                          = 1;
		fMCMCNLag                             = 1;
		fMCMCNIterationsPreRunMin             = 1500;
		fMCMCNIterationsPreRunMax             = 10000;
		fMCMCNIterationsRun                   = 10000;
		fMCMCNIterationsEfficiencyCheck       = 500;
		fMCMCNIterationsConvergenceCheck      = 1000;
		fMCMCNIterationsClearConvergenceStats = 5000;
		break;
		
	case BCEngineMCMC::kQuick:
		fMCMCNChains                          = 2;
		fMCMCNLag                             = 1;
		fMCMCNIterationsPreRunMin             = 1500;
		fMCMCNIterationsPreRunMax             = 10000;
		fMCMCNIterationsRun                   = 10000;
		fMCMCNIterationsEfficiencyCheck       = 500;
		fMCMCNIterationsConvergenceCheck      = 1000;
		fMCMCNIterationsClearConvergenceStats = 5000;
		break;

	case  BCEngineMCMC::kMedium:
		fMCMCNChains                          = 5;
		fMCMCNLag                             = 1;
		fMCMCNIterationsPreRunMin             = 1500;
		fMCMCNIterationsPreRunMax             = 100000;
		fMCMCNIterationsRun                   = 100000;
		fMCMCNIterationsEfficiencyCheck       = 500;
		fMCMCNIterationsConvergenceCheck      = 1000;
		fMCMCNIterationsClearConvergenceStats = 5000;
		break;

	case  BCEngineMCMC::kHigh:
		fMCMCNChains                          = 10;
		fMCMCNLag                             = 10;
		fMCMCNIterationsPreRunMin             = 5000;
		fMCMCNIterationsPreRunMax             = 1000000;
		fMCMCNIterationsRun                   = 1000000;
		fMCMCNIterationsEfficiencyCheck       = 1000;
		fMCMCNIterationsConvergenceCheck      = 1000;
		fMCMCNIterationsClearConvergenceStats = 5000;
		break;

	case  BCEngineMCMC::kVeryHigh:
		fMCMCNChains                          = 10;
		fMCMCNLag                             = 10;
		fMCMCNIterationsPreRunMin             = 10000;
		fMCMCNIterationsPreRunMax             = 10000000;
		fMCMCNIterationsRun                   = 10000000;
		fMCMCNIterationsEfficiencyCheck       = 1000;
		fMCMCNIterationsConvergenceCheck      = 1000;
		fMCMCNIterationsClearConvergenceStats = 5000;
		break;
	}

	// re-initialize
	MCMCInitialize();
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetPrecision(const BCEngineMCMC * other) {
	fMCMCNChains                          = other -> MCMCGetNChains();
	fMCMCNLag                             = other -> MCMCGetNLag();
	fMCMCNIterationsPreRunMin             = other -> MCMCGetNIterationsPreRunMin();
	fMCMCNIterationsPreRunMax             = other -> MCMCGetNIterationsPreRunMax();
	fMCMCNIterationsRun                   = other -> MCMCGetNIterationsRun();
	fMCMCNIterationsEfficiencyCheck       = other -> MCMCGetNIterationsEfficiencyCheck();
	fMCMCNIterationsConvergenceCheck      = other -> MCMCGetNIterationsConvergenceCheck();
	fMCMCNIterationsClearConvergenceStats = other -> MCMCGetNIterationsClearConvergenceStats();
	fMCMCRValueCriterion                  = other -> MCMCGetRValueCriterion();
	fMCMCRValueParametersCriterion        = other -> MCMCGetRValueParametersCriterion();
	fMCMCEfficiencyMin                    = other -> MCMCGetMinimumEfficiency();
	fMCMCEfficiencyMax                    = other -> MCMCGetMaximumEfficiency();

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
   fMCMCNIterationsConvergenceCheck          = other.fMCMCNIterationsConvergenceCheck;
   fMCMCNIterationsClearConvergenceStats     = other.fMCMCNIterationsClearConvergenceStats;
   fMCMCNIterationsConvergenceGlobal         = other.fMCMCNIterationsConvergenceGlobal;
   fMCMCNIterationsPreRunMax                 = other.fMCMCNIterationsPreRunMax;
   fMCMCNIterationsRun                       = other.fMCMCNIterationsRun;
   fMCMCNIterationsPreRunMin                 = other.fMCMCNIterationsPreRunMin;
   fMCMCNTrialsTrue                          = other.fMCMCNTrialsTrue;
   fMCMCNTrials                              = other.fMCMCNTrials;
   fMCMCFlagWriteChainToFile                 = other.fMCMCFlagWriteChainToFile;
   fMCMCFlagWritePreRunToFile                = other.fMCMCFlagWritePreRunToFile;
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
   fMCMCObservables                          = other.fMCMCObservables;
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

	 fMCMCTree       = 0;
	 fMCMCTreeLoaded = false;
	 fMCMCTreeReuseObservables = true;
	 fParameterTree  = 0;
	 fMCMCOutputFile = 0;
	 fMCMCOutputFilename      = other.fMCMCOutputFilename;
	 fMCMCOutputFileOption    = other.fMCMCOutputFileOption;
	 fMCMCOutputFileAutoclose = other.fMCMCOutputFileAutoclose;

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
void BCEngineMCMC::SetNbins(unsigned int nbins)
{
	fParameters.SetNBins(nbins);
	fObservables.SetNBins(nbins);
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChain(bool flag) {
	if (flag)
		BCLog::OutError("BCEngineMCMC::WriteMarkovChain: To turn on output use WriteMarkovChain(filename,option).");
	fMCMCFlagWriteChainToFile = false;
	fMCMCFlagWritePreRunToFile = false;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChain(std::string filename, std::string option, bool autoclose) {
	if (filename.empty()) {
		BCLog::OutError("BCEngineMCMC::WriteMarkovChain: You must specify the filename when turning on Markov chain output.");
		return WriteMarkovChain(false);
	}
	fMCMCOutputFilename = filename;
	fMCMCOutputFileOption = option;
	fMCMCOutputFileAutoclose = autoclose;
	fMCMCFlagWriteChainToFile = true;
	fMCMCFlagWritePreRunToFile = true;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarginalizedDistributions(std::string filename, std::string option) {
	// remember current directory
	TDirectory * dir = gDirectory;

	TFile * fOut = TFile::Open(filename.c_str(),option.c_str());
	if (!fOut) {
		BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: Could not open output file.");
		return;
	}
	if (!fOut->IsWritable()) {
		BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: File must be opened in writeable mode.");
		return;
	}
	for (unsigned i=0; i<GetNVariables(); ++i) {
		if (MarginalizedHistogramExists(i))
			fOut -> WriteTObject(GetMarginalizedHistogram(i));
		for (unsigned j=0; j<GetNVariables(); ++j)
			if (MarginalizedHistogramExists(i,j))
				fOut -> WriteTObject(GetMarginalizedHistogram(i,j));
	}
	fOut -> Write();
	fOut -> Close();
	gDirectory = dir;
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
	if ( index<GetNParameters() ) // Marginalization of model parameter
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for parameter %s", GetParameter(index)->GetName().data()));
	else if ( index-GetNParameters()<GetNObservables() ) // Marginalization of user-defined observable
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for user-observable %s", GetObservable(index-GetNParameters())->GetName().data()));
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
		if ( i<GetNParameters() ) // Marginalization of model parameter
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for parameter %s", GetParameter(i)->GetName().data()));
		else if ( i-GetNParameters()<GetNObservables() ) // Marginalization of user-defined observable
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for user-observable %s", GetObservable(i-GetNParameters())->GetName().data()));
	}
	if (!fH1Marginalized[j]) {
		if ( j<GetNParameters() ) // Marginalization of model parameter
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for parameter %s", GetParameter(j)->GetName().data()));
		else if ( j-GetNParameters()<GetNObservables() ) // Marginalization of user-defined observable
			BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distributions not stored for user-observable %s", GetObservable(j-GetNParameters())->GetName().data()));
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
	if (GetBestFitParameters().size() == GetNParameters())
		hprob -> SetGlobalMode(GetBestFitParameter(index));
	
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
	if (GetBestFitParameters().size() == GetNParameters())
		hprob -> SetGlobalMode(GetBestFitParameter(i),GetBestFitParameter(j));
		
	return hprob;
}

// ---------------------------------------------------------
BCVariable * BCEngineMCMC::GetVariable(unsigned index) const {
	if (index<GetNParameters())
		return GetParameter(index);
	if (index-GetNParameters()<GetNObservables())
		return GetObservable(index-GetNParameters());
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
double BCEngineMCMC::GetBestFitParameter(unsigned index)
{

	if (index > GetNParameters() + GetNObservables())
		return -std::numeric_limits<double>::infinity();

	std::vector<double> pars = GetBestFitParameters();

	if (pars.size() != GetNParameters()) {
		BCLog::OutError("BCEngineMCMC::GetBestFitParameter : MCMC not yet run, returning center of the range.");
		return GetParameter(index)->GetRangeCenter();
	}

	// parameter
	if ( index < GetNParameters() )
		return pars[index];
	
	// user-define observable
	CalculateObservables(pars);
	return GetObservable(index-GetNParameters()) -> Value();
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
	 while (it+GetNParameters()<=x0s.end()) {
		 fMCMCInitialPosition.push_back(std::vector<double>(it,it+GetNParameters()));
		 it += GetNParameters();
	 }																		
	 MCMCSetFlagInitialPosition(2);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetFlagFillHistograms(bool flag)
{
	fParameters.FillHistograms(flag);
	fObservables.FillHistograms(flag);	
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
void BCEngineMCMC::InitializeMarkovChainTree(bool replacetree, bool replacefile)
{
	if (fMCMCTree and replacetree) {
		delete fMCMCTree;
		fMCMCTree = 0;
	}
	if (fParameterTree and replacetree) {
		delete fParameterTree;
		fParameterTree = 0;
	}
	if (fMCMCOutputFile and replacefile) {
		fMCMCOutputFile -> Close();
		delete fMCMCOutputFile;
		fMCMCOutputFile = 0;
	}

	TDirectory * dir = gDirectory;
	
	// create file
	if (!fMCMCOutputFile)
		fMCMCOutputFile = TFile::Open(fMCMCOutputFilename.c_str(),fMCMCOutputFileOption.c_str());
	// if failed
	if (!fMCMCOutputFile) {
		BCLog::OutError("BCEngineMCMC::InitializeMarkovChainTree: Could not create new file.");
		WriteMarkovChain(false);
		return;
	}
	// if file mode not writeable
	if (fMCMCOutputFile and !fMCMCOutputFile->IsWritable()) {
		BCLog::OutError("BCEngineMCMC::InitializeMarkovChainTree: File must be opened in a writeable mode.");
		delete fMCMCOutputFile;
		fMCMCOutputFile = 0;
		WriteMarkovChain(false);
		return;
	}

	// change to output file (directory)
	fMCMCOutputFile -> cd();

	if (fMCMCTree)
		// Add existing MCMC tree to output file
		fMCMCOutputFile -> Add(fMCMCTree);
	else {			 
		// or create new MCMC tree (under umbrella of output file)
		fMCMCTree = new TTree(TString::Format("%s_mcmc",GetSafeName().data()),TString::Format("%s_mcmc",GetSafeName().data()));
		fMCMCTree -> Branch("Chain",          &fMCMCTree_Chain,     "chain/i");
		fMCMCTree -> Branch("Iteration",      &fMCMCTree_Iteration, "iteration/i");
		fMCMCTree -> Branch("Phase",          &fMCMCPhase,          "phase/I");
		fMCMCTree -> Branch("LogProbability", &fMCMCTree_Prob,      "log(probability)/D");
		fMCMCTree_Parameters.assign(GetNParameters(),0);
		for (unsigned j=0; j<GetNParameters(); ++j) {
			fMCMCTree -> Branch(GetParameter(j)->GetSafeName().data(),&fMCMCTree_Parameters[j],(GetParameter(j)->GetSafeName()+"/D").data());
			fMCMCTree -> SetAlias(TString::Format("Parameter%i",j),GetParameter(j)->GetSafeName().data());
		}
		fMCMCTree_Observables.assign(GetNObservables(),0);
		for (unsigned j=0; j<GetNObservables(); ++j) {
			fMCMCTree -> Branch(GetObservable(j)->GetSafeName().data(),&fMCMCTree_Observables[j],(GetObservable(j)->GetSafeName()+"/D").data());
			fMCMCTree -> SetAlias(TString::Format("Observable%i",j),GetObservable(j)->GetSafeName().data());
		}
	}

	if (fParameterTree)
		// add existing parameter tree to output file
		fMCMCOutputFile -> Add(fParameterTree);
	else {
		// or create new parameter tree (under umbrella of output file)
		fParameterTree = new TTree(TString::Format("%s_parameters",GetSafeName().data()),TString::Format("%s_parameters",GetSafeName().data()));
		bool p_parameter, p_fill, p_fixed;
		unsigned p_index, p_precision, p_nbins;
		char p_name[200], p_safename[200], p_latexname[200];
		double p_lowerlimit, p_upperlimit, p_fixedvalue;
		fParameterTree -> Branch("parameter",&p_parameter,"parameter/O");
		fParameterTree -> Branch("index",&p_index,"index/i");
		fParameterTree -> Branch("name",p_name,"name/C");
		fParameterTree -> Branch("safe_name",p_safename,"safe_name/C");
		fParameterTree -> Branch("latex_name",p_latexname,"latex_name/C");
		fParameterTree -> Branch("lower_limit",&p_lowerlimit,"lower_limit/D");
		fParameterTree -> Branch("upper_limit",&p_upperlimit,"upper_limit/D");
		fParameterTree -> Branch("precision",&p_precision,"precision/i");
		fParameterTree -> Branch("nbins",&p_nbins,"nbins/i");
		fParameterTree -> Branch("fill",&p_fill,"fill/O");
		fParameterTree -> Branch("fixed",&p_fixed,"fixed/O");
		fParameterTree -> Branch("fixed_value",&p_fixedvalue,"fixed_value/D");
		for (unsigned i=0; i<GetNVariables(); ++i) {
			p_parameter  = (i<GetNParameters());
			p_index      = (p_parameter) ? i : i - GetNParameters();
			strcpy(p_name,      GetVariable(i)->GetName().data());
			strcpy(p_safename,  GetVariable(i)->GetSafeName().data());
			strcpy(p_latexname, GetVariable(i)->GetLatexName().data());
			p_lowerlimit = GetVariable(i) -> GetLowerLimit();
			p_upperlimit = GetVariable(i) -> GetUpperLimit();
			p_precision  = GetVariable(i) -> GetPrecision();
			p_nbins      = GetVariable(i) -> GetNbins();
			p_fill       = GetVariable(i) -> FillHistograms();
			p_fixed      = p_parameter and GetParameter(i)->Fixed();
			p_fixedvalue = (p_parameter) ? GetParameter(i)->GetFixedValue() : 0;
			fParameterTree -> Fill();
		}
		fParameterTree -> AutoSave("SaveSelf");
	}

	// return to old directory
	gDirectory = dir;
}

// --------------------------------------------------------
void BCEngineMCMC::UpdateParameterTree() {
	if (!fParameterTree)
		return;

	unsigned nchains = MCMCGetNChains();
	std::vector<double> scale(MCMCGetNChains(),0);
	std::vector<double> eff(MCMCGetNChains(),0);

	// check for branch existences
	TBranch * b_nchains = fParameterTree -> GetBranch("nchains");
	TBranch * b_scale   = fParameterTree -> GetBranch("scale");

	// if nchains branch doesn't exist, create it
	if (b_nchains == 0)
		b_nchains = fParameterTree -> Branch("nchains",&nchains,"nchains/i");
	// else set 0, so as not to fill
	else
		b_nchains = 0;

	// if scale branch doesn't exist, create it
	if (b_scale == 0)
		b_scale = fParameterTree -> Branch("scale",&(scale.front()),TString::Format("scale[%d]/D",MCMCGetNChains()));
	// else set 0, so as not to fill
	else
		b_scale = 0;

	// create next effiency branch
	unsigned i = 0;
	while (fParameterTree->GetBranch(TString::Format("efficiency_%d",i)))
		++i;
	TBranch * b_eff = fParameterTree -> Branch(TString::Format("efficiency_%d",i),&(eff.front()),TString::Format("efficiency_%d[%d]/D",i,MCMCGetNChains()));

	for (unsigned n=0; n<fParameterTree->GetEntries(); ++n) {
		if (b_nchains)
			b_nchains -> Fill();

		for (unsigned j=0; j<nchains; ++j) {
			scale[j] = (n<GetNParameters()) ? fMCMCTrialFunctionScaleFactor[j][n] : -1;
			eff[j]   = (n<GetNParameters()) ? fMCMCEfficiencies[j][n] : -1;
		}

		if (b_scale)
			b_scale -> Fill();

		b_eff -> Fill();
	}
	fParameterTree -> AutoSave("SaveSelf");
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidMCMCTree(TTree * tree, bool checkObservables) {
	if (!tree)
		return false;
	if (!tree->GetBranch("Chain"))
		return false;
	// if (!tree->GetBranch("Iteration"))
	// 	return false;
	if (!tree->GetBranch("Phase"))
		return false;
	if (!tree->GetBranch("LogProbability"))
		return false;
	for (unsigned i=0; i<GetNParameters(); ++i)
		if (!tree->GetBranch(GetParameter(i)->GetSafeName().data()))
			return false;
	if (checkObservables)
		for (unsigned i=0; i<GetNObservables(); ++i)
			if (!tree->GetBranch(GetObservable(i)->GetSafeName().data()))
				return false;
	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidParameterTree(TTree * tree) {
	if (!tree)
		return false;
	if (!tree->GetBranch("parameter"))
		return false;
	if (!tree->GetBranch("index"))
		return false;
	if (!tree->GetBranch("name"))
		return false;
	// if (!tree->GetBranch("safe_name"))
	// 	return false;
	if (!tree->GetBranch("latex_name"))
	 	return false;
	if (!tree->GetBranch("lower_limit"))
		return false;
	if (!tree->GetBranch("upper_limit"))
		return false;
	if (!tree->GetBranch("precision"))
	 	return false;
	if (!tree->GetBranch("nbins"))
	 	return false;
	if (!tree->GetBranch("fill"))
	 	return false;
	if (!tree->GetBranch("fixed"))
		return false;
	if (!tree->GetBranch("fixed_value"))
		return false;
	// if (!tree->GetBranch("nchains"))
	// 	return false;
	// if (!tree->GetBranch("scale"))
	// 	return false;
	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::LoadParametersFromTree(TTree * partree, bool reuseObservables) {
	bool p_fill, p_fixed;
	unsigned p_precision, p_nbins;
	char p_name[200], p_latexname[200];
	double p_lowerlimit, p_upperlimit, p_fixedvalue;
	partree -> SetBranchAddress("name",p_name);
	partree -> SetBranchAddress("latex_name",p_latexname);
	partree -> SetBranchAddress("lower_limit",&p_lowerlimit);
	partree -> SetBranchAddress("upper_limit",&p_upperlimit);
	partree -> SetBranchAddress("precision",&p_precision);
	partree -> SetBranchAddress("nbins",&p_nbins);
	partree -> SetBranchAddress("fill",&p_fill);
	partree -> SetBranchAddress("fixed",&p_fixed);
	partree -> SetBranchAddress("fixed_value",&p_fixedvalue);
	partree -> BuildIndex("parameter","index");

	// load parameters
	unsigned i = 0;
	while ( partree->GetEntryNumberWithIndex(1,i) >= 0 ) {
		partree -> GetEntryWithIndex(1,i);
		BCParameter * Par = new BCParameter(p_name,p_lowerlimit,p_upperlimit,p_latexname);
		if (p_fixed)
			Par -> Fix(p_fixedvalue);
		Par -> SetPrecision(p_precision);
		Par -> FillHistograms(p_fill);
		Par -> SetNbins(p_nbins);
		AddParameter(Par);
		++i;
	}

	// load user-defined observables
	if (!reuseObservables)
		return true;
	fObservables.Clear(true);
	i = 0;
	while ( partree->GetEntryNumberWithIndex(0,i) >= 0 ) {
		partree -> GetEntryWithIndex(0,i);
		BCObservable * Obs = new BCObservable(p_name,p_lowerlimit,p_upperlimit,p_latexname);
		Obs -> SetPrecision(p_precision);
		Obs -> FillHistograms(p_fill);
		Obs -> SetNbins(p_nbins);
		AddObservable(Obs);
		++i;
	}
	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::ParameterTreeMatchesModel(TTree * partree, bool checkObservables) {
	bool p_fixed;
	char p_name[200];
	double p_lowerlimit, p_upperlimit, p_fixedvalue;
	partree -> SetBranchAddress("name",p_name);
	partree -> SetBranchAddress("lower_limit",&p_lowerlimit);
	partree -> SetBranchAddress("upper_limit",&p_upperlimit);
	partree -> SetBranchAddress("fixed",&p_fixed);
	partree -> SetBranchAddress("fixed_value",&p_fixedvalue);
	partree -> BuildIndex("parameter","index");

	// check parameters
	for (unsigned i=0; i<GetNParameters(); ++i) {
		if ( partree->GetEntryNumberWithIndex(1,i) < 0 ) {
			BCLog::OutError("BCEngineMCMC::ParameterTreeMatchesModel : Parameter tree contains too few entries.");
			return false;
		}
		partree -> GetEntryWithIndex(1,i);
		if ( !GetParameter(i)->IsNamed(p_name) ) {
			BCLog::OutError(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Parameter[%d]'s names do not match.",i));
			return false;
		}
		if ( GetParameter(i)->GetLowerLimit() != p_lowerlimit )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Lower limit of parameter \"%s\" does not match.",p_name));
		if ( GetParameter(i)->GetUpperLimit() != p_upperlimit )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Upper limit of parameter \"%s\" does not match.",p_name));
		if ( GetParameter(i)->Fixed() != p_fixed )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Fixed status of parameter \"%s\" does not match.",p_name));
		if ( GetParameter(i)->Fixed() and GetParameter(i)->GetFixedValue() != p_fixedvalue )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Fixed value of parameter \"%s\" does not match.",p_name));
	}
	if (!checkObservables)
		return true;
	// check observables
	for (unsigned i=0; i<GetNObservables(); ++i) {
		if ( partree->GetEntryNumberWithIndex(0,i) < 0 ) {
			BCLog::OutError("BCEngineMCMC::ParameterTreeMatchesModel : Parameters tree contains too few entries.");
			return false;
		}
		partree -> GetEntryWithIndex(0,i);
		if ( !GetObservable(i)->IsNamed(p_name) ) {
			BCLog::OutError(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Observable[%d]'s names do not match.",i));
			return false;
		}
		if ( GetObservable(i)->GetLowerLimit() != p_lowerlimit )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Lower limit of observable \"%s\" does not match.",p_name));
		if ( GetObservable(i)->GetUpperLimit() != p_upperlimit )
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Upper limit of observable \"%s\" does not match.",p_name));
	}
	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::LoadMCMC(std::string filename, std::string mcmcTreeName, std::string parameterTreeName, bool reuseObservables) {
	// save current directory
	TDirectory * dir = gDirectory;

	TFile * inputfile = TFile::Open(filename.c_str(),"READ");
	if (!inputfile or inputfile->IsZombie()) {
		BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC: Could not open file %s.",filename.data()));
		gDirectory = dir;
		return false;
	}

	// set tree names if empty
	if ( mcmcTreeName.empty() )		// default mcmc tree name
		mcmcTreeName = TString::Format("%s_mcmc",GetSafeName().data());
	if ( parameterTreeName.empty() ) // default parameter tree name
		parameterTreeName = TString::Format("%s_parameters",GetSafeName().data());

	TTree * mcmcTree = NULL;
	inputfile -> GetObject(mcmcTreeName.data(),mcmcTree);
	if ( !mcmcTree )
		BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC : %s does not contain a tree named %s",filename.data(),mcmcTreeName.data()));


	TTree * parTree = NULL;
	inputfile -> GetObject(parameterTreeName.data(),parTree);
	if ( !parTree )
		BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC : %s does not contain a tree named %s",filename.data(),mcmcTreeName.data()));
	
	gDirectory = dir;
	return LoadMCMC(mcmcTree,parTree,reuseObservables);
}

// --------------------------------------------------------
bool BCEngineMCMC::LoadMCMC(TTree * mcmcTree, TTree * parTree, bool reuseObservables) {
	fMCMCTreeLoaded = false;
	fMCMCTreeReuseObservables = reuseObservables;

	if (!mcmcTree or !parTree)
		return false;

	// load parameter tree
	if ( !ValidParameterTree(parTree) ) {
		BCLog::OutError("BCEngineMCMC::LoadMCMC : invalid parameter tree");
		return false;
	}
	if (fParameterTree)
		delete fParameterTree;
	fParameterTree = parTree;

	// if parameters is empty
	if ( GetNParameters() == 0 )
		LoadParametersFromTree(fParameterTree,fMCMCTreeReuseObservables);
	// else check parameter tree
	else if ( !ParameterTreeMatchesModel(fParameterTree,fMCMCTreeReuseObservables) ) {
		BCLog::OutError("BCEngineMCMC::LoadMCMC : Parameter tree does not match model.");
		return false;
	}

	// check mcmc tree
	if ( !ValidMCMCTree(mcmcTree,fMCMCTreeReuseObservables) ) {
		BCLog::OutError("BCEngineMCMC::LoadMCMC : invalid MCMC tree");
		return false;
	}
	if (fMCMCTree)
		delete fMCMCTree;
	fMCMCTree = mcmcTree;

	fMCMCTreeLoaded = true;
	return true;
}

// --------------------------------------------------------
void BCEngineMCMC::Remarginalize(bool autorange) {
	if (!ValidMCMCTree(fMCMCTree))
		return;

	fMCMCTree -> SetBranchAddress("Chain",          &fMCMCTree_Chain);
	fMCMCTree -> SetBranchAddress("Iteration",      &fMCMCTree_Iteration);
	fMCMCTree -> SetBranchAddress("Phase",          &fMCMCPhase);
	fMCMCTree -> SetBranchAddress("LogProbability", &fMCMCTree_Prob);

	fMCMCTree_Parameters.assign(GetNParameters(),0);
	for (unsigned i=0; i<GetNParameters(); ++i)
		fMCMCTree -> SetBranchAddress(GetParameter(i)->GetSafeName().data(),&fMCMCTree_Parameters[i]);
	if (fMCMCTreeReuseObservables) {
		fMCMCTree_Observables.assign(GetNObservables(),0);
		for (unsigned i=0; i<GetNObservables(); ++i)
			fMCMCTree -> SetBranchAddress(GetObservable(i)->GetSafeName().data(),&fMCMCTree_Observables[i]);
	}

	// find out how many chains used to generate tree
	fMCMCNChains = 0;
	for (int n=0; n<fMCMCTree->GetEntries(); ++n) {
		fMCMCTree -> GetEntry(n);
		if (fMCMCNChains>0 and fMCMCTree_Chain==0)
			break;
		if (fMCMCTree_Chain+1>fMCMCNChains)
			fMCMCNChains = fMCMCTree_Chain+1;
	}

	MCMCInitialize();
	MCMCInitializeMarkovChains();
	
	if (autorange) {
		// find min and max
		std::vector<double> XMin(GetNVariables(),+std::numeric_limits<double>::infinity());
		std::vector<double> XMax(GetNVariables(),-std::numeric_limits<double>::infinity());
		for (int n=0; n<fMCMCTree->GetEntries(); ++n) {
			fMCMCTree -> GetEntry(n);

			if (fMCMCPhase<=0)
				continue;

			for (unsigned i=0; i<fMCMCTree_Parameters.size(); ++i) {
				if (fMCMCTree_Parameters[i]<XMin[i])
					XMin[i] = fMCMCTree_Parameters[i];
				if (fMCMCTree_Parameters[i]>XMax[i])
					XMax[i] = fMCMCTree_Parameters[i];
			}
			if (fMCMCTreeReuseObservables) {
				for (unsigned i=0; i<fMCMCTree_Observables.size(); ++i) {
					if (fMCMCTree_Observables[i]<XMin[i])
						XMin[i] = fMCMCTree_Observables[i];
					if (fMCMCTree_Observables[i]>XMax[i])
						XMax[i] = fMCMCTree_Observables[i];
				}
			} else {
				CalculateObservables(fMCMCTree_Parameters);
				for (unsigned i=fMCMCTree_Parameters.size(); i<XMin.size(); ++i) {
					if (((BCObservable*)GetVariable(i))->Value() < XMin[i])
						XMin[i] = ((BCObservable*)GetVariable(i))->Value();
					if (((BCObservable*)GetVariable(i))->Value() > XMax[i])
						XMax[i] = ((BCObservable*)GetVariable(i))->Value();
				}
			}
		}
		// recreate histgrams
		for (unsigned i=0; i<GetNVariables(); ++i) {			
			if (MarginalizedHistogramExists(i)) {
				TString name  = fH1Marginalized[i] -> GetName();
				TString title = fH1Marginalized[i] -> GetTitle();
				int nbins     = fH1Marginalized[i] -> GetNbinsX();
				delete fH1Marginalized[i];
				fH1Marginalized[i] = new TH1D(name,title,nbins,XMin[i],XMax[i]);
				fH1Marginalized[i] -> SetStats(false);
			}
			for (unsigned j=0; j<GetNVariables(); ++j)
				if (MarginalizedHistogramExists(i,j)) {
					TString name  = fH2Marginalized[i][j] -> GetName();
					TString title = fH2Marginalized[i][j] -> GetTitle();
					int nbinsx    = fH2Marginalized[i][j] -> GetNbinsX();
					int nbinsy    = fH2Marginalized[i][j] -> GetNbinsY();
					delete fH2Marginalized[i][j];
					fH2Marginalized[i][j] = new TH2D(name,title,nbinsx,XMin[i],XMax[i],nbinsy,XMin[j],XMax[j]);
					fH2Marginalized[i][j] -> SetStats(false);
				}
		}
	}
	
	for (unsigned n=0; n<fMCMCTree->GetEntries(); ++n) {
		fMCMCTree -> GetEntry(n);

		if (fMCMCTree_Prob > fMCMCLogMaximum) {
			fMCMCBestFitParameters = fMCMCTree_Parameters;
			fMCMCLogMaximum = fMCMCTree_Prob;
		}

		if (fMCMCPhase <= 0)
			continue;

		fMCMCx[fMCMCTree_Chain]    = fMCMCTree_Parameters;
		fMCMCprob[fMCMCTree_Chain] = fMCMCTree_Prob;
		MCMCCurrentPointInterface(fMCMCx[fMCMCTree_Chain], fMCMCTree_Chain, true);
		if (fMCMCTreeReuseObservables)
			fMCMCObservables[fMCMCTree_Chain] = fMCMCTree_Observables;

		if (fMCMCTree_Chain==fMCMCNChains-1) {
			MCMCIterationInterface();
			if ( !fMCMCTreeReuseObservables)
				EvaluateObservables();
			if ( !fH1Marginalized.empty() or !fH2Marginalized.empty() )
				MCMCInChainFillHistograms();
		}

	}
	
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunction(unsigned ichain, std::vector<double> &x)
{
	// call MCMCTrialFunctionSingle() for all parameters by default
	for (unsigned i = 0; i < GetNParameters(); ++i)
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
	for (unsigned i = 0; i < GetNParameters(); ++i)
		if (GetParameter(i)->Fixed())
			x[i] = GetParameter(i) -> GetFixedValue();
		else
			x[i] = fMCMCx[chain][i] + x[i]*(GetParameter(i)->GetUpperLimit()-GetParameter(i)->GetLowerLimit());

	// check if the point is in the correct volume.
	for (unsigned i = 0; i < GetNParameters(); ++i)
		if (!GetParameter(i)->IsValid(x[i]))
			return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned ichain, unsigned ipar, std::vector<double> &x)
{
  // copy the old point into the new
	x = fMCMCx[ichain];
	
  // check if parameter is fixed
  if (GetParameter(ipar)->Fixed()) {
    x[ipar] = GetParameter(ipar)->GetFixedValue();
    return true; // assume that value is inside allowed region
  }

  // get unscaled random point in the dimension of the chosen
  // parameter. this point might not be in the correct volume.
  double proposal = MCMCTrialFunctionSingle(ichain, ipar);

  // modify the parameter under study
  x[ipar] += proposal * (GetParameter(ipar)->GetUpperLimit() - GetParameter(ipar)->GetLowerLimit());

  // check if the point is in the correct volume.
	return GetParameter(ipar) -> IsValid(x[ipar]);
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
			BCLog::OutDebug(Form("Log(likelihood) evaluated to nan or inf in chain %i while varying parameter %s to %.3e",chain,GetParameter(parameter)->GetName().data(),fMCMCThreadLocalStorage[chain].xLocal[parameter]));
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
				for (unsigned i = 0; i < GetNParameters(); ++i)
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

		for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar)

			if ( !GetParameter(ipar)->Fixed() ) {

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
		for (unsigned j = 0; j < GetNParameters(); ++j)
			if (TH1 * h = fH1Marginalized[j])
				h -> Fill(fMCMCx[i][j]);
		// User-defined Observables
		for (unsigned j = 0; j < GetNObservables(); ++j)
			if (TH1 * h = fH1Marginalized[j+GetNParameters()])
				h -> Fill(fMCMCObservables[i][j]);


		////////////////////////////////////////
		// fill each 2-dimensional histogram (if supposed to be filled)

		for (unsigned j = 0; j < GetNParameters(); ++j) {
			// Parameter vs Parameter
			for (unsigned k = j+1; k < GetNParameters(); ++k)
				if (TH2D * h = fH2Marginalized[j][k])
					h -> Fill(fMCMCx[i][j],fMCMCx[i][k]);
			// User-defined Observable vs Parameter
			for (unsigned k = 0; k < GetNObservables(); ++k)
				if (TH2D * h = fH2Marginalized[j][k+GetNParameters()])
					h -> Fill(fMCMCx[i][j],fMCMCObservables[i][k]);
		}
		// User-defined Observable vs User-defined Observable
		for (unsigned j = 0; j < GetNObservables(); ++j)
			for (unsigned k = j+1; k < GetNObservables(); ++k)
				if (TH2D * h = fH2Marginalized[j+GetNParameters()][k+GetNParameters()])
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
	for (unsigned iparameter = 0; iparameter < GetNParameters(); ++iparameter) {
		if ( GetParameter(iparameter) -> Fixed() )
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
	if (!fMCMCTree)
		return;
	// loop over all chains
	for (fMCMCTree_Chain = 0; fMCMCTree_Chain < fMCMCNChains; ++fMCMCTree_Chain) {
		fMCMCTree_Iteration = fMCMCNIterations[fMCMCTree_Chain];
		fMCMCTree_Prob = fMCMCprob[fMCMCTree_Chain];
		fMCMCTree_Parameters = fMCMCx[fMCMCTree_Chain];
		fMCMCTree_Observables = fMCMCObservables[fMCMCTree_Chain];
		fMCMCTree->Fill();
	}
}

//---------------------------------------------------------
void BCEngineMCMC::MCMCCloseOutputFile() {
	if (!fMCMCOutputFile or !fMCMCOutputFile->IsOpen())
		return;
	fMCMCOutputFile -> Write(0,TObject::kWriteDelete);
	fMCMCOutputFile -> Close();
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
	if (fMCMCFlagWritePreRunToFile)
		InitializeMarkovChainTree();
	
	// reset run statistics
	MCMCResetRunStatistics();

	// perform run
	BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsPreRunMax));


	//////////////////////////////////////////////////
	// Adjust scales until all parameters are in correct efficiency range in all chains

	fMCMCPhase = BCEngineMCMC::MCMCPreRunEfficiencyCheck;
	fMCMCCurrentIteration = 0;
	fMCMCNTrials = 0;
	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(GetNParameters(),0));
	fMCMCEfficiencies.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
	bool allEfficient = false;
	bool inefficientScalesAdjustable = true;
	while (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax and !allEfficient and inefficientScalesAdjustable) {

		MCMCGetNewPointMetropolis();
		++fMCMCNTrials;

		if (fMCMCFlagWritePreRunToFile) {
			EvaluateObservables();
			MCMCInChainWriteChains();
		}

		if ( fMCMCNTrials != fMCMCNIterationsEfficiencyCheck)
			continue;

		allEfficient = true;
		inefficientScalesAdjustable = false;
		for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain) {
			for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar) {

				if (GetParameter(ipar)->Fixed())
					continue;

				fMCMCEfficiencies[ichain][ipar] = 1. * fMCMCNTrialsTrue[ichain][ipar] / fMCMCNIterationsEfficiencyCheck;

				if ( fMCMCEfficiencies[ichain][ipar] < fMCMCEfficiencyMin ) {
					// if efficiency too low ...
					
					if (allEfficient)			// print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					
					allEfficient = false;

					fMCMCTrialFunctionScaleFactor[ichain][ipar] /= (fMCMCEfficiencies[ichain][ipar] < 0.5*fMCMCEfficiencyMin) ? 4 : 2;
					
					if ( fMCMCTrialFunctionScaleFactor[ichain][ipar] > fMCMCScaleFactorLowerLimit ) {
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale decreased to %6.2f %%", fParameters.MaxNameLength(), GetParameter(ipar)->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					}	else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorLowerLimit;
						BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale now at lower limit (%6.2f %%)",	fParameters.MaxNameLength(), GetParameter(ipar)->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCScaleFactorLowerLimit));
					}
					
				} else if (fMCMCEfficiencies[ichain][ipar] > fMCMCEfficiencyMax ) {
					// if efficiency too high ...

					if (allEfficient)		 // print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					
					allEfficient = false;

					fMCMCTrialFunctionScaleFactor[ichain][ipar] *= 2;

					if ( fMCMCTrialFunctionScaleFactor[ichain][ipar] < fMCMCScaleFactorUpperLimit ) {
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale increased to %6.2f %%", fParameters.MaxNameLength(), GetParameter(ipar)->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, 100*fMCMCTrialFunctionScaleFactor[ichain][ipar]));
						inefficientScalesAdjustable = true;
					} else {
						fMCMCTrialFunctionScaleFactor[ichain][ipar] = fMCMCScaleFactorUpperLimit;
						BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale now at upper limit (%6.2f %%)", fParameters.MaxNameLength(), GetParameter(ipar)->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCEfficiencies[ichain][ipar], ichain, fMCMCScaleFactorUpperLimit));														 
					}
				}
			}
		}
		fMCMCNTrials = 0;
		fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(GetNParameters(),0));
	}
	if (allEfficient)
		BCLog::OutDetail(Form("     * Efficiency status: Efficiencies within predefined range after %i iterations.",fMCMCCurrentIteration));
	else if (!inefficientScalesAdjustable)
		BCLog::OutWarning(Form("     * Efficiency status: Some efficiencies outside predefined range, but scales are at limits after %i iterations.",fMCMCCurrentIteration));
	else
		BCLog::OutDetail(Form("     * Efficiency status: Some efficiencies outside predefined range, but maximum number of iterations (%i) reached.",fMCMCNIterationsPreRunMax));


	// continue measuring efficiency
	unsigned NTrialsForEff = fMCMCNTrials;

	if (fMCMCNChains > 1) {
		//////////////////////////////////////////////////
		// Run until all chains have converged

		unsigned nIterationsCheckConvergence = fMCMCNIterationsConvergenceCheck;
		if ( fMCMCNIterationsClearConvergenceStats > 0 and nIterationsCheckConvergence > fMCMCNIterationsClearConvergenceStats )
			nIterationsCheckConvergence = fMCMCNIterationsClearConvergenceStats;

		fMCMCNTrials = fMCMCNIterationsClearConvergenceStats;
	
		fMCMCNIterationsConvergenceGlobal = -1;
		
		if (fMCMCCurrentIteration >= (int)fMCMCNIterationsPreRunMax) {
      BCLog::OutWarning(" Convergence never checked !");
      BCLog::OutWarning("   Increase maximum number of iterations in the pre-run /MCMCSetNIterationsMax()/");
		}

		fMCMCPhase = BCEngineMCMC::MCMCPreRunConvergenceCheck;
		while ( fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax and fMCMCNIterationsConvergenceGlobal < 0 ) {
		
			// Clear information (also on first iteration)
			if ( fMCMCNTrials == fMCMCNIterationsClearConvergenceStats ) {
				fMCMCxMean = fMCMCx;
				fMCMCxVar.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
				fMCMCprobMean = fMCMCprob;
				fMCMCprobVar.assign(fMCMCNChains,0);
				fMCMCNTrials = 0;
			}

			MCMCGetNewPointMetropolis();
			++fMCMCNTrials;
			++NTrialsForEff;

			if (fMCMCFlagWritePreRunToFile) {
				EvaluateObservables();
				MCMCInChainWriteChains();
			}

			MCMCInChainUpdateStatistics();

			if ( fMCMCNTrials % nIterationsCheckConvergence != 0 and fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax )
				continue;

			MCMCInChainTestConvergenceAllChains();

			if ( fMCMCNIterationsConvergenceGlobal > 0 )
				continue;

			BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));
			BCLog::OutDetail(Form("       - %-*s : R-Value",fParameters.MaxNameLength(),"Parameter"));
				
			for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar) {

				if ( GetParameter(ipar)->Fixed() )
					continue;

				if( fMCMCRValueParameters[ipar]-1 < fMCMCRValueParametersCriterion )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else if ( fMCMCRValueParameters.at(ipar) != std::numeric_limits<double>::max() )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f <--",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else
					BCLog::OutDetail(TString::Format("         %-*s :  MAX_DOUBLE <--",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data()));
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
				BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, but all scales are adjusted.", fMCMCNChains, fMCMCNIterationsPreRunMax));
			else
				BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, and could not adjust all scales.", fMCMCNChains, fMCMCNIterationsPreRunMax));
	}

		
	//////////////////////////////////////////////////
	// Run until pre-run iteration min criteria met
	if ( fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin ) {

		unsigned N = fMCMCNIterationsPreRunMin - fMCMCCurrentIteration;
		unsigned nwrite = UpdateFrequency(N);

		BCLog::OutDetail(Form(" Current iteration (%d) is below minimum for pre-run. Running %d more iterations.",fMCMCCurrentIteration,N));

		fMCMCPhase = BCEngineMCMC::MCMCPreRunFulfillMinimum;
		unsigned n = 0;
		while ( fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin ) {
			MCMCGetNewPointMetropolis();
			++NTrialsForEff;
			++n;

			if (fMCMCFlagWritePreRunToFile) {
				EvaluateObservables();
				MCMCInChainWriteChains();
			}

			if ( n % nwrite == 0) {
				BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, 100.*n/N));
				if (fMCMCFlagWritePreRunToFile and fMCMCTree)
					fMCMCTree -> AutoSave("SaveSelf");
			}
		}
		BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations in pre-run.", fMCMCCurrentIteration));
	}

	// print scale factors and efficiencies
	std::vector<double> scalefactors (GetNParameters(),0);
	std::vector<double> efficiencies (GetNParameters(),0);

	BCLog::OutDetail(Form(" --> Average scale factors and efficiencies (measured in %d iterations):",(NTrialsForEff==0) ? fMCMCNIterationsEfficiencyCheck : NTrialsForEff));
	BCLog::OutDetail(Form("       - %-*s : Scale factor    Efficiency",fParameters.MaxNameLength(),"Parameter"));
	for (unsigned i = 0; i < GetNParameters(); ++i) {
		if (GetParameter(i)->Fixed())
			continue;
		for (unsigned j = 0; j < fMCMCNChains; ++j) {
			efficiencies[i] += ( NTrialsForEff==0 ) ? fMCMCEfficiencies[j][i] / fMCMCNChains : 1.*fMCMCNTrialsTrue[j][i]/NTrialsForEff/fMCMCNChains;
			scalefactors[i] += fMCMCTrialFunctionScaleFactor[j][i] / fMCMCNChains;
		}
		BCLog::OutDetail(Form("         %-*s :     %6.02f %%        %4.1f %%",fParameters.MaxNameLength(),GetParameter(i)->GetName().data(), 100.*scalefactors[i], 100.*efficiencies[i]));
	}

	// reset current iteration
	fMCMCCurrentIteration = -1;

	// reset current chain
	fMCMCCurrentChain = -1;

	if (fMCMCFlagWritePreRunToFile)
		UpdateParameterTree();

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
	if (fMCMCFlagPreRun) {
		MCMCMetropolisPreRun();
		if (!fMCMCFlagWritePreRunToFile and fMCMCFlagWriteChainToFile)
			InitializeMarkovChainTree();
	}
	else {
		BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Not running prerun. This can cause trouble if the data have changed.");
		if (fMCMCFlagWriteChainToFile)
			InitializeMarkovChainTree();
	}

	// print to screen
	BCLog::OutSummary( "Run Metropolis MCMC ...");

	// reset run statistics
	MCMCResetRunStatistics();

	// set phase and cycle number
	fMCMCPhase = BCEngineMCMC::MCMCMainRun;

	BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

	unsigned nwrite = UpdateFrequency(fMCMCNIterationsRun);

	// start the run
	fMCMCCurrentIteration = 0;
	while ( fMCMCCurrentIteration < (int)fMCMCNIterationsRun ) {

		MCMCGetNewPointMetropolis();

		if ( fMCMCCurrentIteration % nwrite == 0 ) {
			BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration)/(double)fMCMCNIterationsRun*100.));
			if (fMCMCFlagWriteChainToFile and fMCMCTree)
				fMCMCTree -> AutoSave("SaveSelf");
		}
		
		if (fMCMCCurrentIteration % fMCMCNLag != 0) // apply lag
			continue;
			
		MCMCIterationInterface();		// user action (overloadable)
		
		EvaluateObservables();

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
	for (unsigned i = 0; i < GetNParameters(); ++i) {
		if (GetParameter(i)->Fixed())
			continue;
		BCLog::OutDetail(Form("         %-*s :     %4.1f %%",fParameters.MaxNameLength(),GetParameter(i)->GetName().data(), 100.*efficiencies[i]));
	}

	if (fMCMCFlagWriteChainToFile)
		UpdateParameterTree();

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

	if (fMCMCFlagWriteChainToFile and fMCMCOutputFileAutoclose)
		MCMCCloseOutputFile();

	return 1;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCResetRunStatistics()
{
	fMCMCNTrials     = 0;

	fMCMCNIterations.assign(fMCMCNChains,0);
	fMCMCprobMean.assign(fMCMCNChains,0);
	fMCMCprobVar.assign(fMCMCNChains,0);
	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(GetNParameters(),0));

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
int BCEngineMCMC::AddObservable(const char * name, double min, double max, const char * latexname)
{
	return AddObservable(new BCObservable(name, min, max, latexname));
}

// --------------------------------------------------------
int BCEngineMCMC::AddObservable(BCObservable * obs)
{
	return fObservables.Add(obs);
}

// --------------------------------------------------------
void BCEngineMCMC::EvaluateObservables() {
	for (unsigned i = 0; i < fMCMCNChains; ++i ) {
		CalculateObservables(fMCMCx[i]);
		for (unsigned j = 0; j < GetNObservables(); ++j)
			fMCMCObservables[i][j] = GetObservable(j) -> Value();
	}
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

	fMCMCNTrialsTrue.assign(fMCMCNChains,std::vector<int>(GetNParameters(), 0));
	fMCMCNTrials = 0;
	fMCMCxMax.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
	fMCMCxMean.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
	fMCMCxVar.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));

	fMCMCRValueParameters.assign(GetNParameters(), 0);

	SyncThreadStorage();

   if (fMCMCTrialFunctionScaleFactorStart.size() == 0 or fMCMCTrialFunctionScaleFactorStart.size()!=GetNParameters())
		 fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,std::vector<double>(GetNParameters(), 1.0));
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
				 if (fMCMCInitialPosition[j].size() == GetNParameters()) {
					 for (unsigned i = 0; i < GetNParameters(); ++i)
						 if (!GetParameter(i)->IsValid(fMCMCInitialPosition[j][i])) {
							 BCLog::OutError(TString::Format("BCEngine::MCMCInitialize : Initial position of parameter %u (\"%s\") is out of boundaries in chain %u.",i,GetParameter(i)->GetName().data(),j).Data());
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
			 for (unsigned i = 0; i < GetNParameters(); ++i)
				 if (GetParameter(i)->Fixed())
					 for (unsigned j = 0; j < fMCMCNChains; ++j)
						 fMCMCx[j][i] = GetParameter(i)->GetFixedValue();
		 }
	 }
	 
	 if (fMCMCFlagInitialPosition == 0) { // center of the interval
		 fMCMCx.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
		 for (unsigned j = 0; j < fMCMCNChains; ++j)
			 for (unsigned i = 0; i < GetNParameters(); ++i)
				 fMCMCx[j][i] = (GetParameter(i)->Fixed()) ? GetParameter(i)->GetFixedValue() : 0.5 * (GetParameter(i)->GetLowerLimit() + GetParameter(i)->GetUpperLimit());

	 } else { // random number (default)
		 fMCMCx.assign(fMCMCNChains,std::vector<double>(GetNParameters(),0));
		 for (unsigned j = 0; j < fMCMCNChains; ++j)
			 for (unsigned i = 0; i < GetNParameters(); ++i)
				 fMCMCx[j][i] = (GetParameter(i)->Fixed()) ? GetParameter(i)->GetFixedValue() : GetParameter(i)->GetLowerLimit() + fMCMCThreadLocalStorage[j].rng->Rndm() * (GetParameter(i)->GetUpperLimit()-GetParameter(i)->GetLowerLimit());
	 }
	 
	 // initialize user-defined observables
	 fMCMCObservables.assign(fMCMCNChains,std::vector<double>(GetNObservables(),0));

   // define 1-dimensional histograms for marginalization
   bool fillAny = false;
	 fH1Marginalized.assign(GetNParameters()+GetNObservables(),NULL);

   for(unsigned i = 0; i < GetNParameters(); ++i)
		 if (GetParameter(i)->FillHistograms() && !GetParameter(i)->Fixed()) {
			 fH1Marginalized[i] = GetParameter(i) -> CreateH1(TString::Format("h1_%s_parameter_%i", GetSafeName().data() ,i));
			 fH1Marginalized[i] -> SetStats(kFALSE);
			 fillAny = true;
		 }
   for(unsigned i = 0; i < GetNObservables(); ++i)
		 if (GetObservable(i)->FillHistograms()) {
			 fH1Marginalized[i+GetNParameters()] = GetObservable(i) -> CreateH1(TString::Format("h1_%s_observable_%i", GetSafeName().data() ,i));
			 fH1Marginalized[i+GetNParameters()] -> SetStats(kFALSE);
			 fillAny = true;
		 }
	 
   // if filling no histograms, set H1 vector to zero size, implies no 2D histograms either
   if (!fillAny)
		 fH1Marginalized.clear();
   else {
		 fH2Marginalized.assign(GetNParameters()+GetNObservables(),std::vector<TH2D*>(GetNParameters()+GetNObservables(),0));
		 // define 2-dimensional histograms for marginalization
		 for(unsigned i = 0; i < GetNParameters(); ++i)
			 if (GetParameter(i)->FillHistograms() and !GetParameter(i)->Fixed()) {
				 // paramater vs parameter
				 for (unsigned j = i + 1; j < GetNParameters(); ++j)
					 if (GetParameter(j)->FillHistograms() and !GetParameter(j)->Fixed()) {
						 fH2Marginalized[i][j] = GetParameter(i) -> CreateH2(Form("h2_%s_parameters_%i_vs_%i", GetSafeName().data(), i, j), GetParameter(j));
						 fH2Marginalized[i][j] -> SetStats(kFALSE);
					 }
				 // user-defined observable vs parameter
				 for (unsigned j = 0; j < GetNObservables(); ++j)
					 if (GetObservable(j)->FillHistograms()) {
						 fH2Marginalized[i][j+GetNParameters()] = GetParameter(i) -> CreateH2(Form("h2_%s_par_%i_vs_obs_%i", GetSafeName().data(), i, j), GetObservable(j));
						 fH2Marginalized[i][j+GetNParameters()] -> SetStats(kFALSE);
					 }
			 }
		 // user-defined observable vs user-defined observable
		 for(unsigned i = 0; i < GetNObservables(); ++i)
			 if (GetObservable(i)->FillHistograms())
				 for (unsigned j = i + 1; j < GetNObservables(); ++j)
					 if (GetObservable(j)->FillHistograms()) {
						 fH2Marginalized[i+GetNParameters()][j+GetNParameters()] = GetObservable(i) -> CreateH2(Form("h2_%s_observables_%i_vs_%i", GetSafeName().data(), i, j), GetObservable(j));
						 fH2Marginalized[i+GetNParameters()][j+GetNParameters()] -> SetStats(kFALSE);
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
      GetParameter(i)->PrintSummary();

	 if (GetNObservables()>0) {
		 BCLog::OutSummary(Form("Number of observables : %u", GetNObservables()));
		 BCLog::OutSummary("Observables:");

		 // observable summary
		 for (unsigned i = 0; i < GetNObservables(); i++)
			 GetObservable(i)->PrintSummary();
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
	if ( P.size() != GetNParameters() )
		return;

	for (unsigned i = 0; i < GetNParameters(); ++i)
		output(TString::Format("          %-*s :   % .*g", fParameters.MaxNameLength(), GetParameter(i)->GetName().data(),GetParameter(i)->GetPrecision(),P[i]));
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintAllMarginalized1D(const char * filebase) {
	if (fH1Marginalized.size() == 0) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : No marginalized distributions stored.");
		return 0;
	}

	int nh = 0;
	for (unsigned i = 0; i < fH1Marginalized.size(); ++i) {
		if (!MarginalizedHistogramExists(i))
			continue;
		GetMarginalized(i) -> Print(Form("%s_1D_%s.pdf", filebase, GetVariable(i)->GetSafeName().data()));
		nh++;
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
	for (unsigned i = 0; i<fH2Marginalized.size(); ++i)
		for (unsigned j = 1; j<fH2Marginalized[i].size(); ++i) {
			if (MarginalizedHistogramExists(i,j))
				continue;
			GetMarginalized(i,j) -> Print(Form("%s_2D_%s_%s",filebase,GetVariable(i)->GetSafeName().data(),GetVariable(j)->GetSafeName().data()));
			nh++;
		}

	return nh;
}

// ---------------------------------------------------------
int BCEngineMCMC::PrintAllMarginalized(std::string filename, std::string options1d, std::string options2d, unsigned int hdiv, unsigned int vdiv) {
	if (fH1Marginalized.empty() and fH2Marginalized.empty()) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : Marginalized distributions not stored.");
		return 0;
	}

	 // Find nonempty H1's
   std::vector<BCH1D *> h1;
   for (unsigned i = 0; i < GetNVariables(); ++i)
		 if ( MarginalizedHistogramExists(i) ) {
			 if (GetMarginalizedHistogram(i)->Integral()==0) { // histogram was never filled in range
				 BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 1D Marginalized histogram for \"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data()));
				 continue;
			 }
			 h1.push_back(GetMarginalized(i));
			 if (!h1.back()) // BCH1D doesn't exist
				 h1.pop_back();
		 }

	 // Find nonempty H2's
   std::vector<BCH2D *> h2;
	 // fill h2 vector in order (par vs par; obs vs par; obs vs obs)

	 // parameter vs parameter
	 for (unsigned i = 0; i<GetNParameters(); ++i)
		 for (unsigned j = i+1; j<GetNParameters(); ++j)
			 if ( MarginalizedHistogramExists(i,j) ) {
				 if (fH2Marginalized[i][j]->Integral()==0) { // histogram was never filled in range
					 BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 2D Marginalized histogram for \"%s\":\"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data(),GetVariable(i)->GetName().data()));
					 continue;
				 }
				 h2.push_back(GetMarginalized(i,j));
				 if (!h2.back()) // BCH2D doesn't exist
					 h2.pop_back();
			 }
	 
		 // user-defined observable vs parameter
	 for (unsigned i = 0; i<GetNParameters(); ++i)
		 for (unsigned j = GetNParameters(); j<GetNVariables(); ++j)
			 if ( MarginalizedHistogramExists(i,j) ) {
				 if (GetMarginalizedHistogram(i,j)->Integral()==0) { // histogram was never filled in range
					 BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 2D Marginalized histogram for \"%s\":\"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data(),GetVariable(i)->GetName().data()));
					 continue;
				 }
				 h2.push_back(GetMarginalized(i,j));
				 // check if histogram exists
				 if (!h2.back())
					 h2.pop_back();
			 }

	 // user-defined observable vs user-defined observable
	 for (unsigned i = GetNParameters(); i<GetNVariables(); ++i)
		 for (unsigned j = i+1; j<GetNVariables(); ++j)
			 if ( MarginalizedHistogramExists(i,j) ) {
				 if (GetMarginalizedHistogram(i,j)->Integral()==0) { // histogram was never filled in range
					 BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 2D Marginalized histogram for \"%s\":\"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data(),GetVariable(i)->GetName().data()));
					 continue;
				 }
				 h2.push_back(GetMarginalized(i,j));
				 // check if histogram exists
				 if (!h2.back())
					 h2.pop_back();
			 }

	 if (h1.empty() and h2.empty()) {
		 BCLog::OutWarning("BCEngineMCMC::PrintAllMarginalized : No marginalizations to print");
		 return 0;
	 }

	 // set default drawing options
	 if (options1d.empty())
		 options1d = "BTsiB3CS1D0pdf0Lmeanmode";
	 if (options2d.empty())
		 options2d = "BTfB3CS1meangmode";

   // if file has no extension or if it's not ".pdf" or ".ps", make it ".pdf"
   if ( (filename.find_last_of(".") == std::string::npos) or
				((filename.substr(filename.find_last_of(".")) != ".pdf") and	(filename.substr(filename.find_last_of(".")) != ".ps")))
		 filename += ".pdf";

	 int c_width  = 297*4;
	 int c_height = 210*4;
	 if (hdiv < vdiv)
		 std::swap(c_width,c_height);

   const unsigned nplots = h1.size() + h2.size();

   // give out warning if too many plots
   BCLog::OutSummary(Form("Printing all marginalized distributions (%lu x 1D + %lu x 2D = %u) into file %s", h1.size(), h2.size(), nplots, filename.c_str()));
   if (nplots > 100)
		 BCLog::OutDetail("This can take a while...");
	 
   // setup the canvas and file
   TCanvas c("c", "canvas", c_width, c_height);
   c.Divide(hdiv, vdiv);
	 c.Print((filename+"[").data());

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
	c_par -> cd();

	if (npar<=0) // all parameters on one page, all user-defined observables on the next
		npar = std::max<int> (GetNParameters(),GetNObservables());

	int return_val = 1;

	// parameters first
	for (unsigned i = 0; i<GetNParameters(); i += npar)
		if (DrawParameterPlot(i,std::min<int>(npar,GetNParameters()-i), interval_content, quantiles)) {
			c_par->Print(filename);
			c_par->Clear();
		}
		else
			return_val = 0;
	
	// then user-defined observables
	for (unsigned i = GetNParameters(); i<GetNVariables(); i += npar)
		if(DrawParameterPlot(i,std::min<int>(npar,GetNVariables()-i), interval_content, quantiles)) {
			c_par -> Print(filename);
			c_par -> Clear();
		}
		else
			return_val = 0;

	c_par -> Print(Form("%s]",filename));
	return return_val;
}

// ---------------------------------------------------------
int BCEngineMCMC::DrawParameterPlot(unsigned i0, unsigned npar, double interval_content, std::vector<double> quantiles) {

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
	TLegend * legend = new TLegend();
	legend -> SetBorderSize(0);
	legend -> SetFillColor(kWhite);
	legend -> SetNColumns(2);
	legend -> SetTextAlign(12);
	legend -> SetTextFont(62);
	legend -> SetTextSize(0.02*gPad->GetWNDC());

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
		graph_mean->SetMarkerSize(1*gPad->GetWNDC());
		graph_mean->Draw("SAMEP");

		legend -> AddEntry(graph_mean, "Mean and RMS", "LEP");
		legend -> AddEntry(graph_intervals, TString::Format("Smallest %.0f%% interval and local mode",100.*interval_content), "FL");
	}

	// Global Modes
	if (!x_i_bf.empty()) {
		TGraph * graph_mode = new TGraph(x_i_bf.size(), x_i_bf.data(), global_mode.data());
		graph_mode->SetMarkerColor(kRed);
		graph_mode->SetMarkerStyle(20);
		graph_mode->SetMarkerSize(1*gPad->GetWNDC());
		graph_mode->Draw("SAMEP");
		legend->AddEntry(graph_mode, "Global mode", "P");
	}

	gPad->SetTopMargin(0.02);

	// place legend on top of histogram
	legend->SetX1NDC(gPad->GetLeftMargin());
	legend->SetX2NDC(1. - gPad->GetRightMargin());
	double y1 = gPad->GetTopMargin() + legend->GetTextSize()*legend->GetNRows();
	legend->SetY1NDC(1.-y1);
	legend->SetY2NDC(1. - gPad->GetTopMargin());
	legend->Draw("SAME");

	gPad -> SetTopMargin(y1+0.01);

	gPad->RedrawAxis();

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

				hh -> SetStats(false);
				hh -> GetXaxis() -> SetLabelSize(0);
				hh -> GetYaxis() -> SetLabelSize(0);
				hh -> GetXaxis() -> SetTitleSize(0);
				hh -> GetYaxis() -> SetTitleSize(0);

				if (bh1)
					bh1 -> Draw("BTsiB3CS1D0");
				else
					bh2 -> Draw("BTfB3CS1nL");

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
		fMCMCThreadLocalStorage.push_back(MCMCThreadLocalStorage(GetNParameters()));
		// each chains gets a different seed. fRandom always returns same seed after the fixing done above
		fMCMCThreadLocalStorage.back().rng->SetSeed(fRandom->GetSeed() + fMCMCThreadLocalStorage.size());
	}

	// remove storage until equal to number of chain
	while (fMCMCThreadLocalStorage.size() > fMCMCNChains)
		fMCMCThreadLocalStorage.pop_back();

	// update parameter size for each chain
	for (unsigned i = 0 ; i < fMCMCThreadLocalStorage.size(); ++i)
		fMCMCThreadLocalStorage[i].xLocal.assign(GetNParameters(), 0.0);
}
