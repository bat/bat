/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCEngineMCMC.h"

#include "BCAux.h"
#include "BCMath.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCVariable.h"

#include <TH1.h>
#include <TH2.h>
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
#include <TF1.h>
#include <TObject.h>
#include <TKey.h>
#include <TList.h>
#include <TDecompChol.h>
#include <TVectorD.h>
#include <TError.h>

#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>

#include <typeinfo>

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const char * name)
	: fMCMCFlagWriteChainToFile(false)
	, fMCMCFlagWritePreRunToFile(false)
	, fMCMCOutputFile(0)
	, fMCMCOutputFilename("")
	, fMCMCOutputFileOption("")
	, fMCMCScaleFactorLowerLimit(0)
	, fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max())
	, fMCMCAutoSetTrialFunctionScaleFactors(true)
	// , fMultivariateProposalFunctionTuningScheduleParameter(0.5)
	, fMultivariateProposalFunctionEpsilon(1e-3)
	, fMultivariateProposalFunctionScaleMultiplier(1.5)
	, fMCMCFlagPreRun(true)
	, fMCMCInitialPositionExpansionFactor(1.1)
	, fMCMCEfficiencyMin(0.15)
	, fMCMCEfficiencyMax(0.50)
	, fMCMCFlagInitialPosition(BCEngineMCMC::kMCMCInitRandomUniform)
	, fMCMCMultivariateProposalFunction(false)
	, fMCMCPhase(BCEngineMCMC::MCMCUnsetPhase)
	, fCorrectRValueForSamplingVariability(false)
	, fMCMCRValueCriterion(0.1)
	, fMCMCRValueParametersCriterion(0.1)
	, fRandom(new TRandom3())
	, fMCMCTree(0)
	, fMCMCTreeLoaded(false)
	, fMCMCTreeReuseObservables(true)
	, fParameterTree(0)
	, fBCH1DdrawingOptions(new BCH1D)
	, fBCH2DdrawingOptions(new BCH2D)
	, fRescaleHistogramRangesAfterPreRun(false)
	, fHistogramRescalePadding(0.1)
{
	SetName(name);
	MCMCSetPrecision(BCEngineMCMC::kMedium);
	MCMCSetRandomSeed(0);
	PartnerUp(&fParameters,&fObservables);
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(std::string filename, std::string name, bool reuseObservables)
	: fMCMCFlagWriteChainToFile(false)
	, fMCMCFlagWritePreRunToFile(false)
	, fMCMCOutputFile(0)
	, fMCMCOutputFilename("")
	, fMCMCOutputFileOption("")
	, fMCMCScaleFactorLowerLimit(0)
	, fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max())
	, fMCMCAutoSetTrialFunctionScaleFactors(true)
	// , fMultivariateProposalFunctionTuningScheduleParameter(0.5)
	, fMultivariateProposalFunctionEpsilon(1e-3)
	, fMultivariateProposalFunctionScaleMultiplier(1.5)
	, fMCMCFlagPreRun(true)
	, fMCMCInitialPositionExpansionFactor(1.1)
	, fMCMCEfficiencyMin(0.15)
	, fMCMCEfficiencyMax(0.50)
	, fMCMCFlagInitialPosition(BCEngineMCMC::kMCMCInitRandomUniform)
	, fMCMCMultivariateProposalFunction(false)
	, fMCMCPhase(BCEngineMCMC::MCMCUnsetPhase)
	, fCorrectRValueForSamplingVariability(false)
	, fMCMCRValueCriterion(0.1)
	, fMCMCRValueParametersCriterion(0.1)
	, fRandom(new TRandom3())
	, fMCMCTree(0)
	, fMCMCTreeLoaded(false)
	, fMCMCTreeReuseObservables(true)
	, fParameterTree(0)
	, fBCH1DdrawingOptions(new BCH1D)
	, fBCH2DdrawingOptions(new BCH2D)
	, fRescaleHistogramRangesAfterPreRun(false)
	, fHistogramRescalePadding(0.1)
{
	SetName(name);
	MCMCSetPrecision(BCEngineMCMC::kMedium);
	MCMCSetRandomSeed(0);
	LoadMCMC(filename,"","",reuseObservables);
	PartnerUp(&fParameters,&fObservables);
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC & other)
	: fBCH1DdrawingOptions(new BCH1D)
	, fBCH2DdrawingOptions(new BCH2D)
{
	Copy(other);
}

// ---------------------------------------------------------
BCEngineMCMC::~BCEngineMCMC() {
	MCMCCloseOutputFile();

	// delete random number generator
	delete fRandom;
	
	// delete 1-d marginalized distributions
	for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
		delete fH1Marginalized[i];
	
	// delete 2-d marginalized distributions
	for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
		for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
			delete fH2Marginalized[i][j];

	delete fBCH1DdrawingOptions;
	delete fBCH2DdrawingOptions;
}


// ---------------------------------------------------------
void BCEngineMCMC::SetName(const char * name) {
	fName = name;
	fSafeName = name;
	fSafeName.erase(std::remove_if(fSafeName.begin(),fSafeName.end(),::isspace),fSafeName.end());
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCSetPrecision(BCEngineMCMC::Precision precision) {

	// all precision levels want a pre-run:
	fMCMCFlagPreRun = true;

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
	fMCMCFlagPreRun                       = other -> MCMCGetFlagPreRun();
	fMCMCMultivariateProposalFunction     = other -> MCMCGetMultivariateProposalFunction();
	// fMultivariateProposalFunctionTuningScheduleParameter = other -> MCMCGetMultivariateProposalFunctionTuningScheduleParameter();
	fMultivariateProposalFunctionEpsilon  = other -> MCMCGetMultivariateProposalFunctionEpsilon();
	fMultivariateProposalFunctionScaleMultiplier = other -> MCMCGetMultivariateProposalFunctionScaleMultiplier();
}

// ---------------------------------------------------------
void BCEngineMCMC::Copy(const BCEngineMCMC & other) {
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
	fMCMCFlagWriteChainToFile                 = other.fMCMCFlagWriteChainToFile;
	fMCMCFlagWritePreRunToFile                = other.fMCMCFlagWritePreRunToFile;
	fMCMCTrialFunctionScaleFactor             = other.fMCMCTrialFunctionScaleFactor;
	fMCMCTrialFunctionScaleFactorStart        = other.fMCMCTrialFunctionScaleFactorStart;
	fMCMCAutoSetTrialFunctionScaleFactors     = other.fMCMCAutoSetTrialFunctionScaleFactors;
	fMCMCFlagPreRun                           = other.fMCMCFlagPreRun;
	fMCMCFlagRun                              = other.fMCMCFlagRun;
	fMCMCInitialPosition                      = other.fMCMCInitialPosition;
	fMCMCInitialPositionExpansionFactor       = other.fMCMCInitialPositionExpansionFactor;
	fMCMCEfficiencyMin                        = other.fMCMCEfficiencyMin;
	fMCMCEfficiencyMax                        = other.fMCMCEfficiencyMax;
	fMCMCScaleFactorLowerLimit                = other.fMCMCScaleFactorLowerLimit;
	fMCMCScaleFactorUpperLimit                = other.fMCMCScaleFactorUpperLimit;
	fMCMCFlagInitialPosition                  = other.fMCMCFlagInitialPosition;
	fMCMCMultivariateProposalFunction         = other.fMCMCMultivariateProposalFunction;
	fMCMCPhase                                = other.fMCMCPhase;
	fMCMCx                                    = other.fMCMCx;
	fMCMCObservables                          = other.fMCMCObservables;
	fMCMCStatistics                           = other.fMCMCStatistics;
	fMCMCStatistics_AllChains                 = other.fMCMCStatistics_AllChains;
	fMCMCprob                                 = other.fMCMCprob;
	fMCMCLogLikelihood                        = other.fMCMCLogLikelihood;
	fMCMCLogLikelihood_Provisional            = other.fMCMCLogLikelihood_Provisional;
	fMCMCLogPrior                             = other.fMCMCLogPrior;
	fMCMCLogPrior_Provisional                 = other.fMCMCLogPrior_Provisional;
	fCorrectRValueForSamplingVariability      = other.fCorrectRValueForSamplingVariability;
	fMCMCRValueCriterion                      = other.fMCMCRValueCriterion ;
	fMCMCRValueParametersCriterion            = other.fMCMCRValueParametersCriterion;
	fMCMCRValue                               = other.fMCMCRValue;
	fMCMCRValueParameters                     = other.fMCMCRValueParameters;
	fRandom                                   = (other.fRandom) ? new TRandom3(*other.fRandom) : NULL;
	fMCMCThreadLocalStorage                   = other.fMCMCThreadLocalStorage;
	fRescaleHistogramRangesAfterPreRun        = other.fRescaleHistogramRangesAfterPreRun;
	fHistogramRescalePadding                  = other.fHistogramRescalePadding;

	// multivariate proposal function shtuff
	// fMultivariateProposalFunctionTuningScheduleParameter = other.fMultivariateProposalFunctionTuningScheduleParameter;
	fMultivariateProposalFunctionEpsilon = other.fMultivariateProposalFunctionEpsilon;
	fMultivariateProposalFunctionScaleMultiplier = other.fMultivariateProposalFunctionScaleMultiplier;
	fMultivariateProposalFunctionCovariance = other.fMultivariateProposalFunctionCovariance;
	fMultivariateProposalFunctionCholeskyDecomposition = other.fMultivariateProposalFunctionCholeskyDecomposition;
	fMultivariateProposalFunctionTuningSteps = other.fMultivariateProposalFunctionTuningSteps;

	fParameters = other.fParameters;
	fObservables = other.fObservables;
	PartnerUp(&fParameters,&fObservables);

	// clear existing histograms
	for (unsigned i=0; i<fH1Marginalized.size(); ++i)
		delete fH1Marginalized[i];
	for (unsigned i=0; i<fH2Marginalized.size(); ++i)
		for (unsigned j=0; j<fH2Marginalized[i].size(); ++j)
			delete fH2Marginalized[i][j];

	fH1Marginalized = std::vector<TH1*>(other.fH1Marginalized.size(),NULL);
	for (unsigned i = 0; i < other.fH1Marginalized.size(); ++i)
		if (other.fH1Marginalized[i])
			fH1Marginalized[i] = dynamic_cast<TH1*>(other.fH1Marginalized[i]->Clone());
	// fH1Marginalized[i] = new TH1D(*(other.fH1Marginalized[i]));
	 
	if (!other.fH2Marginalized.empty() and !other.fH2Marginalized.front().empty()) {
		fH2Marginalized = std::vector<std::vector<TH2*> > (other.fH2Marginalized.size(), std::vector<TH2*>(other.fH2Marginalized.front().size(),NULL));
		for (unsigned i = 0; i < other.fH2Marginalized.size(); ++i) {
			fH2Marginalized[i].assign(other.fH2Marginalized[i].size(),NULL);
			for (unsigned j = 0; j < other.fH2Marginalized[i].size(); ++j)
				if (other.fH2Marginalized[i][j])
					fH2Marginalized[i][j] = dynamic_cast<TH2*>(other.fH2Marginalized[i][j]->Clone());
			// fH2Marginalized[i][j] = new TH2D(*(other.fH2Marginalized[i][j]));
		}
	} else
		fH2Marginalized = std::vector<std::vector<TH2*> >();
	 
	fMCMCTree       = 0;
	fMCMCTreeLoaded = false;
	fMCMCTreeReuseObservables = true;
	fParameterTree  = 0;
	fMCMCOutputFile = 0;
	fMCMCOutputFilename      = other.fMCMCOutputFilename;
	fMCMCOutputFileOption    = other.fMCMCOutputFileOption;

	fLocalModes            = other.fLocalModes;
	fBCH1DdrawingOptions -> CopyOptions(*(other.fBCH1DdrawingOptions));
	fBCH2DdrawingOptions -> CopyOptions(*(other.fBCH2DdrawingOptions));
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChain(bool flag) {
	if (flag)
		BCLog::OutError("BCEngineMCMC::WriteMarkovChain: To turn on output use WriteMarkovChain(filename,option).");
	fMCMCFlagWriteChainToFile = false;
	fMCMCFlagWritePreRunToFile = false;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChain(std::string filename, std::string option) {
	if (filename.empty()) {
		BCLog::OutError("BCEngineMCMC::WriteMarkovChain: You must specify the filename when turning on Markov chain output.");
		return WriteMarkovChain(false);
	}
	fMCMCOutputFilename = filename;
	fMCMCOutputFileOption = option;
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
TH1 * const BCEngineMCMC::GetMarginalizedHistogram(unsigned index) const {
	if ( index >= fH1Marginalized.size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.",index));
		return 0;
	}

	if (fH1Marginalized[index])
		return fH1Marginalized[index];
	
	// else output warning
	if ( index<GetNVariables() ) // Marginalization of model parameter
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for %s %s", GetVariable(index)->GetPrefix().data(),GetVariable(index)->GetName().data()));
	return 0;
}

// --------------------------------------------------------
TH2 * const BCEngineMCMC::GetMarginalizedHistogram(unsigned i, unsigned j) const {
	if (i == j) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Called with identical indices %u.", i));
		return NULL;
	}
	
	if ( i >= fH2Marginalized.size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", i));
		return NULL;
	}
	if ( j >= fH2Marginalized[i].size() ) {
		BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", j));
		return NULL;
	}

	if (fH2Marginalized[i][j])
		return fH2Marginalized[i][j];

	if (i<GetNVariables() and j<GetNVariables())
		BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram : marginal distribution not stored for %s %s vs %s %s",
													 GetVariable(i)->GetPrefix().data(),GetVariable(i)->GetName().data(),
													 GetVariable(j)->GetPrefix().data(),GetVariable(j)->GetName().data()));
	return NULL;
}

// --------------------------------------------------------
BCH1D * BCEngineMCMC::GetMarginalized(unsigned index) const {
	TH1 * h = GetMarginalizedHistogram(index);

	if ( !h )
		return NULL;
	
	BCH1D * hprob = new BCH1D((TH1*)h);
	
	// set global mode if available
	if (index<GetGlobalMode().size())
		hprob -> SetGlobalMode(GetGlobalMode()[index]);
	
	return hprob;
}
	
// --------------------------------------------------------
BCH2D * BCEngineMCMC::GetMarginalized(unsigned i, unsigned j) const {
	TH2 * h = GetMarginalizedHistogram(i,j);

	if (!h)
		return NULL;

	BCH2D * hprob = new BCH2D(h);

	// set global mode if available
	if (i<GetGlobalMode().size() and j<GetGlobalMode().size())
		hprob -> SetGlobalMode(GetGlobalMode()[i],GetGlobalMode()[j]);
		
	return hprob;
}

// ---------------------------------------------------------
const std::vector<double> & BCEngineMCMC::GetLocalModes(bool force_recalculation) {
	if (fLocalModes.empty() or force_recalculation) {
		fLocalModes.clear();
		for (unsigned i=0; i<GetNVariables(); ++i)
			if (i<GetNParameters() and GetParameter(i)->Fixed())
				fLocalModes.push_back(GetParameter(i)->GetFixedValue());
			else if (fH1Marginalized[i]) {
				fLocalModes.push_back(fH1Marginalized[i]->GetBinCenter(fH1Marginalized[i]->GetMaximumBin()));
			} else {
				if (i<GetNParameters()) {
					BCLog::OutError("BCEngineMCMC::GetLocalModes : unfixed parameter is missing marginalized information. returning empty vector.");
					fLocalModes.clear();
					break;
				} else {
					BCLog::OutWarning("BCEngineMCMC::GetLocalModes : user-defined observable is missing marginalized information. returning only local modes of parameters.");
					fLocalModes.resize(GetNParameters());
					break;
				}
			}
	}
	return fLocalModes;
}

// ---------------------------------------------------------
std::vector<double> BCEngineMCMC::GetBestFitParameters() const {
	if (GetGlobalMode().empty() or GetGlobalMode().size()<GetNParameters())
		return std::vector<double>();
	return std::vector<double>(GetGlobalMode().begin(),GetGlobalMode().begin()+GetNParameters());
}

// ---------------------------------------------------------
std::vector<double> BCEngineMCMC::GetBestFitObservables() const {
	if (GetGlobalMode().empty() or GetGlobalMode().size()<GetNVariables())
		return std::vector<double>();
	return std::vector<double>(GetGlobalMode().begin()+GetNParameters(),GetGlobalMode().begin()+GetNVariables());
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetInitialPositions(const std::vector<double> & x0s) {
   fMCMCInitialPosition.clear();
	 for (std::vector<double>::const_iterator it=x0s.begin(); it+GetNParameters()<=x0s.end(); it += GetNParameters())
		 fMCMCInitialPosition.push_back(std::vector<double>(it,it+GetNParameters()));
	 MCMCSetFlagInitialPosition(BCEngineMCMC::kMCMCInitUserDefined);
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCSetRandomSeed(unsigned seed) {
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
void BCEngineMCMC::InitializeMarkovChainTree(bool replacetree, bool replacefile) {
	if (replacetree) {
		delete fMCMCTree;
		fMCMCTree = 0;
	}
	if (replacetree) {
		delete fParameterTree;
		fParameterTree = 0;
	}
	if (replacefile) {
		if (fMCMCOutputFile)
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
		fMCMCTree -> Branch("Chain",          &fMCMCTree_Chain,       "chain/i");
		fMCMCTree -> Branch("Iteration",      &fMCMCCurrentIteration, "iteration/i");
		fMCMCTree -> Branch("Phase",          &fMCMCPhase,            "phase/I");
		fMCMCTree -> Branch("LogProbability", &fMCMCTree_Prob,        "log(probability)/D");
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
		bool p_parameter, p_fill_1d, p_fill_2d, p_fixed;
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
		fParameterTree -> Branch("fill_1d",&p_fill_1d,"fill_1d/O");
		fParameterTree -> Branch("fill_2d",&p_fill_2d,"fill_2d/O");
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
			p_fill_1d    = GetVariable(i) -> FillH1();
			p_fill_2d    = GetVariable(i) -> FillH2();
			p_fixed      = p_parameter and GetParameter(i)->Fixed();
			p_fixedvalue = (p_parameter) ? GetParameter(i)->GetFixedValue() : 0;
			fParameterTree -> Fill();
		}
		fParameterTree -> AutoSave("SaveSelf");
		// fParameterTree -> ResetBranchAddresses();
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
			if (!fMCMCMultivariateProposalFunction)
				eff[j]   = (n<GetNParameters()) ? fMCMCStatistics[j].efficiency[n] : -1;
			else if (!fMCMCStatistics[j].efficiency.empty())
				eff[j] = fMCMCStatistics[j].efficiency.front();
			else 
				eff[j] = -1;
		}

		if (b_scale)
			b_scale -> Fill();

		b_eff -> Fill();
	}
	fParameterTree -> AutoSave("SaveSelf");
	// fParameterTree -> ResetBranchAddresses();
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidMCMCTree(TTree * tree, bool checkObservables) const {
	if (!tree)
		return false;

	if (!(tree->GetBranch("Chain")))
		return false;

	if (!(tree->GetBranch("Phase")))
		return false;

	// if (!(tree->GetBranch("Iteration")))
	// 	return false;
	// if (!(tree->GetBranch("LogProbability")))
	// 	return false;
	// if (!(tree->GetBranch("LogLikelihood")))
	// 	return false;
	// if (!(tree->GetBranch("LogPrior")))
	// 	return false;
	
	unsigned nvar = checkObservables ? GetNObservables() : GetNParameters();
	for (unsigned i=0; i<nvar; ++i)
		if (!(tree->GetBranch(GetVariable(i)->GetSafeName().data())))
			return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidParameterTree(TTree * tree) const {
	if (!tree)
		return false;

	if (!(tree->GetBranch("parameter")))
		return false;

	if (!(tree->GetBranch("index")))
		return false;
	
	if (!(tree->GetBranch("name")))
		return false;

	if (!(tree->GetBranch("lower_limit")))
		return false;

	if (!(tree->GetBranch("upper_limit")))
		return false;

	// if (!(tree->GetBranch("safe_name")))
	// 	return false;
	// if (!(tree->GetBranch("latex_name")))
	// 	return false;
	// if (!(tree->GetBranch("precission")))
	// 	return false;
	// if (!(tree->GetBranch("nbins")))
	// 	return false;
	// if (!(tree->GetBranch("fill_1d")))
	// 	return false;
	// if (!(tree->GetBranch("fill_2d")))
	// 	return false;
	// if (!(tree->GetBranch("fixed")))
	// 	return false;
	// if (!(tree->GetBranch("fixed_value")))
	// 	return false;
	// if (!(tree->GetBranch("nchain")))
	// 	return false;
	// if (!(tree->GetBranch("scale")))
	// 	return false;

	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::LoadParametersFromTree(TTree * partree, bool reuseObservables) {
	bool p_fill_1d = true;
	bool p_fill_2d = true;
	bool p_fixed = false;
	unsigned p_precision = 6;
	unsigned p_nbins = 100;
	char p_name[200];
	char p_latexname[200] = "";
	double p_lowerlimit;
	double p_upperlimit;
	double p_fixedvalue = 0;

	// absolutely necessary branches
	partree -> SetBranchAddress("name",p_name);
	partree -> SetBranchAddress("lower_limit",&p_lowerlimit);
	partree -> SetBranchAddress("upper_limit",&p_upperlimit);

	// not entirely necessary branches
	if (partree -> GetBranch("latex_name"))
		partree -> SetBranchAddress("latex_name",p_latexname);
	if (partree -> GetBranch("precision")) 
		partree -> SetBranchAddress("precision",&p_precision);
	if (partree -> GetBranch("nbins"))
		partree -> SetBranchAddress("nbins",&p_nbins);
	if (partree -> GetBranch("fill_1d"))
		partree -> SetBranchAddress("fill_1d",&p_fill_1d);
	if (partree -> GetBranch("fill_2d"))
		partree -> SetBranchAddress("fill_2d",&p_fill_2d);
	if (partree -> GetBranch("fixed"))
		partree -> SetBranchAddress("fixed",&p_fixed);
	if (partree -> GetBranch("fixed_value"))
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
		Par -> FillHistograms(p_fill_1d,p_fill_2d);
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
		Obs -> FillHistograms(p_fill_1d,p_fill_2d);
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
	partree -> BuildIndex("parameter","index");

	bool has_fixed = true;
	if (partree->GetBranch("fixed")) 
		partree -> SetBranchAddress("fixed",&p_fixed);
	else
		has_fixed = false;

	bool has_fixed_value = true;
	if (partree->GetBranch("fixed_value"))
		partree -> SetBranchAddress("fixed_value",&p_fixedvalue);
	else
		has_fixed_value = false;

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
		if ( has_fixed and GetParameter(i)->Fixed() != p_fixed ) {
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Fixed status of parameter \"%s\" does not match. Fixing it.",p_name));
			GetParameter(i) -> Fix(p_fixedvalue);
		}
		if ( has_fixed and GetParameter(i)->Fixed() and has_fixed_value and GetParameter(i)->GetFixedValue() != p_fixedvalue ) {
			BCLog::OutWarning(TString::Format("BCEngineMCMC::ParameterTreeMatchesModel : Fixed value of parameter \"%s\" does not match. Updating it.",p_name));
			GetParameter(i) -> Fix(p_fixedvalue);
		}
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

	// set model name if empty
	if (fName.empty()) {
		// check mcmcTreeName and parameterTreeName for default BAT name scheme [modelname]_mcmc/parameters:
		if ( mcmcTreeName.find_last_of("_") != std::string::npos and mcmcTreeName.substr(mcmcTreeName.find_last_of("_"))=="_mcmc"
				 and parameterTreeName.find_last_of("_") != std::string::npos and parameterTreeName.substr(parameterTreeName.find_last_of("_"))=="_parameters") {
			fName = mcmcTreeName.substr(0,mcmcTreeName.find_last_of("_"));
		}
		// else look through file for trees named according to BAT scheme
		else {
			TList * LoK = inputfile -> GetListOfKeys();
			std::vector<std::string> mcmc_names;
			std::vector<std::string> parameter_names;
			for (int i=0; i<LoK->GetEntries(); ++i) {
				TKey * k = (TKey*)(LoK->At(i));
				if (strcmp(k->GetClassName(),"TTree")!=0)
					continue;
				std::string treeName(k->GetName());
				if (treeName.find_last_of("_") == std::string::npos) 
					continue;
				if (treeName.substr(treeName.find_last_of("_"))=="_mcmc")
					mcmc_names.push_back(treeName.substr(0,treeName.find_last_of("_")));
				else if (treeName.substr(treeName.find_last_of("_"))=="_parameters")
					parameter_names.push_back(treeName.substr(0,treeName.find_last_of("_")));
			}

			std::vector<std::string> model_names;
			for (unsigned i=0; i<mcmc_names.size(); ++i)
				for (unsigned j=0; j<parameter_names.size(); ++j)
					if (mcmc_names[i] == parameter_names[j])
						model_names.push_back(mcmc_names[i]);

			if (model_names.empty()) {
				BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC : %s contains no matching MCMC and Parameter trees.",filename.data()));
				return false;
			}

			if (model_names.size()>1) {
				BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC : %s contains more than one model, please select one by providing a model name:",filename.data()));
				for (unsigned i=0; i<model_names.size(); ++i)
					BCLog::OutError(TString::Format("BCEngineMCMC::LoadMCMC : \"%s\"",model_names[i].data()));
				return false;
			}
			
			SetName(model_names[0]);
		}
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
	delete fMCMCTree;
	fMCMCTree = mcmcTree;

	fMCMCTreeLoaded = true;
	return true;
}

// --------------------------------------------------------
void BCEngineMCMC::Remarginalize(bool autorange) {
	if (!ValidMCMCTree(fMCMCTree,fMCMCTreeReuseObservables))
		return;

	fMCMCTree -> SetBranchAddress("Chain",          &fMCMCTree_Chain);
	fMCMCTree -> SetBranchAddress("Phase",          &fMCMCPhase);

	if (fMCMCTree -> GetBranch("LogProbability"))
		fMCMCTree -> SetBranchAddress("LogProbability", &fMCMCTree_Prob);

	bool has_iteration = true;
	if (fMCMCTree -> GetBranch("Iteration"))
			fMCMCTree -> SetBranchAddress("Iteration",      &fMCMCTree_Iteration);
	else
		has_iteration = false;

	if (fMCMCTree -> GetBranch("LogLikelihood"))
		fMCMCTree -> SetBranchAddress("LogLikelihood",&fMCMCTree_LogLikelihood);
	else
		fMCMCTree_LogLikelihood = -std::numeric_limits<double>::infinity();
	if (fMCMCTree -> GetBranch("LogPrior"))
		fMCMCTree -> SetBranchAddress("LogPrior",&fMCMCTree_LogPrior);
	else
		fMCMCTree_LogPrior = -std::numeric_limits<double>::infinity();
	

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

	if (autorange) {
		std::vector<double> XMin;
		std::vector<double> XMax;
		if (fMCMCStatistics_AllChains.n_samples > 0) {
			XMin = fMCMCStatistics_AllChains.minimum;
			XMax = fMCMCStatistics_AllChains.maximum;
		} else {
			MCMCInitialize();
			
			// find min and max
			XMin.assign(GetNVariables(),+std::numeric_limits<double>::infinity());
			XMax.assign(GetNVariables(),-std::numeric_limits<double>::infinity());
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
						if (fMCMCTree_Observables[i]<XMin[fMCMCTree_Parameters.size()+i])
							XMin[fMCMCTree_Parameters.size()+i] = fMCMCTree_Observables[i];
						if (fMCMCTree_Observables[i]>XMax[fMCMCTree_Parameters.size()+i])
							XMax[fMCMCTree_Parameters.size()+i] = fMCMCTree_Observables[i];
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
		}
		if (!XMin.empty() and !XMax.empty()) {
			// store mins and maxes, and change
			std::vector<double> xmin;
			std::vector<double> xmax;
			for (unsigned i=0; i<GetNVariables(); ++i) {
				xmin.push_back(GetVariable(i)->GetLowerLimit());
				xmax.push_back(GetVariable(i)->GetUpperLimit());
				GetVariable(i) -> SetLimits(XMin[i],XMax[i]);
			}
			CreateHistograms();
			// restore mins and maxes
			for (unsigned i=0; i<GetNVariables(); ++i)
				GetVariable(i) -> SetLimits(xmin[i],xmax[i]);
			
		}
	}
	

	fMCMCStatistics.assign(fMCMCNChains,BCEngineMCMC::MCMCStatistics(GetNParameters(),GetNObservables()));
	fMCMCStatistics_AllChains.Init(GetNParameters(),GetNObservables());
	fMCMCTree_Prob = -std::numeric_limits<double>::infinity();

	bool in_main_run = false;
	
	for (unsigned n=0; n<fMCMCTree->GetEntries(); ++n) {
		fMCMCTree -> GetEntry(n);
		if (!has_iteration)
			fMCMCTree_Iteration = n;

		fMCMCx[fMCMCTree_Chain]    = fMCMCTree_Parameters;
		fMCMCprob[fMCMCTree_Chain] = fMCMCTree_Prob;
		fMCMCLogLikelihood[fMCMCTree_Chain] = fMCMCTree_LogLikelihood;
		fMCMCLogPrior[fMCMCTree_Chain] = fMCMCTree_LogPrior;

		if (fMCMCTreeReuseObservables)
			fMCMCObservables[fMCMCTree_Chain] = fMCMCTree_Observables;
		else
			EvaluateObservables(fMCMCTree_Chain);

		if (!in_main_run and fMCMCPhase > 0) {
			for (unsigned c=0; c<fMCMCNChains; ++c)
				fMCMCStatistics[c].Reset(false,true);
			in_main_run = true;
		}

		fMCMCStatistics[fMCMCTree_Chain].Update(fMCMCTree_Prob,fMCMCTree_Parameters,fMCMCTree_Observables);
		
		if (fMCMCPhase <= 0)
			continue;

		MCMCCurrentPointInterface(fMCMCx[fMCMCTree_Chain], fMCMCTree_Chain, true);

		if (fMCMCTree_Chain==fMCMCNChains-1) {
			MCMCIterationInterface();
			if ( !fH1Marginalized.empty() or !fH2Marginalized.empty() )
				MCMCInChainFillHistograms();
		}

	}

	// combine statistics
	for (unsigned c=0; c<fMCMCNChains; ++c)
		fMCMCStatistics_AllChains += fMCMCStatistics[c];
}

// --------------------------------------------------------
double BCEngineMCMC::CalculateEvidence(double epsilon) {
	return 0;
}

// --------------------------------------------------------
bool BCEngineMCMC::UpdateCholeskyDecompositions() {
	if (fMultivariateProposalFunctionCovariance.size() != fMCMCNChains)
		return false;

	// double a = pow((double)fMultivariateProposalFunctionTuningSteps,-fMultivariateProposalFunctionTuningScheduleParameter);
	// Set covariance matricies
	unsigned I = 0;
	for (unsigned i=0; i<GetNParameters(); ++i) {
		if (GetParameter(i)->Fixed())
			continue;
		unsigned J = I;
		for (unsigned j=i; j<GetNParameters(); ++j) {
			if (GetParameter(j)->Fixed())
				continue;
			for (unsigned c=0; c<fMCMCNChains; ++c) {
				fMultivariateProposalFunctionCovariance[c][I][J] = /*(1-a) * fMultivariateProposalFunctionCovariance[c][I][J]	+ a **/ fMCMCStatistics[c].covariance[i][j];
				fMultivariateProposalFunctionCovariance[c][J][I] = fMultivariateProposalFunctionCovariance[c][I][J];
			}
			++J;
		}
		++I;
	}

	// create decomposer
	TDecompChol CholeskyDecomposer;
	// Update cholesky decompositions
	for (unsigned c=0; c<fMCMCNChains; ++c) {
		
		// try cholesky decomposition
		CholeskyDecomposer.SetMatrix(fMultivariateProposalFunctionCovariance[c]*fMCMCTrialFunctionScaleFactor[c][0]);
		if (CholeskyDecomposer.Decompose())
			fMultivariateProposalFunctionCholeskyDecomposition[c].Transpose(CholeskyDecomposer.GetU());
		
		else {
			// try with covariance + epsilon*1
			BCLog::OutDetail("BCEngineMCMC:UpdateCholeskyDecompositions : Cholesky decomposition failed! Adding epsilon*I and trying again.");
			TMatrixDSym U(fMultivariateProposalFunctionCovariance[c]*fMCMCTrialFunctionScaleFactor[c][0]);
			for (int i=0; i<U.GetNrows(); ++i)
				U[i][i] *= (1+fMultivariateProposalFunctionEpsilon);
			CholeskyDecomposer.SetMatrix(U);
			if (CholeskyDecomposer.Decompose())
				fMultivariateProposalFunctionCholeskyDecomposition[c].Transpose(CholeskyDecomposer.GetU());
			
			else {
				// diagonalize
				BCLog::OutDetail("BCEngineMCMC::UpdateCholeskyDecompositions : Cholesky decomposition failed! Setting off-diagonal elements of covariance to zero");
				for (int i=0; i<fMultivariateProposalFunctionCholeskyDecomposition[c].GetNrows(); ++i) {
					fMultivariateProposalFunctionCholeskyDecomposition[c][i][i] = sqrt(fMultivariateProposalFunctionCovariance[c][i][i]*fMCMCTrialFunctionScaleFactor[c][0]);
					for (int j=i+1; j<fMultivariateProposalFunctionCholeskyDecomposition[c].GetNcols(); ++j) {
						fMultivariateProposalFunctionCholeskyDecomposition[c][i][j] = 0;
						fMultivariateProposalFunctionCholeskyDecomposition[c][j][i] = 0;
					}
				}
			}
		}
	}
	return true;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCTrialFunction(unsigned ichain, std::vector<double> &x) {
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
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned chain, std::vector<double> &x) {
	x = fMCMCx[chain];

	// generate N-Free N(0,1) random values
	TVectorD y(GetNFreeParameters());
	for (int i=0; i<y.GetNrows(); ++i)
		y[i] = fMCMCThreadLocalStorage[chain].rng -> Gaus(0,1);

	// multiply by cholesky decomposition
	y *= fMultivariateProposalFunctionCholeskyDecomposition[chain];
	
	// add values into x
	int I = 0;
	for (unsigned i=0; i<GetNParameters() and I<y.GetNrows(); ++i)
		if (!GetParameter(i)->Fixed()) {
			x[i] += y[I];
			++I;
		}

	// return whether point is within limits, ignoring fixed parameters
	return GetParameters().IsWithinLimits(x,true);
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCGetProposalPointMetropolis(unsigned ichain, unsigned ipar, std::vector<double> &x) {
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
  x[ipar] += proposal * GetParameter(ipar)->GetRangeWidth();

  // check if the point is in the correct volume.
	return GetParameter(ipar) -> IsWithinLimits(x[ipar]);
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
				// increase efficiency
				fMCMCStatistics[chain].efficiency[parameter] += (1.-fMCMCStatistics[chain].efficiency[parameter])/(fMCMCStatistics[chain].n_samples_efficiency+1.);
				// copy the point
				fMCMCx[chain][parameter] = fMCMCThreadLocalStorage[chain].xLocal[parameter];
				// save the probability of the point
				fMCMCprob[chain] = p1;
				fMCMCLogLikelihood[chain] = fMCMCLogLikelihood_Provisional[chain];
				fMCMCLogPrior[chain] = fMCMCLogPrior_Provisional[chain];
				 
				// execute user code
				MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, true);
				return true;
			} else {
				// decrease efficiency
				fMCMCStatistics[chain].efficiency[parameter] *= 1.*fMCMCStatistics[chain].n_samples_efficiency/(fMCMCStatistics[chain].n_samples_efficiency+1.);
			}
		} else {						// new log(likelihood) was not a finite number
			BCLog::OutDebug(Form("Log(likelihood) evaluated to nan or inf in chain %i while varying parameter %s to %.3e",chain,GetParameter(parameter)->GetName().data(),fMCMCThreadLocalStorage[chain].xLocal[parameter]));
			// decrease efficiency
			fMCMCStatistics[chain].efficiency[parameter] *= 1.*fMCMCStatistics[chain].n_samples_efficiency/(fMCMCStatistics[chain].n_samples_efficiency+1.);
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
				// increase efficiency
				fMCMCStatistics[chain].efficiency[0] += (1-fMCMCStatistics[chain].efficiency[0])/(fMCMCStatistics[chain].n_samples_efficiency+1);
				
				// copy the point
				fMCMCx[chain] = fMCMCThreadLocalStorage[chain].xLocal;
				// save the probability of the point
				fMCMCprob[chain] = p1;
				fMCMCLogLikelihood[chain] = fMCMCLogLikelihood_Provisional[chain];
				fMCMCLogPrior[chain] = fMCMCLogPrior_Provisional[chain];

				// execute user code
				MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].xLocal, chain, true);
				return true;
			} else {
				// decrease efficiency
				fMCMCStatistics[chain].efficiency[0] *= 1.*fMCMCStatistics[chain].n_samples_efficiency/(fMCMCStatistics[chain].n_samples_efficiency+1);
			}
		} else { // new log(likelihood) was not a finite number
			BCLog::OutDebug("Log(likelihood) evaluated to nan or inf at");
			// decrease efficiency
			fMCMCStatistics[chain].efficiency[0] *= 1.*fMCMCStatistics[chain].n_samples_efficiency/(fMCMCStatistics[chain].n_samples_efficiency+1);
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

	if ( !fMCMCMultivariateProposalFunction ) { // run over pars one at a time

		for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar) {
			if ( GetParameter(ipar)->Fixed() )
				continue;

			//loop over chains
			unsigned chunk = 1; (void) chunk;
			unsigned ichain;    (void) ichain;
#pragma omp parallel for shared(chunk) private(ichain) schedule(static, chunk)
			for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain)
				return_value *= MCMCGetNewPointMetropolis(ichain,ipar);
				
			fMCMCCurrentChain = -1;
		}

	} else {											// run over all pars at once

		//loop over chains
		unsigned chunk = 1; (void) chunk;
		unsigned ichain;    (void) ichain;
#pragma omp parallel for shared(chunk) private(ichain) schedule(static, chunk)
		for (unsigned ichain = 0; ichain < fMCMCNChains; ++ichain)
			return_value *= MCMCGetNewPointMetropolis(ichain);

		fMCMCCurrentChain = -1;
	}

	// increase number of iterations used in each chain for calculating efficiencies
	for (unsigned c=0; c<fMCMCNChains; ++c)
		fMCMCStatistics[c].n_samples_efficiency += 1;

	++fMCMCCurrentIteration;
	return return_value;
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainFillHistograms() {
	// loop over chains
	for (unsigned c = 0; c < fMCMCNChains; ++c) {
		////////////////////////////////////////
		// fill each 1-dimensional histogram that exists
		for (unsigned j=0; j<GetNVariables() and j<fH1Marginalized.size(); ++j)
			if (dynamic_cast<TH1*>(fH1Marginalized[j])!=NULL)
				fH1Marginalized[j] -> Fill((j<GetNParameters()) ? fMCMCx[c][j] : fMCMCObservables[c][j-GetNParameters()]);
		
		////////////////////////////////////////
		// fill each 2-dimensional histogram that exists
		for (unsigned j=0; j<GetNVariables() and j<fH2Marginalized.size(); ++j)
			for (unsigned k=0; k<GetNVariables() and k<fH2Marginalized[j].size(); ++k)
				if (dynamic_cast<TH2*>(fH2Marginalized[j][k])!=NULL)
					fH2Marginalized[j][k] -> Fill((j<GetNParameters()) ? fMCMCx[c][j] : fMCMCObservables[c][j-GetNParameters()],
																				(k<GetNParameters()) ? fMCMCx[c][k] : fMCMCObservables[c][k-GetNParameters()]);
	}
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInChainWriteChains() {
	if (!fMCMCTree)
		return;
	// loop over all chains
	for (fMCMCTree_Chain=0; fMCMCTree_Chain < fMCMCNChains; ++fMCMCTree_Chain) {
		fMCMCTree_Prob          = fMCMCprob[fMCMCTree_Chain];
		fMCMCTree_LogLikelihood = fMCMCLogLikelihood[fMCMCTree_Chain];
		fMCMCTree_LogPrior      = fMCMCLogPrior[fMCMCTree_Chain];
		fMCMCTree_Parameters    = fMCMCx[fMCMCTree_Chain];
		fMCMCTree_Observables   = fMCMCObservables[fMCMCTree_Chain];
		fMCMCTree -> Fill();
	}
}

//---------------------------------------------------------
void BCEngineMCMC::MCMCCloseOutputFile() {
	if ( !fMCMCOutputFile or !fMCMCOutputFile->IsOpen() )
		return;
	if ( fMCMCOutputFile -> IsWritable() )
		fMCMCOutputFile -> Write(0,TObject::kWriteDelete);
	fMCMCOutputFile -> Close();
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCMetropolisPreRun() {
	// print on screen
	BCLog::OutSummary("Pre-run Metropolis MCMC...");
	
	// initialize Markov chain
	MCMCInitialize();

	if (fMCMCFlagWritePreRunToFile)
		InitializeMarkovChainTree();
	
	// perform run
	BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsPreRunMax));

	const int old_error_ignore_level = gErrorIgnoreLevel;

	if (fMCMCMultivariateProposalFunction) {
		// multivariate proposal function

		// suppress ROOT errors about Cholesky decomposition
		gErrorIgnoreLevel = kBreak;

		// initialize proposal function scale factors to 2.38^2 / number of dimensions
		fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,std::vector<double>(1,2.38*2.38/GetNFreeParameters()));

		// initialize covariance matrices as diag(var_0, var_1, var_2, ...)
		// initialize cholesky decomposition as diag(std_0, std_1, std_2, ...)
		// for free parameters only
		TMatrixDSym S0(GetNFreeParameters());
		TMatrixD    CD0(GetNFreeParameters(),GetNFreeParameters());
		unsigned I = 0;
		for (unsigned i=0; i<GetNParameters(); ++i)
			if (!(GetParameter(i)->Fixed())) {
				S0[I][I] = GetParameter(i) -> GetPriorVariance();
				CD0[I][I] = sqrt(fMCMCTrialFunctionScaleFactor[0][0]*S0[I][I]);
				++I;
			}
		fMultivariateProposalFunctionCovariance.assign(fMCMCNChains,S0);
		fMultivariateProposalFunctionCholeskyDecomposition.assign(fMCMCNChains,CD0);
		fMultivariateProposalFunctionTuningSteps = 0;
	} else if (fMCMCAutoSetTrialFunctionScaleFactors and GetParameters().ArePriorsSet(true)) {
		std::vector<double> temp;
		for (unsigned i=0; i<GetNParameters(); ++i)
			temp.push_back((GetParameter(i)->Fixed() or GetParameter(i)->GetRangeWidth()==0) ? 1 : GetParameter(i)->GetPriorStandardDeviation()/GetParameter(i)->GetRangeWidth());
		fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,temp);
	}

	//////////////////////////////////////////////////
	// Adjust scales until all parameters are in correct efficiency range in all chains
	bool allEfficient = false;
	bool inefficientScalesAdjustable = true;
	fMCMCCurrentIteration = 0;
	fMCMCPhase = BCEngineMCMC::MCMCPreRunEfficiencyCheck;

	// Cholesky Decomposer for multivariate proposal function
	TDecompChol CholeskyDecomposer;

	// do while not yet at max number of pre-run iterations
	// and not all efficiencies are within range
	// but all are scales are still tuneable
	// and make at least two iterations of the multivariate update
	while (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax and ((!allEfficient and inefficientScalesAdjustable) or (fMCMCMultivariateProposalFunction and fMultivariateProposalFunctionTuningSteps<=1))) {
		
		MCMCGetNewPointMetropolis();

		EvaluateObservables();

		for (unsigned c=0; c<fMCMCNChains; ++c)
			fMCMCStatistics[c].Update(fMCMCprob[c],fMCMCx[c],fMCMCObservables[c]);

		if (fMCMCFlagWritePreRunToFile)
			MCMCInChainWriteChains();

		if ( fMCMCStatistics.front().n_samples_efficiency != fMCMCNIterationsEfficiencyCheck)
			continue;

		// check efficiencies and update scale parameters
		allEfficient = true;
		inefficientScalesAdjustable = false;

		for (unsigned c=0; c<fMCMCNChains; ++c) {

			if (fMCMCMultivariateProposalFunction) {
				// multivariate proposal function

				if (fMCMCStatistics[c].efficiency[0] >= fMCMCEfficiencyMin and fMCMCStatistics[c].efficiency[0] <= fMCMCEfficiencyMax)
					// if chain efficiency is in range,
					continue;
					
				if (allEfficient) 		// print header if encountering first bad efficiency
					BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
				allEfficient = false;
				
				if (fMCMCStatistics[c].efficiency[0] < fMCMCEfficiencyMin) {
					// efficiency too low... decrease scale factor
					fMCMCTrialFunctionScaleFactor[c][0] /= fMultivariateProposalFunctionScaleMultiplier;
					
					if (fMCMCTrialFunctionScaleFactor[c][0] > fMCMCScaleFactorLowerLimit) {
						BCLog::OutDetail(Form("         chain %d is below %.0f %% (%4.1f %%). Scale decreased to %.4g", c, 100*fMCMCEfficiencyMin, 100*fMCMCStatistics[c].efficiency[0], fMCMCTrialFunctionScaleFactor[c][0]));
						// still room to tune
						inefficientScalesAdjustable = true;
					} else {
						// no more room to tune
						fMCMCTrialFunctionScaleFactor[c][0] = fMCMCScaleFactorLowerLimit;
						BCLog::OutDetail(Form("         chain %d is below %.0f %% (%4.1f %%). Scale now at lower limit (%.4g)", c, 100*fMCMCEfficiencyMin, 100*fMCMCStatistics[c].efficiency[0], fMCMCScaleFactorLowerLimit));
					}

				} else {
					// efficiency too high... increase scale factor
					fMCMCTrialFunctionScaleFactor[c][0] *= fMultivariateProposalFunctionScaleMultiplier;
					
					if (fMCMCTrialFunctionScaleFactor[c][0] < fMCMCScaleFactorUpperLimit) {
						// still room to tune
						BCLog::OutDetail(Form("         chain %d is above %.0f %% (%4.1f %%). Scale increased to %.4g", c, 100*fMCMCEfficiencyMax, 100*fMCMCStatistics[c].efficiency[0], fMCMCTrialFunctionScaleFactor[c][0]));
						inefficientScalesAdjustable = true;
					} else {
						// no more room to tune
						fMCMCTrialFunctionScaleFactor[c][0] = fMCMCScaleFactorUpperLimit;
						BCLog::OutDetail(Form("         chain %d is above %.0f %% (%4.1f %%). Scale now at upper limit (%.4g)", c, 100*fMCMCEfficiencyMax, 100*fMCMCStatistics[c].efficiency[0], fMCMCScaleFactorUpperLimit));
					}
				}

			} else {
				// factorized proposal function. loop over parameters
				for (unsigned p=0; p<GetNParameters(); ++p) {

					if (GetParameter(p)->Fixed())
						continue;

					if ( fMCMCStatistics[c].efficiency[p] >= fMCMCEfficiencyMin and fMCMCStatistics[c].efficiency[p] <= fMCMCEfficiencyMax)
						// if parameter efficiency is in range,
						continue;

					if (allEfficient)			// print header if first bad efficiency
						BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ",fMCMCCurrentIteration));
					allEfficient = false;
						
					if (fMCMCStatistics[c].efficiency[p] < fMCMCEfficiencyMin) {
						// efficiency too low... decrease scale factor
						fMCMCTrialFunctionScaleFactor[c][p] /= (fMCMCStatistics[c].efficiency[p] < fMCMCEfficiencyMin/2) ? 4 : 2;
						
						if (fMCMCTrialFunctionScaleFactor[c][p] > fMCMCScaleFactorLowerLimit ) {
							// still room to tune
							BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p)->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCStatistics[c].efficiency[p], c, fMCMCTrialFunctionScaleFactor[c][p]));
							inefficientScalesAdjustable = true;
						}	else {
							// no more room to tune
							fMCMCTrialFunctionScaleFactor[c][p] = fMCMCScaleFactorLowerLimit;
							BCLog::OutDetail(Form("         %-*s is below %.0f %% (%4.1f %%) in chain %i. Scale now at lower limit (%.4g %%)",	fParameters.MaxNameLength(), GetParameter(p)->GetName().data(), 100*fMCMCEfficiencyMin, 100*fMCMCStatistics[c].efficiency[p], c, fMCMCScaleFactorLowerLimit));
						}
						
					} else {
						// if efficiency too high ... increase scale factor
						fMCMCTrialFunctionScaleFactor[c][p] *= 2;
						
						if ( fMCMCTrialFunctionScaleFactor[c][p] < fMCMCScaleFactorUpperLimit ) {
							// still room to tune
							BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale increased to %.4g", fParameters.MaxNameLength(), GetParameter(p)->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCStatistics[c].efficiency[p], c, fMCMCTrialFunctionScaleFactor[c][p]));
							inefficientScalesAdjustable = true;
						} else {
							// no more room to tune
							fMCMCTrialFunctionScaleFactor[c][p] = fMCMCScaleFactorUpperLimit;
							BCLog::OutDetail(Form("         %-*s is above %.0f %% (%4.1f %%) in chain %i. Scale now at upper limit (%.4g)", fParameters.MaxNameLength(), GetParameter(p)->GetName().data(), 100*fMCMCEfficiencyMax, 100*fMCMCStatistics[c].efficiency[p], c, fMCMCScaleFactorUpperLimit));														 
						}
					}
				}
			}
		}

		if (fMCMCFlagWritePreRunToFile and fMCMCTree)
			fMCMCTree -> AutoSave("SaveSelf");
		
		if (fMCMCMultivariateProposalFunction)
			++fMultivariateProposalFunctionTuningSteps;

		if ((allEfficient or !inefficientScalesAdjustable) and (!fMCMCMultivariateProposalFunction or fMultivariateProposalFunctionTuningSteps>1))
			// now efficient, or no longer tuneable
			// BUT, not if this is the first iteration of the multivariate scheme
			continue;
		
		if (fMCMCMultivariateProposalFunction)
			UpdateCholeskyDecompositions();
		
		// reset statistics
		for (unsigned c=0; c<fMCMCStatistics.size(); ++c)
			fMCMCStatistics[c].Reset(false,true); // preserve mode information
	}

	// restore ROOT error ignore level
	gErrorIgnoreLevel = old_error_ignore_level;

	if (allEfficient)
		BCLog::OutDetail(Form("     * Efficiency status: Efficiencies within predefined range after %i iterations.",fMCMCCurrentIteration));
	else if (!inefficientScalesAdjustable)
		BCLog::OutWarning(Form("     * Efficiency status: Some efficiencies outside predefined range, but scales are at limits after %i iterations.",fMCMCCurrentIteration));
	else
		BCLog::OutDetail(Form("     * Efficiency status: Some efficiencies outside predefined range, but maximum number of iterations (%i) reached.",fMCMCNIterationsPreRunMax));

	if (fMCMCNChains > 1) {
		//////////////////////////////////////////////////
		// Run until all chains have converged

		unsigned nIterationsCheckConvergence = fMCMCNIterationsConvergenceCheck;

		if ( fMCMCNIterationsClearConvergenceStats > 0 and nIterationsCheckConvergence > fMCMCNIterationsClearConvergenceStats )
			nIterationsCheckConvergence = fMCMCNIterationsClearConvergenceStats;

		fMCMCNIterationsConvergenceGlobal = -1;

		// statistics about inter_chain parameters
		BCEngineMCMC::MCMCStatistics inter_chain(fMCMCStatistics.front().variance.size(),fMCMCStatistics.front().mean.size());
		
		// flag for whether to make one R value check, if already above max number of pre-run iterations
		bool make_one_check = false; 

		if (fMCMCCurrentIteration >= (int)fMCMCNIterationsPreRunMax) {
			make_one_check = (fMCMCStatistics.front().n_samples > 0);
			if (!make_one_check) {
				BCLog::OutWarning(" Convergence never checked !");
				BCLog::OutWarning("   Increase maximum number of iterations in the pre-run via MCMCSetNIterationsPreRunMax()");
			} else if (fMCMCStatistics.front().n_samples < nIterationsCheckConvergence) {
				BCLog::OutWarning(TString::Format(" Convergence will be checked only once and with only %d iterations!",fMCMCStatistics.front().n_samples));
				BCLog::OutWarning("   Increase maximum number of iterations in the pre-run via MCMCSetNIterationsPreRunMax()");
			}
		}

		fMCMCPhase = BCEngineMCMC::MCMCPreRunConvergenceCheck;
		while ( make_one_check or (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax and fMCMCNIterationsConvergenceGlobal < 0) ) {

			// clear statistics
			if ( !make_one_check and fMCMCNIterationsClearConvergenceStats > 0 and fMCMCStatistics.front().n_samples >= fMCMCNIterationsClearConvergenceStats )
				for (unsigned c=0; c<fMCMCNChains; ++c)
					fMCMCStatistics[c].Reset(false,false); // keep mode and efficiency stats

			MCMCGetNewPointMetropolis();
			EvaluateObservables();

			// update means, variances, etc.
			for (unsigned c=0; c<fMCMCNChains; ++c)
				fMCMCStatistics[c].Update(fMCMCprob[c],fMCMCx[c],fMCMCObservables[c]);

			if (fMCMCFlagWritePreRunToFile)
				MCMCInChainWriteChains();

			if ( !make_one_check and fMCMCStatistics.front().n_samples % nIterationsCheckConvergence != 0 and fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax )
				continue;

			make_one_check = false;

			// calculate R values according to Brooks & Gelman,
			// "General Methods for Monitoring Convergence of Iterative Simulations, 1998

			fMCMCStatistics_AllChains.Reset();
			inter_chain.Reset();
			for (unsigned c=0; c<fMCMCNChains; ++c) {
				fMCMCStatistics_AllChains += fMCMCStatistics[c];
				inter_chain.Update(fMCMCStatistics[c].probability_variance,fMCMCStatistics[c].variance,fMCMCStatistics[c].mean);
			}

			// parameter R values
			for (unsigned i=0; i<GetNParameters(); ++i) {
				if (GetParameter(i)->Fixed())
					continue;
				if (inter_chain.mean[i] == 0) { // mean of variance = 0
					if (fMCMCStatistics_AllChains.variance[i] == 0) { // variance of all samples = 0
						BCLog::OutDebug("BCEngineMCMC::MCMCMetropolisPreRun : all samples in all chains identical!");
						fMCMCRValueParameters[i] = 1;
					} else {
						BCLog::OutDebug("BCEngineMCMC::MCMCMetropolisPreRun : all chain variances zero!");
						fMCMCRValueParameters[i] = std::numeric_limits<double>::infinity();
						}
				} else {
					fMCMCRValueParameters[i] = sqrt( fMCMCStatistics_AllChains.variance[i] / inter_chain.mean[i] ); // variance(all samples) / mean(chain variances)
				}
			}

			if (inter_chain.probability_mean == 0) { // mean of probability variances = 0
				if (fMCMCStatistics_AllChains.probability_variance == 0) { // variance of all probabilities = 0
					BCLog::OutDebug("BCEngineMCMC::MCMCMetropolisPreRun : all likelihoods in all samples identical!");
					fMCMCRValue = 1;
				} else {
					BCLog::OutDebug("BCEngineMCMC::MCMCMetropolisPreRun : variance of all likelihoods in all chains zero!");
					fMCMCRValue = std::numeric_limits<double>::infinity();
				}
			} else {
				fMCMCRValue = sqrt( fMCMCStatistics_AllChains.probability_variance / inter_chain.probability_mean ); // variance(all samples) / mean(chain variances)
			}

			if (fCorrectRValueForSamplingVariability) {
				double m = static_cast<double>(fMCMCNChains);
				double n = static_cast<double>(fMCMCStatistics.front().n_samples);
				double N = (n-1)/n;
				double M = (m+1)/m;
				double K1 = N*N/m;
				double K2 = M * 2/(m-1);
				double K3 = 2*N*M/m;

				// correct parameters (i>=0) and posterior (i==-1) for initial sampling variability
				for (int i=-1; i<(int)GetNParameters(); ++i) {
					if (i>=0 and GetParameter(i)->Fixed())
						continue;
					double W = (i>=0) ? fMCMCStatistics_AllChains.variance[i]                 : fMCMCStatistics_AllChains.probability_variance; // variance of all samples
					double S = (i>=0) ? inter_chain.mean[i]                                   : inter_chain.probability_mean; // mean of chain variances
					double Z = (i>=0) ? inter_chain.variance[i]                               : inter_chain.probability_variance;	// variance of chain variances
					double X = (i>=0) ? inter_chain.mean[inter_chain.efficiency.size()+i]     : 0; // mean of chain means
					double B = (i>=0) ? inter_chain.variance[inter_chain.efficiency.size()+i] : 0; // variance of chain means
					// calculate X and B for the probability
					if (i<0) {
						for (unsigned c=0; c<fMCMCNChains; ++c)
							X += fMCMCStatistics[c].probability_mean;
						X /= fMCMCNChains;
						for (unsigned c=0; c<fMCMCNChains; ++c)
							B+= (fMCMCStatistics[c].probability_mean-X)*(fMCMCStatistics[c].probability_mean-X);
						B /= fMCMCNChains - 1.;
					}

					double V = N*W + M*B;

					double X2 = 0;
					for (unsigned c=0; c<fMCMCNChains; ++c)
						X2 += (i>=0) ? fMCMCStatistics[c].mean[i] * fMCMCStatistics[c].mean[i] : fMCMCStatistics[c].probability_mean * fMCMCStatistics[c].probability_mean;
					X2 /= fMCMCNChains;
					double cov_s_x = 0;
					double cov_s_x2 = 0;
					for (unsigned c=0; c<fMCMCNChains; ++c) {
						double delta_s = (i>=0) ? fMCMCStatistics[c].variance[i]-S : fMCMCStatistics[c].probability_variance-S;
						cov_s_x  += (i>=0) ? delta_s*(fMCMCStatistics[c].mean[i]-X)                             : delta_s*(fMCMCStatistics[c].probability_mean-X);
						cov_s_x2 += (i>=0) ? delta_s*(fMCMCStatistics[c].mean[i]*fMCMCStatistics[c].mean[i]-X2) : delta_s*(fMCMCStatistics[c].probability_mean*fMCMCStatistics[c].probability_mean-X2);
					}
					cov_s_x  /= fMCMCNChains - 1.;
					cov_s_x2 /= fMCMCNChains - 1.;
					double varV = K1*Z + K2*B*B + K3*(cov_s_x2 - 2.*X*cov_s_x);

					double df = 2*V*V/varV;
					if (i>=0)
						fMCMCRValueParameters[i] *= sqrt((df+3.)/(df+1.));
					else
						fMCMCRValue *= sqrt((df+3.)/(df+1.));
				}
			}
		
			// check R values
			fMCMCNIterationsConvergenceGlobal = fMCMCCurrentIteration;
			for (unsigned i=0; i<GetNParameters() and fMCMCNIterationsConvergenceGlobal>0; ++i)
				if ( ! GetParameter(i)->Fixed() and (fMCMCRValueParameters[i]-1) > fMCMCRValueParametersCriterion )
					fMCMCNIterationsConvergenceGlobal = -1;
			if (fMCMCNIterationsConvergenceGlobal > 0 and (fMCMCRValue-1) > fMCMCRValueCriterion)
				fMCMCNIterationsConvergenceGlobal = -1;

			if ( fMCMCNIterationsConvergenceGlobal > 0 )
				continue;

			BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));
			BCLog::OutDetail(Form("       - %-*s : R-Value",fParameters.MaxNameLength(),"Parameter"));
				
			for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar) {

				if ( GetParameter(ipar)->Fixed() )
					continue;

				if( fMCMCRValueParameters[ipar]-1 < fMCMCRValueParametersCriterion )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else if ( fMCMCRValueParameters.at(ipar) != std::numeric_limits<double>::infinity() )
					BCLog::OutDetail(TString::Format("         %-*s :  %.06f <--",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data(), fMCMCRValueParameters.at(ipar)));
				else
					BCLog::OutDetail(TString::Format("         %-*s :  INFINITY <--",fParameters.MaxNameLength(),GetParameter(ipar)->GetName().data()));
			}

			if( fMCMCRValue-1 < fMCMCRValueCriterion )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f", fMCMCRValue));
			else if ( fMCMCRValue != std::numeric_limits<double>::max() )
				BCLog::OutDetail(Form("       - Log-Likelihood :  %.06f <--", fMCMCRValue));
			else
				BCLog::OutDetail("       - Log-Likelihood :  INFINITY <--");

			if (fMCMCFlagWritePreRunToFile and fMCMCTree)
				fMCMCTree -> AutoSave("SaveSelf");
			
			if (fMCMCMultivariateProposalFunction)
				UpdateCholeskyDecompositions();
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

		if (fMCMCFlagWritePreRunToFile and fMCMCTree)
			fMCMCTree -> AutoSave("SaveSelf");
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
			EvaluateObservables();
			++n;
			
			for (unsigned c=0; c<fMCMCNChains; ++c)
				fMCMCStatistics[c].Update(fMCMCprob[c],fMCMCx[c],fMCMCObservables[c]);

			if (fMCMCFlagWritePreRunToFile)
				MCMCInChainWriteChains();

			if ( n % nwrite == 0) {
				BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, 100.*n/N));
				if (fMCMCFlagWritePreRunToFile and fMCMCTree)
					fMCMCTree -> AutoSave("SaveSelf");
			}
		}
		BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations in pre-run.", fMCMCCurrentIteration));
	}

	// combine statistics:
	fMCMCStatistics_AllChains.Reset();
	for (unsigned c=0; c<fMCMCNChains; ++c)
		fMCMCStatistics_AllChains += fMCMCStatistics[c];

	if (fMCMCMultivariateProposalFunction) {
		BCLog::OutDetail(Form(" --> Scale factors and efficiencies (measured in %d iterations):",fMCMCStatistics.front().n_samples_efficiency));
		BCLog::OutDetail("       - Chain : Scale factor    Efficiency");
		for (unsigned c=0; c<fMCMCNChains; ++c)
			BCLog::OutDetail(Form("         %-3d :         % 6.4g         %4.1f %%",c, fMCMCTrialFunctionScaleFactor[c][0], 100.*fMCMCStatistics[c].efficiency[0]));
	} else {
		BCLog::OutDetail(Form(" --> Average scale factors and efficiencies (measured in %d iterations):",fMCMCStatistics.front().n_samples_efficiency));
		BCLog::OutDetail(Form("       - %-*s : Scale factor    Efficiency",fParameters.MaxNameLength(),"Parameter"));
		// print scale factors and efficiencies
		std::vector<double> scalefactors (GetNParameters(),0);
		for (unsigned i = 0; i < GetNParameters(); ++i) {
			if (GetParameter(i)->Fixed())
				continue;
			for (unsigned j = 0; j < fMCMCNChains; ++j)
				scalefactors[i] += fMCMCTrialFunctionScaleFactor[j][i] / fMCMCNChains;
			BCLog::OutDetail(Form("         %-*s :          % 6.4g %%        %4.1f %%",fParameters.MaxNameLength(),GetParameter(i)->GetName().data(), 100.*scalefactors[i], 100.*fMCMCStatistics_AllChains.efficiency[i]));
		}
	}

	// reset current iteration
	fMCMCCurrentIteration = -1;

	// reset current chain
	fMCMCCurrentChain = -1;
	
	if (fMCMCFlagWritePreRunToFile) {
		// UpdateParameterTree();
		if (fMCMCTree)
			fMCMCTree -> AutoSave("SaveSelf");
	}

	// rescale histograms
	if (fRescaleHistogramRangesAfterPreRun)
		CreateHistograms(true);

	// no error
	return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCMetropolis() {
  // check the number of free parameters
  if (GetNFreeParameters() <= 0) {
    BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Number of free parameters <= 0. Do not run Metropolis.");
    return false;
  }

	// check if prerun should be performed
	if (fMCMCFlagPreRun) {
		if (!MCMCMetropolisPreRun())
			return false;
		if (!fMCMCFlagWritePreRunToFile and fMCMCFlagWriteChainToFile)
			InitializeMarkovChainTree();
	}
	else {
		BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Not running prerun. This can cause trouble if the data have changed.");
		if (fMCMCFlagWriteChainToFile)
			InitializeMarkovChainTree();
	}

	// make sure enough statistics containers exist
	fMCMCStatistics.resize(fMCMCNChains,BCEngineMCMC::MCMCStatistics(GetNParameters(),GetNObservables()));

	// reset statistics
	for (unsigned c=0; c<fMCMCNChains; ++c)
		fMCMCStatistics[c].Reset(false,false); // keep mode and efficiency information

	// print to screen
	BCLog::OutSummary( "Run Metropolis MCMC ...");

	// set phase and cycle number
	fMCMCPhase = BCEngineMCMC::MCMCMainRun;

	BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

	unsigned nwrite = UpdateFrequency(fMCMCNIterationsRun);

	// start the run
	fMCMCCurrentIteration = 0;
	while ( fMCMCCurrentIteration < (int)fMCMCNIterationsRun ) {

		MCMCGetNewPointMetropolis();
		EvaluateObservables();

		if ( fMCMCCurrentIteration % nwrite == 0 ) {
			BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration)/(double)fMCMCNIterationsRun*100.));
			if (fMCMCFlagWriteChainToFile and fMCMCTree)
				fMCMCTree -> AutoSave("SaveSelf");
		}
		
		if (fMCMCCurrentIteration % fMCMCNLag != 0) // apply lag
			continue;

		MCMCIterationInterface();		// user action (overloadable)

		for (unsigned c=0; c<fMCMCNChains; ++c)
			fMCMCStatistics[c].Update(fMCMCprob[c],fMCMCx[c],fMCMCObservables[c]);

		// fill histograms
		if ( !fH1Marginalized.empty() or !fH2Marginalized.empty() )
			MCMCInChainFillHistograms();

		// write chain to file
		if ( fMCMCFlagWriteChainToFile )
			MCMCInChainWriteChains();

	} // end run

	BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

	// set all-chain stats to first-chain stats
	fMCMCStatistics_AllChains.Reset();
	// loop over remaining chains
	for (unsigned c=0; c<fMCMCStatistics.size(); ++c)
		fMCMCStatistics_AllChains += fMCMCStatistics[c];

	// print efficiencies
	if (fMCMCMultivariateProposalFunction) {
		BCLog::OutDetail(Form(" --> Efficiencies (measured in %d iterations):",fMCMCStatistics.front().n_samples_efficiency));
		BCLog::OutDetail("       - Chain : Efficiency");
		for (unsigned c=0; c<fMCMCNChains; ++c)
			BCLog::OutDetail(Form("         %-3d :       %4.1f %%",c, 100.*fMCMCStatistics[c].efficiency[0]));
	} else {
		BCLog::OutDetail(Form(" --> Average efficiencies (measured in %d iterations):",fMCMCStatistics_AllChains.n_samples_efficiency/fMCMCNChains));
		BCLog::OutDetail(Form("       - %-*s : Efficiency",fParameters.MaxNameLength(),"Parameter"));
		for (unsigned i = 0; i < GetNParameters(); ++i) {
			if (GetParameter(i)->Fixed())
				continue;
			BCLog::OutDetail(Form("         %-*s :     %4.1f %%",fParameters.MaxNameLength(),GetParameter(i)->GetName().data(), 100.*fMCMCStatistics_AllChains.efficiency[i]));
		}
	}

	if (fMCMCFlagWriteChainToFile)
		UpdateParameterTree();

	BCLog::OutDetail(" --> Global mode from MCMC:");
	BCLog::OutDebug(Form(" --> Posterior value: %g", fMCMCStatistics_AllChains.probability_mode));
	PrintParameters(fMCMCStatistics_AllChains.mode,BCLog::OutDetail);

	// reset counter
	fMCMCCurrentIteration = -1;

	// reset current chain
	fMCMCCurrentChain = -1;

	// set flags
	fMCMCFlagRun = true;

	return true;
}

// --------------------------------------------------------
void BCEngineMCMC::EvaluateObservables() {
	for (unsigned c=0; c<fMCMCNChains; ++c)
		EvaluateObservables(c);
}

// --------------------------------------------------------
void BCEngineMCMC::EvaluateObservables(unsigned chain) {
	if (chain>fMCMCNChains)
		return;
	fMCMCCurrentChain = chain;
	CalculateObservables(fMCMCx[fMCMCCurrentChain]);
	for (unsigned j = 0; j < GetNObservables(); ++j)
		fMCMCObservables[fMCMCCurrentChain][j] = GetObservable(j) -> Value();
}

// --------------------------------------------------------
void BCEngineMCMC::ResetResults() {
	// reset variables
	fMCMCCurrentIteration = -1;
	fMCMCCurrentChain = -1;
	fMCMCStatistics.clear();
	fMCMCStatistics_AllChains.Clear();
	fMCMCTrialFunctionScaleFactor.clear();
	fMCMCx.clear();
	fMCMCObservables.clear();
	fMCMCprob.clear();
	fMCMCLogLikelihood.clear();
	fMCMCLogLikelihood_Provisional.clear();
	fMCMCLogPrior.clear();
	fMCMCLogPrior_Provisional.clear();
	fMCMCNIterationsConvergenceGlobal = -1;
	fMCMCRValueParameters.clear();
	fMCMCRValue = std::numeric_limits<double>::max();

	for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
		delete fH1Marginalized[i];
	
	for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
		for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
			delete fH2Marginalized[i][j];
	
	// clear plots
	fH1Marginalized.clear();
	fH2Marginalized.clear();

	// reset flags
	fMCMCFlagRun = false;

	fLocalModes.clear();
}

// --------------------------------------------------------
bool BCEngineMCMC::MCMCInitialize() {
	// resource allocation must be done only by one thread

	fMCMCNIterationsConvergenceGlobal = -1;

	// free memory for vectors
	fMCMCNIterations.assign(fMCMCNChains, 0);
	fMCMCStatistics.assign(fMCMCNChains,BCEngineMCMC::MCMCStatistics(GetNParameters(),GetNObservables()));
	fMCMCStatistics_AllChains.Init(GetNParameters(),GetNObservables());
	fMCMCprob.assign(fMCMCNChains, -std::numeric_limits<double>::infinity());
	fMCMCLogLikelihood.assign(fMCMCNChains, -std::numeric_limits<double>::infinity());
	fMCMCLogLikelihood_Provisional.assign(fMCMCNChains, -std::numeric_limits<double>::infinity());
	fMCMCLogPrior.assign(fMCMCNChains, -std::numeric_limits<double>::infinity());
	fMCMCLogPrior_Provisional.assign(fMCMCNChains, -std::numeric_limits<double>::infinity());

	fMCMCRValueParameters.assign(GetNParameters(), std::numeric_limits<double>::infinity());
	fMCMCRValue = std::numeric_limits<double>::max();
	
	// clear info about local modes
	fLocalModes.clear();

	SyncThreadStorage();

   if (fMCMCTrialFunctionScaleFactorStart.size() == 0 or fMCMCTrialFunctionScaleFactorStart.size()!=GetNParameters())
		 fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,std::vector<double>(GetNParameters(), 1.0));
   else
		 fMCMCTrialFunctionScaleFactor.assign(fMCMCNChains,fMCMCTrialFunctionScaleFactorStart);

	 // initialize markov chain positions
	 switch (fMCMCFlagInitialPosition) {

		 // keep previous values
	 case kMCMCInitContinue : {
		 // check position vector size
		 if (fMCMCx.size() != fMCMCNChains) {
			 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Number of chains has been changed; cannot continue previous chains.");
			 return false;
		 }
		 // else do nothing --- continue chains
		 break;
	 }

		 // use range centers
	 case kMCMCInitCenter : {
		 fMCMCx.assign(fMCMCNChains,GetParameters().GetRangeCenters(true));
		 break;
	 }

		 // uniformly distribute all coordinates in provided ranges
	 case kMCMCInitRandomUniform : {
		 fMCMCx.clear();
		 for (unsigned c=0; c<fMCMCNChains; ++c)
			 fMCMCx.push_back(GetParameters().GetUniformRandomValues(fMCMCThreadLocalStorage[c].rng,true));
		 break;
	 }

		 // use user-defined starting points
	 case kMCMCInitUserDefined : {
		 // check initial position vector size
		 if (fMCMCInitialPosition.size() < fMCMCNChains) {
			 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Too few initial positions provided.");
			 return false;
		 }
		 // copy positions
		 fMCMCx.assign(fMCMCInitialPosition.begin(),fMCMCInitialPosition.begin()+fMCMCNChains);
		 // set fixed values then check whether initial positions are within bounds
		 // check whether initial positions are within bounds
		 // also checks that initial position vectors are correct size
		 for (unsigned c=0; c<fMCMCNChains; ++c) {
			 GetParameters().ApplyFixedValues(fMCMCx[c]);
			 if (!GetParameters().IsWithinLimits(fMCMCx[c],true)) {
				 BCLog::OutError("BCEngineMCMC::MCMCInitialize : User-defined initial point is out of bounds.");
				 fMCMCx.clear();
				 return false;
			 }
		 }
		 break;
	 }

		 // randomly distribute according to factorized priors
	 case kMCMCInitRandomPrior : {
		 if (!GetParameters().ArePriorsSet(true)) {
			 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Not all unfixed parameters have priors set.");
			 return false; 
		 }
		 fMCMCx.clear();
		 for (unsigned c=0; c<fMCMCNChains; ++c) {
			 fMCMCx.push_back(GetParameters().GetRandomValuesAccordingToPriors(fMCMCThreadLocalStorage[c].rng,true));
			 // check new point
			 if (!GetParameters().IsWithinLimits(fMCMCx[c],true)) {
				 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Could not generate random point within limits.");
				 fMCMCx.clear();
				 return false;
			 }
		 }
		 break;
	 }

		 // use mean and variance of factorized priors
		 // to distribute values according to normal distribution,
		 // allowing for over-distribution of parameters.
	 case kMCMCInitRandomGaussPrior : {
		 if (!GetParameters().ArePriorsSet(true)) {
			 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Not all unfixed parameters have priors set.");
			 return false; 
		 }
		 fMCMCx.clear();
		 for (unsigned c=0; c<fMCMCNChains; ++c) {
			 // add new point
			 fMCMCx.push_back(GetParameters().GetRandomValuesAccordingToGaussiansOfPriors(fMCMCThreadLocalStorage[c].rng,true,fMCMCInitialPositionExpansionFactor));
			 // check new point
			 if (!GetParameters().IsWithinLimits(fMCMCx[c],true)) {
				 BCLog::OutError("BCEngineMCMC::MCMCInitialize : Could not generate random point within limits.");
				 fMCMCx.clear();
				 return false;
			 }
		 }
		 break;
	 }

	 default:
		 fMCMCx.clear();
		 BCLog::OutError("BCEngineMCMC::MCMCInitialize : No MCMC position initialization scheme specified.");
		 return false;
	 }
	 
	 if (fMCMCx.empty())
		 return false;

	 // initialize user-defined observables
	 fMCMCObservables.assign(fMCMCNChains,std::vector<double>(GetNObservables(),0));

	 CreateHistograms(false);

	 // set that a main run has not been made
	 fMCMCFlagRun = false;

	 return true;
}

// ------------------------------------------------------------
void BCEngineMCMC::CreateHistograms(bool rescale_ranges) {
	// clear existing histograms
	for (unsigned i=0; i<fH1Marginalized.size(); ++i)
		delete fH1Marginalized[i];
	fH1Marginalized.assign(GetNVariables(),NULL);

	for (unsigned i=0; i<fH2Marginalized.size(); ++i)
		for (unsigned j=0; j<fH2Marginalized[i].size(); ++j)
			delete fH2Marginalized[i][j];
	fH2Marginalized.assign(GetNVariables(),std::vector<TH2*>(GetNVariables(),NULL));
	
	// store old bounds, rescale if rescaling:
	std::vector<std::pair<double,double> > original_bounds;
	original_bounds.reserve(GetNVariables());
	for (unsigned i=0; i<GetNVariables(); ++i) {
		original_bounds.push_back(std::make_pair(GetVariable(i)->GetLowerLimit(),GetVariable(i)->GetUpperLimit()));
		if (rescale_ranges) {
			if (i < fMCMCStatistics_AllChains.minimum.size() and std::isfinite(fMCMCStatistics_AllChains.minimum[i]))
				GetVariable(i) -> SetLowerLimit(fMCMCStatistics_AllChains.minimum[i]);
			if (i < fMCMCStatistics_AllChains.maximum.size() and std::isfinite(fMCMCStatistics_AllChains.maximum[i]))
				GetVariable(i) -> SetUpperLimit(fMCMCStatistics_AllChains.maximum[i]);
			if (fHistogramRescalePadding>0) {
				// calculate enlargement factor
				double range_width_rescale = GetVariable(i) -> GetRangeWidth() * fHistogramRescalePadding;
				// push bounds out, but not beyond original parameter bounds
				GetVariable(i) -> SetLowerLimit(std::max<double>(GetVariable(i)->GetLowerLimit()-range_width_rescale,original_bounds.back().first));
				GetVariable(i) -> SetUpperLimit(std::min<double>(GetVariable(i)->GetUpperLimit()+range_width_rescale,original_bounds.back().second));
			}
		}
	}

	// define 1-dimensional histograms for marginalization
	int filling = 0;

	for (unsigned i=0; i<GetNVariables(); ++i)
		if (GetVariable(i)->FillH1()) {
			if (i<GetNParameters()) {	// parameter
				if (!GetParameter(i)->Fixed()) {
					fH1Marginalized[i] = GetVariable(i) -> CreateH1(TString::Format("h1_%s_parameter_%i", GetSafeName().data() ,i));
					++filling;
				}
			} else {									// user-defined observable
				fH1Marginalized[i] = GetVariable(i) -> CreateH1(TString::Format("h1_%s_observable_%i", GetSafeName().data() ,i-GetNParameters()));
				++filling;
			}
		}
	
	if (filling==0)	// if filling no 1D histograms, clear vector
		fH1Marginalized.clear();

	// define 2D histograms for marginalization
	filling = 0;

	// parameter i as abscissa:
	for(unsigned i = 0; i < GetNParameters(); ++i) {
		if (GetParameter(i)->Fixed())
			continue;
		
		// parameter j as ordinate
		for (unsigned j = 0; j < GetNParameters(); ++j)
			if (!(GetParameter(j)->Fixed()) and GetParameters().FillH2(i,j)) {
				fH2Marginalized[i][j] = GetParameter(i) -> CreateH2(Form("h2_%s_parameters_%i_vs_%i", GetSafeName().data(), i, j), GetParameter(j));
				++filling;
			}
		// user-defined observable j as ordinate
		for (unsigned j = 0; j < GetNObservables(); ++j)
			if (GetParameters().FillH2Partner(i,j)) {
				fH2Marginalized[i][j+GetNParameters()] = GetParameter(i) -> CreateH2(Form("h2_%s_par_%i_vs_obs_%i", GetSafeName().data(), i, j), GetObservable(j));
				++filling;
			}
	}

	// user-defined observable i as abscissa
	for(unsigned i = 0; i < GetNObservables(); ++i) {

		// user-defined observable j as ordinate
		for (unsigned j = 0; j < GetNObservables(); ++j)
			if (GetObservables().FillH2(i,j)) {
				fH2Marginalized[i+GetNParameters()][j+GetNParameters()] = GetObservable(i) -> CreateH2(Form("h2_%s_observables_%i_vs_%i", GetSafeName().data(), i, j), GetObservable(j));
				++filling;
			}

		// parameter j as ordinate
		for (unsigned j = 0; j < GetNParameters(); ++j)
			if (!(GetParameter(j)->Fixed()) and GetObservables().FillH2Partner(i,j)) {
				fH2Marginalized[i+GetNParameters()][j] = GetObservable(i) -> CreateH2(Form("h2_%s_obs_%i_vs_par_%i", GetSafeName().data(), i, j), GetParameter(j));
				++filling;
			}
	}
	
	// for (unsigned i=0; i<fH2Marginalized.size(); ++i)
	// 	for (unsigned j=0; j<fH2Marginalized[i].size(); ++j)
	// 		if (dynamic_cast<TH2*>(fH2Marginalized[i][j])!=NULL)
	// 			fH2Marginalized[i][j] -> SetDirectory(0);

	if (filling==0)	// if filling no 2D histograms, clear vector
		fH2Marginalized.clear();
		
	// restore bounds
	for (unsigned i=0; i<GetNVariables(); ++i)
		GetVariable(i) -> SetLimits(original_bounds[i].first,original_bounds[i].second);
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintSummary() {
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
	if ( GetGlobalMode().empty() )
		return;
	 
	BCLog::OutSummary(Form("Log of the maximum posterior: %f", GetLogMaximum()));
	BCLog::OutSummary("Best fit results:");

	for (unsigned i = 0; i < GetNVariables(); i++) {
		if (i < GetNParameters() and GetParameter(i)->Fixed() )
			BCLog::OutSummary(Form(" %s = %.*f (fixed)",  GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(),GetParameter(i)->GetFixedValue()));
		else
			BCLog::OutSummary(Form(" %s = %.*f (global)", GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(),GetGlobalMode()[i]));
		 
		if ( (GetLocalModes().size()==GetNParameters() or GetLocalModes().size()==GetNVariables()) and i<GetLocalModes().size() )
			BCLog::OutSummary(Form(" %s = %.*f (marginalized)", GetVariable(i)->GetName().data(), GetVariable(i)->GetPrecision(), GetLocalModes()[i]));
	}
}

// ---------------------------------------------------------
bool BCEngineMCMC::PrintResults(const char * file) const {
	// open file
	std::ofstream ofi(file);

	// check if file is open
	if (!ofi.is_open()) {
		BCLog::OutError(Form("Couldn't open file %s.", file));
		return false;
	}

	PrintSummaryToStream(ofi);
	PrintBestFitToStream(ofi);
	PrintMarginalizationToStream(ofi);
	
	// close file
	ofi.close();
	return true;
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintSummaryToStream(std::ofstream & ofi) const {
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
void BCEngineMCMC::PrintBestFitToStream(std::ofstream & ofi) const {
	if (GetBestFitParameters().empty()) {
		ofi << "No best fit information available." << std::endl << std::endl;
		return;
	}
		
	ofi << " Best Fit Results" << std::endl
			<< " ===========================" << std::endl
			<< " Log of the maximum posterior: " << GetLogMaximum() << std::endl
			<< " List of parameters and global mode:" << std::endl;

	for (unsigned i = 0; i < GetNVariables(); ++i) {
		ofi << TString::Format(" (%d) %10s \"%*s\" : %.*f", i, GetVariable(i)->GetPrefix().data(),
													 GetMaximumParameterNameLength(), GetVariable(i)->GetName().data(),
													 GetVariable(i)->GetPrecision(),GetGlobalMode()[i]);
		if (i<GetNParameters() and GetParameter(i)->Fixed())
			ofi << " (fixed)";
		ofi << std::endl;
	}
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintMarginalizationToStream(std::ofstream & ofi) const {
	if (fMCMCFlagRun)
		ofi << " Results of the marginalization" << std::endl
				<< " ==============================" << std::endl;

	// give warning if MCMC did not converge
	if (MCMCGetNIterationsConvergenceGlobal()<=0 && fMCMCFlagRun)
		ofi << " WARNING: the Markov Chain did not converge!" << std::endl
				<< " Be cautious using the following results!" << std::endl
				<< std::endl;
	
	ofi << " List of parameters and properties of the marginalized" << std::endl
			<< " distributions:" << std::endl;

	for (unsigned i = 0; i < GetNVariables(); ++i) {
		if ( ! GetVariable(i)->FillH1())
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
		GetMarginalized(i) -> PrintToStream(ofi,"      ",prec,std::vector<double>(1,0.68));
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
		
		for (unsigned i = 0; i < GetNParameters(); ++i)
			if (GetParameter(i)->Fixed())
				ofi << TString::Format("  (%d) Parameter \"%*s\" : (fixed)", i, GetMaximumParameterNameLength(false), GetParameter(i)->GetName().data()) << std::endl;
			else {
				double eff = 0;
				for (unsigned j = 0; j < fMCMCEfficiencies.size(); ++j)
					eff += fMCMCEfficiencies[j][i] / fMCMCEfficiencies.size();
				ofi << TString::Format("  (%d) Parameter \"%*s\" : %5.2f %%", i, GetMaximumParameterNameLength(false), GetParameter(i)->GetName().data(),eff*100) << std::endl;
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
void BCEngineMCMC::PrintParameters(std::vector<double> const & P, void (*output)(const char *) ) const {
	if ( P.size() != GetNParameters() and P.size() != GetNVariables() )
		return;

	for (unsigned i = 0; i < P.size(); ++i)
		output(TString::Format("          %-*s :   % .*g", GetMaximumParameterNameLength(P.size()>GetNParameters()), GetVariable(i)->GetName().data(),GetVariable(i)->GetPrecision(),P[i]));
}

// ---------------------------------------------------------
unsigned BCEngineMCMC::PrintAllMarginalized(std::string filename, unsigned hdiv, unsigned vdiv) const {
	if (fH1Marginalized.empty() and fH2Marginalized.empty()) {
		BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : Marginalized distributions not stored.");
		return 0;
	}

	BCAux::ForceToBePDF(filename);
	if (filename.empty())
		return 0;

	// Find nonempty H1's
	std::vector<BCH1D *> h1;
	h1.reserve(GetNVariables());
	for (unsigned i = 0; i < GetNVariables(); ++i)
		if ( MarginalizedHistogramExists(i) ) {
			if (GetMarginalizedHistogram(i)->Integral()==0) { // histogram was never filled in range
				BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 1D Marginalized histogram for \"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data()));
				continue;
			}
			h1.push_back(GetMarginalized(i));
			if (!h1.back()) // BCH1D doesn't exist
				h1.pop_back();
			else
				h1.back() -> CopyOptions(*fBCH1DdrawingOptions);
		}

	// create vector for ordering 2D hists properly
	std::vector<std::pair<unsigned,unsigned> > H2Coords;
	H2Coords.reserve(GetNVariables()*GetNVariables()-1);
	// par vs par
	for (unsigned i=0; i<GetNParameters(); ++i)
		for (unsigned j=0; j<GetNParameters(); ++j)
			if (i!=j) H2Coords.push_back(std::make_pair(i,j));
	// obs vs par
	for (unsigned i=0; i<GetNParameters(); ++i)
		for (unsigned j=GetNParameters(); j<GetNVariables(); ++j)
			H2Coords.push_back(std::make_pair(i,j));
	// par vs obs
	for (unsigned i=GetNParameters(); i<GetNVariables(); ++i)
		for (unsigned j=0; j<GetNParameters(); ++j)
			H2Coords.push_back(std::make_pair(i,j));
	// obs vs obs
	for (unsigned i=GetNParameters(); i<GetNVariables(); ++i)
		for (unsigned j=GetNParameters(); j<GetNVariables(); ++j)
			if (i!=j) H2Coords.push_back(std::make_pair(i,j));

	// Find nonempty H2's
	std::vector<BCH2D *> h2;
	h2.reserve(H2Coords.size());
	
	for (unsigned k=0; k<H2Coords.size(); ++k) {
		unsigned i = H2Coords[k].first;
		unsigned j = H2Coords[k].second;
		if ( MarginalizedHistogramExists(i,j) ) {
			if (fH2Marginalized[i][j]->Integral()==0) { // histogram was never filled in range
				BCLog::OutWarning(TString::Format("BCEngineMCMC::PrintAllMarginalized : 2D Marginalized histogram for \"%s\":\"%s\" is empty; printing is skipped.",GetVariable(i)->GetName().data(),GetVariable(i)->GetName().data()));
				continue;
			}
			h2.push_back(GetMarginalized(i,j));
			if (!h2.back()) // BCH2D doesn't exist
				h2.pop_back();
			else
				h2.back() -> CopyOptions(*fBCH2DdrawingOptions);
		}
	}
	
	if (h1.empty() and h2.empty()) {
		BCLog::OutWarning("BCEngineMCMC::PrintAllMarginalized : No marginalizations to print");
		return 0;
	}

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
		
		h1[i] -> Draw();
		
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
		
		h2[i] -> Draw();
		
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
unsigned BCEngineMCMC::PrintParameterPlot(std::string filename, int npar, double interval_content, std::vector<double> quantiles, bool rescale_ranges) const {

	BCAux::ForceToBePDF(filename);
	if (filename.empty())
		return 0;

	TCanvas * c_par = new TCanvas("c_parplot_init");
	c_par -> Print(Form("%s[",filename.data()));
	c_par -> cd();
	c_par -> SetTicky(1);
	c_par -> SetFrameLineWidth(0);
	c_par -> SetFrameLineColor(0);

	if (npar<=0) // all parameters on one page, all user-defined observables on the next
		npar = std::max<int> (GetNParameters(),GetNObservables());

	unsigned pages_printed = 0;

	// parameters first
	for (unsigned i = 0; i<GetNParameters(); i += npar)
		if (DrawParameterPlot(i,std::min<int>(npar,GetNParameters()-i), interval_content, quantiles,rescale_ranges)) {
			c_par->Print(filename.data());
			c_par->Clear();
			++pages_printed;
		}
	
	// then user-defined observables
	for (unsigned i = GetNParameters(); i<GetNVariables(); i += npar)
		if(DrawParameterPlot(i,std::min<int>(npar,GetNVariables()-i), interval_content, quantiles, rescale_ranges)) {
			c_par -> Print(filename.data());
			c_par -> Clear();
			++pages_printed;
		}

	c_par -> Print(Form("%s]",filename.data()));
	return pages_printed>0;
}

// ---------------------------------------------------------
bool BCEngineMCMC::DrawParameterPlot(unsigned i0, unsigned npar, double interval_content, std::vector<double> quantiles, bool rescale_ranges) const {

	// if npar==0, print all remaining observables
	unsigned i1 = (npar>0 && i0+npar<=GetNVariables()) ? i0+npar : GetNVariables();

	if (i1 <= i0) {
		BCLog::OutError(Form("BCSummaryTool::PrintParameterPlot : invalid parameter range [%d, %d)",i0,i1));
		return false;
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

	// store old ranges, and rescale, if rescaling
	std::vector<double> original_min;
	std::vector<double> original_max;
	original_min.reserve(i1-i0);
	original_max.reserve(i1-i0);
	if (rescale_ranges) {
		for (unsigned i=i0; i<i1; ++i) {
			if (i<fMCMCStatistics_AllChains.minimum.size() and std::isfinite(fMCMCStatistics_AllChains.minimum[i])) {
				original_min.push_back(GetVariable(i)->GetLowerLimit());
				GetVariable(i) -> SetLowerLimit(fMCMCStatistics_AllChains.minimum[i]);
			}
			if (i<fMCMCStatistics_AllChains.maximum.size() and std::isfinite(fMCMCStatistics_AllChains.maximum[i])) {
				original_max.push_back(GetVariable(i)->GetUpperLimit());
				GetVariable(i) -> SetUpperLimit(fMCMCStatistics_AllChains.maximum[i]);
			}
		}
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

		if (i<GetNParameters() and GetParameter(i)->Fixed())
			continue;
		
		// Global Mode
		x_i_bf.push_back(i);
		global_mode.push_back(GetVariable(i)->PositionInRange(GetGlobalMode()[i]));
		mean.push_back(std::numeric_limits<double>::infinity());
		rms.push_back(0);

		if (i < fH1Marginalized.size() and fH1Marginalized[i]) {
			BCH1D * bch1d_temp = GetMarginalized(i);
			if (bch1d_temp) {
				x_i.push_back(i);

				// quantiles
				x_quantiles.insert(x_quantiles.end(),quantiles.size(),i);
				double q[quantiles.size()];
				bch1d_temp -> GetHistogram() -> GetQuantiles(quantiles.size(),q,&quantiles[0]);
				for (unsigned j = 0; j < quantiles.size(); ++j)
					quantile_vals.push_back(GetVariable(i)->PositionInRange(q[j]));
			
				local_mode.push_back(GetVariable(i)->PositionInRange(bch1d_temp->GetLocalMode(0)));
				mean.back() = GetVariable(i)->PositionInRange(bch1d_temp->GetHistogram()->GetMean());
				rms.back() = bch1d_temp->GetHistogram()->GetRMS()/GetVariable(i)->GetRangeWidth();
			
				// smallest interval
				BCH1D::BCH1DSmallestInterval SI = bch1d_temp->GetSmallestIntervals(interval_content);
				if (SI.intervals.empty()) {
					interval_lo.push_back(0);
					interval_hi.push_back(0);
				} else {
					interval_lo.push_back(fabs(SI.intervals.front().mode-SI.intervals.front().xmin)/GetVariable(i)->GetRangeWidth());
					interval_hi.push_back(fabs(SI.intervals.front().mode-SI.intervals.front().xmax)/GetVariable(i)->GetRangeWidth());
				}
			}
		}

		// use chain statistics if they exist:
		if (i<fMCMCStatistics_AllChains.mean.size() and std::isfinite(fMCMCStatistics_AllChains.mean[i]))
			mean.back() = GetVariable(i) -> PositionInRange(fMCMCStatistics_AllChains.mean[i]);
		if (i<fMCMCStatistics_AllChains.variance.size() and std::isfinite(fMCMCStatistics_AllChains.variance[i]))
			rms.back() = sqrt(fMCMCStatistics_AllChains.variance[i])/GetVariable(i)->GetRangeWidth();
	}
	
	if (x_i.empty() and x_i_bf.empty())
		return false;

	/////////////////////////
	// Draw it all

	// Create, label, and draw axes
	TH2D * hist_axes = new TH2D(TString::Format("h2_axes_parplot_%s_%d_%d",GetSafeName().data(),i0,i1), "",//";;Scaled range [a.u.]",
															i1-i0, i0-0.5, i1-0.5, 10, -0.05+1e-3, 1.05-1e-3);
	hist_axes -> SetStats(kFALSE);
	hist_axes -> GetXaxis() -> SetAxisColor(0);
	hist_axes -> GetXaxis() -> SetLabelOffset(0.015);
	hist_axes -> GetXaxis() -> SetLabelSize(std::max<double>(0.01,0.05*std::min<double>(1,4./hist_axes->GetNbinsX())));
	hist_axes -> GetXaxis() -> SetTickLength(0.0);
	hist_axes -> GetYaxis() -> SetLabelSize(0);

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
	for (unsigned i = i0; i < i1; ++i)
		if (i<GetNParameters() and GetParameter(i)->Fixed()) {
			latex -> SetTextAlign(22);
			latex -> DrawLatex((double)i, 0.52, "Fixed at");
			latex -> DrawLatex((double)i, 0.47, Form("%.*g",GetVariable(i)->GetPrecision(),GetParameter(i)->GetFixedValue()));
		} else {
			latex -> SetTextAlign(21);
			latex -> DrawLatex((double)i,  1.015, Form("%+.*g", GetVariable(i)->GetPrecision(),GetVariable(i)->GetUpperLimit()));
			latex -> SetTextAlign(23);
			latex -> DrawLatex((double)i, -0.015, Form("%+.*g", GetVariable(i)->GetPrecision(),GetVariable(i)->GetLowerLimit()));
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
		graph_mean->SetMarkerStyle(21);
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

	gPad -> RedrawAxis();
	gPad -> Update();

	// restore ranges
	for (unsigned i=0; i<original_min.size(); ++i)
		GetVariable(i0+i) -> SetLowerLimit(original_min[i]);
	for (unsigned i=0; i<original_max.size(); ++i)
		GetVariable(i0+i) -> SetUpperLimit(original_max[i]);

	// no error
	return true;
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

	// vector of unfilled values:
	std::vector<std::pair<unsigned,unsigned> > unfilled;

	// fill histogram
	for (unsigned i = 0; i < GetNVariables(); ++i) {
		hist_corr->SetBinContent(i+1, GetNVariables()-i, 1);
		
		double var_i = (i<fMCMCStatistics_AllChains.variance.size()) ? fMCMCStatistics_AllChains.variance[i] : std::numeric_limits<double>::infinity();
		for (unsigned j = i+1; j < GetNVariables(); ++j) {
			double var_j = (j<fMCMCStatistics_AllChains.variance.size()) ? fMCMCStatistics_AllChains.variance[j] : std::numeric_limits<double>::infinity();
			double covar_ij = (i<fMCMCStatistics_AllChains.covariance.size() and j<fMCMCStatistics_AllChains.covariance[i].size()) ? fMCMCStatistics_AllChains.covariance[i][j] : std::numeric_limits<double>::infinity();
			double corr_ij = std::numeric_limits<double>::infinity();
			if (std::isfinite(covar_ij) and std::isfinite(var_i) and std::isfinite(var_j))
				corr_ij = covar_ij / sqrt(var_i*var_j);
			else if (i < fH2Marginalized.size() and j<fH2Marginalized[i].size() and fH2Marginalized[i][j])
				corr_ij = fH2Marginalized[i][j] -> GetCorrelationFactor();
			else if (j < fH2Marginalized.size() and i<fH2Marginalized[j].size() and fH2Marginalized[j][i])
				corr_ij = fH2Marginalized[j][i] -> GetCorrelationFactor();
			if (std::isfinite(corr_ij)) {
				hist_corr -> SetBinContent(i+1, GetNVariables()-j, corr_ij);
				hist_corr -> SetBinContent(j+1, GetNVariables()-i, corr_ij);
			} else {
				unfilled.push_back(std::make_pair(i,j));
			}
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

	 gStyle -> SetPalette(54);
	 gStyle -> SetPaintTextFormat("+.2g");
	 hist_corr -> GetZaxis() -> SetRangeUser(-1,1);
   hist_corr -> Draw("colz");
	 
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
	 	 // for (int j = 1; j <= hist_corr->GetNbinsY(); ++j) {
	 	 // 	 if (i==hist_corr->GetNbinsY()-j+1)
	 	 // 		 bcorr ->SetFillColor(kWhite);
	 	 // 	 else if (hist_corr->GetBinContent(i,j) > 0.5)
	 	 // 		 bcorr -> SetFillColor(kGreen);
	 	 // 	 else if (hist_corr->GetBinContent(i,j) >= -0.5)
	 	 // 		 bcorr -> SetFillColor(kYellow);
	 	 // 	 else 
	 	 // 		 bcorr -> SetFillColor(kRed);
	 	 // 	 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(i), hist_corr->GetYaxis()->GetBinLowEdge(j),
	 	 // 										hist_corr->GetXaxis()->GetBinUpEdge(i),  hist_corr->GetYaxis()->GetBinUpEdge(j));
	 	 // }
	 }

	 // write numbers in
	 hist_corr -> SetMarkerColor(kRed);
   hist_corr -> Draw("text same");

	 // Blank out empty squares
	 bcorr -> SetFillColor(kWhite);
	 for (unsigned i=0; i<unfilled.size(); ++i) {
		 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(unfilled[i].first+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-unfilled[i].second),
											hist_corr->GetXaxis()->GetBinUpEdge (unfilled[i].first+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-unfilled[i].second));
		 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(unfilled[i].second+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-unfilled[i].first),
											hist_corr->GetXaxis()->GetBinUpEdge (unfilled[i].second+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-unfilled[i].first));
	 }
   // for (unsigned i = 0; i < GetNVariables(); ++i)
	 // 	 for (unsigned j = i+1; j < GetNVariables(); ++j)
	 // 		 if (i>=fH2Marginalized.size() or j>=fH2Marginalized[i].size() or !fH2Marginalized[i][j]) {
	 // 			 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(i+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-j),
	 // 												hist_corr->GetXaxis()->GetBinUpEdge (i+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-j));
	 // 			 bcorr -> DrawBox(hist_corr->GetXaxis()->GetBinLowEdge(j+1), hist_corr->GetYaxis()->GetBinLowEdge(GetNVariables()-i),
	 // 												hist_corr->GetXaxis()->GetBinUpEdge (j+1), hist_corr->GetYaxis()->GetBinUpEdge (GetNVariables()-i));
	 // 		 }
	 
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
int BCEngineMCMC::PrintCorrelationPlot(const char * filename, bool include_observables) {
	// Array of indices for which any maginalizations were stored
	std::vector<unsigned> I;
	unsigned n = (include_observables) ? GetNVariables() : GetNParameters();
	for (unsigned i=0; i<n and i<fH1Marginalized.size(); ++i) {
		if (MarginalizedHistogramExists(i))
			I.push_back(i);
		else { 
			for (unsigned j=0; j<n and j<fH2Marginalized[i].size(); ++j)
				if (i!=j and MarginalizedHistogramExists(i,j)) {
					I.push_back(i);
					break;
				}
		}
	}
	
	if (I.empty())
		return 0;
	
	TCanvas * c = new TCanvas("c_correlation_plot");
	c->cd();
	
	double margin = 0.1;
	double padsize = (1 - 2*margin) / I.size();

	// array with pads holding the histograms
	std::vector<std::vector<TPad*> > pad (I.size(), std::vector<TPad*>(I.size(),NULL));
	
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
	TText * text_na = new TText();
	text_na -> SetTextFont(42);
	text_na -> SetTextAlign(22);
	text_na -> SetTextSize(8e-1/I.size());
	text_na -> SetTextColor(kGray);

	// drawing all histograms
	for (unsigned i = 0; i < I.size(); ++i) {
		xlow = i*padsize + margin;
		xup = xlow + padsize;

		for (unsigned j = 0; j < I.size(); ++j) {
			yup = 1. - j*padsize - margin;
			ylow = yup - padsize;

			// preparing the pad
			pad[i][j] =  new TPad(TString::Format("pad_correlation_plots_%d_%d",i,j), "", xlow, ylow, xup, yup);
			pad[i][j] -> SetMargin(marginleft,marginright,marginbottom,margintop);
			pad[i][j] -> SetFillColor(kWhite);
			pad[i][j] -> Draw();
			pad[i][j] -> cd();

			// get the histogram
			BCHistogramBase * bh = 0;

			if (i==j)
				bh = GetMarginalized(I[i]);
			else
				bh = MarginalizedHistogramExists(I[i],I[j]) ? GetMarginalized(I[i],I[j]) : NULL;
			
			if (bh) {
				
				bh -> GetHistogram() -> GetXaxis() -> SetLabelSize(0);
				bh -> GetHistogram() -> GetYaxis() -> SetLabelSize(0);
				bh -> GetHistogram() -> GetXaxis() -> SetTitleSize(0);
				bh -> GetHistogram() -> GetYaxis() -> SetTitleSize(0);

				if (bh->GetHistogram()->GetDimension()==1)
					bh -> CopyOptions(*fBCH1DdrawingOptions);
				else if (bh->GetHistogram()->GetDimension()==2)
					bh -> CopyOptions(*fBCH2DdrawingOptions);
				bh -> SetDrawLegend(false);
				bh -> SetStats(false);
				bh -> Draw();

			} else if (!(MarginalizedHistogramExists(I[j],I[i])) and I[i]>=I[j]) { // if the histogram is not available, draw "N/A"
				// pad[i][j] -> SetFillColor(kWhite);
				text_na -> DrawText(.5,.5,"N/A");
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
		BCLog::OutError(TString::Format("BCEngineMCMC::PrintParameterLatex : Couldn't open file %s",filename));
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
															 texwidth, prec, bch1d->GetHistogram()->GetMean(),
															 texwidth, prec, bch1d->GetHistogram()->GetRMS(),
															 texwidth, prec, GetGlobalMode()[i],
															 texwidth, prec, bch1d -> GetLocalMode(0),
															 texwidth, prec, bch1d -> GetMedian(),
															 texwidth, prec, bch1d -> GetQuantile(0.16),
															 texwidth, prec, bch1d -> GetQuantile(0.84))
						<< std::endl;
			
			// marginalization does not exist
			else
				ofi << TString::Format("        %*s & %*s & %*.*g & %*s & %*s & %*s & %*s\\\\\n",
															 texwidth, blank,
															 texwidth, blank,
															 texwidth, prec, GetGlobalMode()[i],
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
	if (n < 100000)
		return 10000;
	return 100000;
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

// ---------------------------------------------------------
BCEngineMCMC::MCMCStatistics::MCMCStatistics(const BCEngineMCMC::MCMCStatistics & other)
	: n_samples(other.n_samples)
	, mean(other.mean)
	, variance(other.variance)
	, covariance(other.covariance)
	, minimum(other.minimum)
	, maximum(other.maximum)
	, probability_mean(other.probability_mean)
	, probability_variance(other.probability_variance)
	, probability_mode(other.probability_mode)
	, mode(other.mode)
	, n_samples_efficiency(other.n_samples_efficiency)
	, efficiency(other.efficiency)
{
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCStatistics::MCMCStatistics(unsigned n_par, unsigned n_obs)
	:	n_samples(0)
	, mean(n_par+n_obs,0)
	, variance(mean.size(),0)
	, covariance(mean.size(),std::vector<double>(mean.size(),0))
	, minimum(mean.size(),+std::numeric_limits<double>::infinity())
	, maximum(mean.size(),-std::numeric_limits<double>::infinity())
	, probability_mean(0)
	, probability_variance(0)
	, probability_mode(-std::numeric_limits<double>::infinity())
	, mode(mean.size(),0)
	, n_samples_efficiency(0)
	, efficiency(n_par,0.)
{
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCStatistics::Clear(bool clear_mode, bool clear_efficiency) {
	n_samples = 0;
	mean.clear();
	variance.clear();
	covariance.clear();
	minimum.clear();
	maximum.clear();
	probability_mean = 0;
	probability_variance = 0;
	if (clear_mode) {
		probability_mode = -std::numeric_limits<double>::infinity();
		mode.clear();
	}
	if (clear_efficiency) {
		n_samples_efficiency = 0;
		efficiency.clear();
	}
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCStatistics::Init(unsigned n_par, unsigned n_obs) {
	n_samples = 0;
	mean.assign(n_par+n_obs,0);
	variance.assign(mean.size(),0);
	covariance.assign(mean.size(),std::vector<double>(mean.size(),0));
	minimum.assign(mean.size(),+std::numeric_limits<double>::infinity());
	maximum.assign(mean.size(),-std::numeric_limits<double>::infinity());
	probability_mean = 0;
	probability_variance = 0;
	probability_mode = -std::numeric_limits<double>::infinity();
	mode.assign(mean.size(),0);
	n_samples_efficiency = 0;
	efficiency.assign(n_par,0.);
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCStatistics::Reset(bool reset_mode, bool reset_efficiency) {
	n_samples = 0;
	mean.assign(mean.size(),0);
	variance.assign(variance.size(),0);
	covariance.assign(covariance.size(),std::vector<double>(covariance.front().size(),0));
	minimum.assign(minimum.size(),+std::numeric_limits<double>::infinity());
	maximum.assign(maximum.size(),-std::numeric_limits<double>::infinity());
	probability_mean = 0;
	probability_variance = 0;
	if (reset_mode) {
		probability_mode = -std::numeric_limits<double>::infinity();
		mode.assign(mode.size(),0);
	}
	if (reset_efficiency) {
		efficiency.assign(efficiency.size(),0);
		n_samples_efficiency = 0;
	}
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCStatistics::ResetEfficiencies() {
	efficiency.assign(efficiency.size(),0);
	n_samples_efficiency = 0;
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCStatistics & BCEngineMCMC::MCMCStatistics::operator=(const BCEngineMCMC::MCMCStatistics & rhs) {
	n_samples = rhs.n_samples;
	mean = rhs.mean;
	variance = rhs.variance;
	covariance = rhs.covariance;
	minimum = rhs.minimum;
	maximum = rhs.maximum;
	probability_mean = rhs.probability_mean;
	probability_variance = rhs.probability_variance;
	probability_mode = rhs.probability_mode;
	mode = rhs.mode;
	n_samples_efficiency = rhs.n_samples_efficiency;
	efficiency = rhs.efficiency;
	return *this;
}

// ---------------------------------------------------------
void BCEngineMCMC::MCMCStatistics::Update(double prob, const std::vector<double> & par, const std::vector<double> & obs) {
	if (mean.size() != par.size() + obs.size())
		return;

	// increment number of samples
	++n_samples;

	// check mode
	if (prob>probability_mode) {
		for (unsigned i=0; i<par.size(); ++i)
			mode[i] = par[i];
		for (unsigned i=0; i<obs.size(); ++i)
			mode[i+par.size()] = obs[i];
		probability_mode = prob;
	}

	// update probability mean and variance
	double prob_delta = prob - probability_mean;
	probability_mean += prob_delta/n_samples;
	probability_variance += (n_samples > 1) ? prob_delta*prob_delta/n_samples - probability_variance/(n_samples-1) : 0;

	// update parameter means and (co)variances, and maximums and minimums:

	// vector to store difference from current mean
	std::vector<double> delta(mean.size(),0);

	// loop over values
	for (unsigned i=0; i<mean.size(); ++i) {
		// get value
		double x = (i<par.size()) ? par[i] : obs[i-par.size()];
		// store difference to current mean
		delta[i] = x - mean[i];
		// update mean
		mean[i] += delta[i]/n_samples;
		// update variance
		variance[i] += (n_samples > 1) ? delta[i]*delta[i]/n_samples - variance[i]/(n_samples-1.) : 0;
		// update minimum
		if (x<minimum[i])
			minimum[i] = x;
		// update maximum
		if (x>maximum[i])
			maximum[i] = x;
	}
	// update covariances
	if (n_samples>1) {
		for (unsigned i=0; i<mean.size(); ++i)
			for (unsigned j=i; j<mean.size(); ++j)
				covariance[i][j] += delta[i]*delta[j]/n_samples - covariance[i][j]/(n_samples-1);
	}
}

// ---------------------------------------------------------
BCEngineMCMC::MCMCStatistics & BCEngineMCMC::MCMCStatistics::operator+=(const BCEngineMCMC::MCMCStatistics & rhs) {
	// if rhs is empty
	if (rhs.n_samples==0)
		return *this;

	// if this is empty
	if (n_samples==0)
		return operator=(rhs);

	if (mean.size() != rhs.mean.size())
		return *this;

	// check mode:
	if (rhs.probability_mode > probability_mode) {
		probability_mode = rhs.probability_mode;
		mode = rhs.mode;
	}

	double n = n_samples+rhs.n_samples;
	double N = (n>0) ? n_samples*rhs.n_samples/n : 0;

	probability_variance = (probability_variance*(n_samples-1) + rhs.probability_variance*(rhs.n_samples-1) + (rhs.probability_mean-probability_mean)*(rhs.probability_mean-probability_mean)*N)/(n-1);
	probability_mean = (n_samples*probability_mean + rhs.n_samples*rhs.probability_mean) / n;

	// parameter variables:
	for (unsigned i=0; i<mean.size(); ++i) {
		// check minimum
		if (rhs.minimum[i] < minimum[i])
			minimum[i] = rhs.minimum[i];
		// check maximum
		if (rhs.maximum[i] > maximum[i])
			maximum[i] = rhs.maximum[i];
		// combine variances
		variance[i] = (variance[i]*(n_samples-1)+rhs.variance[i]*(rhs.n_samples-1) + (rhs.mean[i]-mean[i])*(rhs.mean[i]-mean[i])*N)/(n-1);
		// combine covariances
		for (unsigned j=0; j<covariance[i].size(); ++j)
			covariance[i][j] = (covariance[i][j]*(n_samples-1)+rhs.covariance[i][j]*(rhs.n_samples-1) + (rhs.mean[i]-mean[i])*(rhs.mean[j]-mean[j])*N)/(n-1);
		// combine means
		mean[i] = (n_samples*mean[i] + rhs.n_samples*rhs.mean[i]) / n;
	}
	// combine n_samples
	n_samples = n;
	
	// combine efficiencies
	double n_eff = n_samples_efficiency + rhs.n_samples_efficiency;
	if (n_eff>0)
		for (unsigned i=0; i<efficiency[i]; ++i)
			efficiency[i] = (n_samples_efficiency*efficiency[i] + rhs.n_samples_efficiency*rhs.efficiency[i]) / (n_eff);

	// combine efficiency samples
	n_samples_efficiency = n_eff;

	return *this;
}
