/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <config.h>

#include "BCEngineMCMC.h"

#include "BCAux.h"
#include "BCGaussianPrior.h"
#include "BCH1D.h"
#include "BCH2D.h"
#include "BCMath.h"
#include "BCPrior.h"
#include "BCSplitGaussianPrior.h"
#include "BCTF1LogPrior.h"
#include "BCTF1Prior.h"
#include "BCTH1Prior.h"
#include "BCVariable.h"
#include "config.h"

#include <TCanvas.h>
#include <TDecompChol.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSeqCollection.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>

#if THREAD_PARALLELIZATION
#include <omp.h>
#endif

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const std::string& name)
    : fMCMCNIterationsConvergenceGlobal(-1),
      fMCMCFlagWriteChainToFile(false),
      fMCMCFlagWritePreRunToFile(false),
      fMCMCOutputFile(0),
      fMCMCOutputFilename(""),
      fMCMCOutputFileOption(""),
      fMCMCScaleFactorLowerLimit(0),
      fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max()),
      fMultivariateCovarianceUpdates(0),
      fMultivariateCovarianceUpdateLambda(0.5),
      fMultivariateEpsilon(0.05),
      fMultivariateScaleMultiplier(1.5),
      fMCMCEfficiencyMin(0.15),
      fMCMCEfficiencyMax(0.35),
      fInitialPositionScheme(BCEngineMCMC::kInitRandomUniform),
      fInitialPositionAttemptLimit(100),
      fMCMCProposeMultivariate(true),
      fMCMCProposalFunctionDof(1.0),
      fMCMCPhase(BCEngineMCMC::kUnsetPhase),
      fCorrectRValueForSamplingVariability(false),
      fMCMCRValueParametersCriterion(1.1),
      fMCMCTree(NULL),
      fMCMCTreeLoaded(false),
      fMCMCTreeReuseObservables(true),
      fParameterTree(NULL),
      fRescaleHistogramRangesAfterPreRun(false),
      fHistogramRescalePadding(0.1)
{
    SetName(name);
    SetPrecision(BCEngineMCMC::kMedium);
    SetRandomSeed(0);
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const std::string& filename, const std::string& name, bool loadObservables)
    : fMCMCNIterationsConvergenceGlobal(-1),
      fMCMCFlagWriteChainToFile(false),
      fMCMCFlagWritePreRunToFile(false),
      fMCMCOutputFile(NULL),
      fMCMCOutputFilename(""),
      fMCMCOutputFileOption(""),
      fMCMCScaleFactorLowerLimit(0),
      fMCMCScaleFactorUpperLimit(std::numeric_limits<double>::max()),
      fMultivariateCovarianceUpdates(0),
      fMultivariateCovarianceUpdateLambda(0.5),
      fMultivariateEpsilon(5.e-2),
      fMultivariateScaleMultiplier(1.5),
      fMCMCEfficiencyMin(0.15),
      fMCMCEfficiencyMax(0.35),
      fInitialPositionScheme(BCEngineMCMC::kInitRandomUniform),
      fInitialPositionAttemptLimit(100),
      fMCMCProposeMultivariate(true),
      fMCMCProposalFunctionDof(1.0),
      fMCMCPhase(BCEngineMCMC::kUnsetPhase),
      fCorrectRValueForSamplingVariability(false),
      fMCMCRValueParametersCriterion(1.1),
      fMCMCTree(NULL),
      fMCMCTreeLoaded(false),
      fMCMCTreeReuseObservables(true),
      fParameterTree(NULL),
      fRescaleHistogramRangesAfterPreRun(false),
      fHistogramRescalePadding(0.1)
{
    SetName(name);
    SetPrecision(BCEngineMCMC::kMedium);
    SetRandomSeed(0);
    LoadMCMC(filename, "", "", loadObservables);
}

// ---------------------------------------------------------
BCEngineMCMC::BCEngineMCMC(const BCEngineMCMC& other)
    : fMCMCThreadLocalStorage(other.fMCMCThreadLocalStorage),
      fChainIndex(other.fChainIndex),
      fName(other.fName),
      fSafeName(other.fSafeName),
      fParameters(other.fParameters),
      fObservables(other.fObservables),
      fMCMCNChains(other.fMCMCNChains),
      fMCMCNLag(other.fMCMCNLag),
      fMCMCCurrentIteration(other.fMCMCCurrentIteration),
      fMCMCNIterationsPreRunCheck(other.fMCMCNIterationsPreRunCheck),
      fMCMCPreRunCheckClear(other.fMCMCPreRunCheckClear),
      fMCMCNIterationsConvergenceGlobal(other.fMCMCNIterationsConvergenceGlobal),
      fMCMCNIterationsPreRunMax(other.fMCMCNIterationsPreRunMax),
      fMCMCNIterationsRun(other.fMCMCNIterationsRun),
      fMCMCNIterationsPreRunMin(other.fMCMCNIterationsPreRunMin),
      fMCMCFlagWriteChainToFile(other.fMCMCFlagWriteChainToFile),
      fMCMCFlagWritePreRunToFile(other.fMCMCFlagWritePreRunToFile),
      fMCMCOutputFile(other.fMCMCOutputFile),
      fMCMCOutputFilename(other.fMCMCOutputFilename),
      fMCMCOutputFileOption(other.fMCMCOutputFileOption),
      fMCMCScaleFactorLowerLimit(other.fMCMCScaleFactorLowerLimit),
      fMCMCScaleFactorUpperLimit(other.fMCMCScaleFactorUpperLimit),
      fMCMCProposalFunctionScaleFactor(other.fMCMCProposalFunctionScaleFactor),
      fMCMCInitialScaleFactors(other.fMCMCInitialScaleFactors),
      fMultivariateProposalFunctionCovariance(other.fMultivariateProposalFunctionCovariance),
      fMultivariateProposalFunctionCholeskyDecomposition(other.fMultivariateProposalFunctionCholeskyDecomposition),
      fMultivariateCovarianceUpdates(other.fMultivariateCovarianceUpdates),
      fMultivariateCovarianceUpdateLambda(other.fMultivariateCovarianceUpdateLambda),
      fMultivariateEpsilon(other.fMultivariateEpsilon),
      fMultivariateScaleMultiplier(other.fMultivariateScaleMultiplier),
      fMCMCInitialPosition(other.fMCMCInitialPosition),
      fMCMCEfficiencyMin(other.fMCMCEfficiencyMin),
      fMCMCEfficiencyMax(other.fMCMCEfficiencyMax),
      fInitialPositionScheme(other.fInitialPositionScheme),
      fInitialPositionAttemptLimit(other.fInitialPositionAttemptLimit),
      fMCMCProposeMultivariate(other.fMCMCProposeMultivariate),
      fMCMCProposalFunctionDof(other.fMCMCProposalFunctionDof),
      fMCMCPhase(other.fMCMCPhase),
      fMCMCStates(other.fMCMCStates),
      fMCMCStatistics(other.fMCMCStatistics),
      fMCMCStatistics_AllChains(other.fMCMCStatistics_AllChains),
      fCorrectRValueForSamplingVariability(other.fCorrectRValueForSamplingVariability),
      fMCMCRValueParametersCriterion(other.fMCMCRValueParametersCriterion),
      fMCMCRValueParameters(other.fMCMCRValueParameters),
      fRandom(other.fRandom),
      fRequestedH2(other.fRequestedH2),
      fMCMCTree(NULL),
      fMCMCTreeLoaded(false),
      fMCMCTreeReuseObservables(other.fMCMCTreeReuseObservables),
      fMCMCTree_State(other.fMCMCTree_State),
      fParameterTree(NULL),
      fLocalModes(other.fLocalModes),
      fBCH1DdrawingOptions(other.fBCH1DdrawingOptions),
      fBCH2DdrawingOptions(other.fBCH2DdrawingOptions),
      fRescaleHistogramRangesAfterPreRun(other.fRescaleHistogramRangesAfterPreRun),
      fHistogramRescalePadding(other.fHistogramRescalePadding),
      fObjectTrash(other.fObjectTrash)
{
    // set again in case user overloads the setter to create custom structures
    SetNChains(other.fMCMCNChains);

    CloneMarginals(other);
}

BCEngineMCMC& BCEngineMCMC::operator=(const BCEngineMCMC& other)
{
    if (this == &other) {
        return *this;
    }

    DeleteMarginals();
    CloseOutputFile();

    // exception safety
    try {
        fMCMCThreadLocalStorage = other.fMCMCThreadLocalStorage;
        fChainIndex = other.fChainIndex;
        fName = other.fName;
        fSafeName = other.fSafeName;
        fParameters = other.fParameters;
        fObservables = other.fObservables;
        SetNChains(other.fMCMCNChains);
        fMCMCNLag = other.fMCMCNLag;
        fMCMCCurrentIteration = other.fMCMCCurrentIteration;
        fMCMCNIterationsPreRunCheck = other.fMCMCNIterationsPreRunCheck;
        fMCMCPreRunCheckClear = other.fMCMCPreRunCheckClear;
        fMCMCNIterationsConvergenceGlobal = other.fMCMCNIterationsConvergenceGlobal;
        fMCMCNIterationsPreRunMax = other.fMCMCNIterationsPreRunMax;
        fMCMCNIterationsRun = other.fMCMCNIterationsRun;
        fMCMCNIterationsPreRunMin = other.fMCMCNIterationsPreRunMin;
        fMCMCFlagWriteChainToFile = other.fMCMCFlagWriteChainToFile;
        fMCMCFlagWritePreRunToFile = other.fMCMCFlagWritePreRunToFile;
        fMCMCOutputFile = other.fMCMCOutputFile;
        fMCMCOutputFilename = other.fMCMCOutputFilename;
        fMCMCOutputFileOption = other.fMCMCOutputFileOption;
        fMCMCScaleFactorLowerLimit = other.fMCMCScaleFactorLowerLimit;
        fMCMCScaleFactorUpperLimit = other.fMCMCScaleFactorUpperLimit;
        fMCMCProposalFunctionScaleFactor = other.fMCMCProposalFunctionScaleFactor;
        fMCMCInitialScaleFactors = other.fMCMCInitialScaleFactors;
        fMultivariateProposalFunctionCovariance = other.fMultivariateProposalFunctionCovariance;
        fMultivariateProposalFunctionCholeskyDecomposition = other.fMultivariateProposalFunctionCholeskyDecomposition;
        fMultivariateCovarianceUpdates = other.fMultivariateCovarianceUpdates;
        fMultivariateCovarianceUpdateLambda = other.fMultivariateCovarianceUpdateLambda;
        fMultivariateEpsilon = other.fMultivariateEpsilon;
        fMultivariateScaleMultiplier = other.fMultivariateScaleMultiplier;
        fMCMCInitialPosition = other.fMCMCInitialPosition;
        fMCMCEfficiencyMin = other.fMCMCEfficiencyMin;
        fMCMCEfficiencyMax = other.fMCMCEfficiencyMax;
        fInitialPositionScheme = other.fInitialPositionScheme;
        fInitialPositionAttemptLimit = other.fInitialPositionAttemptLimit;
        fMCMCProposeMultivariate = other.fMCMCProposeMultivariate;
        fMCMCPhase = other.fMCMCPhase;
        fMCMCStates = fMCMCStates;
        fMCMCStatistics = other.fMCMCStatistics;
        fMCMCStatistics_AllChains = other.fMCMCStatistics_AllChains;
        fCorrectRValueForSamplingVariability = other.fCorrectRValueForSamplingVariability;
        fMCMCRValueParametersCriterion = other.fMCMCRValueParametersCriterion;
        fMCMCRValueParameters = other.fMCMCRValueParameters;
        fRandom = other.fRandom;
        fRequestedH2 = other.fRequestedH2;
        fMCMCTreeReuseObservables = other.fMCMCTreeReuseObservables;
        fMCMCTree_State = other.fMCMCTree_State;
        fLocalModes = other.fLocalModes;
        fBCH1DdrawingOptions = other.fBCH1DdrawingOptions;
        fBCH2DdrawingOptions = other.fBCH2DdrawingOptions;
        fRescaleHistogramRangesAfterPreRun = other.fRescaleHistogramRangesAfterPreRun;
        fHistogramRescalePadding = other.fHistogramRescalePadding;
        fObjectTrash = other.fObjectTrash;

        // don't create file!

        CloneMarginals(other);
    } catch (...) {
        // leave object in sane state but otherwise don't know what to do with exception
        DeleteMarginals();
        throw;
    }

    return *this;
}

// ---------------------------------------------------------
void BCEngineMCMC::CloneMarginals(const BCEngineMCMC& other)
{
    fH1Marginalized = std::vector<TH1*>(other.fH1Marginalized.size(), NULL);
    for (unsigned i = 0; i < other.fH1Marginalized.size(); ++i)
        if (other.fH1Marginalized[i])
            fH1Marginalized[i] = BCAux::OwnClone(other.fH1Marginalized[i]);

    if (!other.fH2Marginalized.empty() && !other.fH2Marginalized.front().empty()) {
        fH2Marginalized = std::vector<std::vector<TH2*> > (other.fH2Marginalized.size(), std::vector<TH2*>(other.fH2Marginalized.front().size(), NULL));
        for (unsigned i = 0; i < other.fH2Marginalized.size(); ++i) {
            fH2Marginalized[i].assign(other.fH2Marginalized[i].size(), NULL);
            for (unsigned j = 0; j < other.fH2Marginalized[i].size(); ++j)
                if (other.fH2Marginalized[i][j])
                    fH2Marginalized[i][j] = BCAux::OwnClone(other.fH2Marginalized[i][j]);
        }
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::DeleteMarginals()
{
    // delete 1-d marginalized distributions
    for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
        delete fH1Marginalized[i];

    // delete 2-d marginalized distributions
    for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
        for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
            delete fH2Marginalized[i][j];
}

// ---------------------------------------------------------
BCEngineMCMC::~BCEngineMCMC()
{
    DeleteMarginals();
    CloseOutputFile();
}

// ---------------------------------------------------------
void BCEngineMCMC::SetName(const std::string& name)
{
    fName = name;
    fSafeName = BCAux::SafeName(name);
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPrior(unsigned index, TF1& f, bool logL)
{
    if (logL)
        fParameters.At(index).SetPrior(new BCTF1LogPrior(f));
    else
        fParameters.At(index).SetPrior(new BCTF1Prior(f));
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPrior(unsigned index, TH1& h, bool interpolate)
{
    fParameters.At(index).SetPrior(new BCTH1Prior(h, interpolate));
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPriorGauss(unsigned index, double mean, double sigma)
{
    fParameters.At(index).SetPrior(new BCGaussianPrior(mean, sigma));
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPriorGauss(unsigned index, double mean, double sigma_below, double sigma_above)
{
    fParameters.At(index).SetPrior(new BCSplitGaussianPrior(mean, sigma_below, sigma_above));
}

// ---------------------------------------------------------
void BCEngineMCMC::SetNLag(unsigned n)
{
    if (n == 0) {
        BCLog::OutError("Invalid lag = 0 given. Set to lag = 1");
        fMCMCNLag = 1;
    } else {
        fMCMCNLag = n;
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPrecision(BCEngineMCMC::Precision precision)
{

    // don't clear means, variances etc. during prerun
    fMCMCPreRunCheckClear = 0;

    // take every iteration
    fMCMCNLag = 1;

    switch (precision) {

        case BCEngineMCMC::kLow:
            SetNChains(1);
            fMCMCNIterationsPreRunCheck           = 500;
            fMCMCNIterationsPreRunMin             = 1500;
            fMCMCNIterationsPreRunMax             = 10000;
            fMCMCNIterationsRun                   = 10000;
            break;

        case BCEngineMCMC::kQuick:
            SetNChains(2);
            fMCMCNIterationsPreRunCheck           = 500;
            fMCMCNIterationsPreRunMin             = 1500;
            fMCMCNIterationsPreRunMax             = 10000;
            fMCMCNIterationsRun                   = 10000;
            break;

        case BCEngineMCMC::kMedium:
            SetNChains(4);
            fMCMCNIterationsPreRunCheck           = 500;
            fMCMCNIterationsPreRunMin             = 1500;
            fMCMCNIterationsPreRunMax             = 100000;
            fMCMCNIterationsRun                   = 100000;
            break;

        case BCEngineMCMC::kHigh:
            SetNChains(8);
            fMCMCNIterationsPreRunCheck           = 1000;
            fMCMCNIterationsPreRunMin             = 5000;
            fMCMCNIterationsPreRunMax             = 1000000;
            fMCMCNIterationsRun                   = 1000000;
            break;

        case BCEngineMCMC::kVeryHigh:
            SetNChains(8);
            fMCMCNIterationsPreRunCheck           = 1000;
            fMCMCNIterationsPreRunMin             = 10000;
            fMCMCNIterationsPreRunMax             = 10000000;
            fMCMCNIterationsRun                   = 10000000;
            break;

        default:
            BCLOG_ERROR("Invalid precision");
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::SetPrecision(const BCEngineMCMC& other)
{
    SetNChains(other.fMCMCNChains);
    fMCMCNLag                             = other.fMCMCNLag;
    fMCMCNIterationsPreRunMin             = other.fMCMCNIterationsPreRunMin;
    fMCMCNIterationsPreRunMax             = other.fMCMCNIterationsPreRunMax;
    fMCMCNIterationsRun                   = other.fMCMCNIterationsRun;
    fMCMCNIterationsPreRunCheck           = other.fMCMCNIterationsPreRunCheck;
    fMCMCPreRunCheckClear                 = other.fMCMCPreRunCheckClear;
    fMCMCRValueParametersCriterion        = other.fMCMCRValueParametersCriterion;
    fMCMCEfficiencyMin                    = other.fMCMCEfficiencyMin;
    fMCMCEfficiencyMax                    = other.fMCMCEfficiencyMax;
    fMCMCProposeMultivariate     = other.fMCMCProposeMultivariate;
    fMCMCProposalFunctionDof  = other.fMCMCProposalFunctionDof;
    fMultivariateEpsilon  = other.fMultivariateEpsilon;
    fMultivariateScaleMultiplier = other.fMultivariateScaleMultiplier;
    fMultivariateCovarianceUpdateLambda = other.fMultivariateCovarianceUpdateLambda;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChainRun(bool flag)
{
    if (flag && fMCMCOutputFilename.empty())
        BCLog::OutError("BCEngineMCMC::WriteMarkovChainRun: First turn on output using WriteMarkovChain(filename, option, main_run, pre_run).");
    fMCMCFlagWriteChainToFile = flag;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChainPreRun(bool flag)
{
    if (flag && fMCMCOutputFilename.empty())
        BCLog::OutError("BCEngineMCMC::WriteMarkovChainPreRun: First turn on output using WriteMarkovChain(filename, option, main_run, pre_run).");
    fMCMCFlagWritePreRunToFile = flag;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarkovChain(const std::string& filename, const std::string& option, bool flag_run, bool flag_prerun)
{
    // if setting both false
    if (!flag_run && !flag_prerun)
        WriteMarkovChain(false);

    if (filename.empty()) {
        BCLog::OutError("BCEngineMCMC::WriteMarkovChain: You must specify the filename when turning on Markov chain output.");
        return WriteMarkovChain(false);
    }

    fMCMCOutputFilename = filename;
    fMCMCOutputFileOption = option;
    fMCMCFlagWriteChainToFile = flag_run;
    fMCMCFlagWritePreRunToFile = flag_prerun;
}

// --------------------------------------------------------
void BCEngineMCMC::WriteMarginalizedDistributions(const std::string& filename, const std::string& option, bool closeExistingFile)
{
    // remember current directory
    TDirectory* dir = gDirectory;

    // look to see if file is already open
    TSeqCollection* listOfFiles = gROOT->GetListOfFiles();
    TFile* fOut = NULL;
    for (int i = 0; i < listOfFiles->GetEntries(); ++i)
        if (listOfFiles->At(i) && filename.compare(listOfFiles->At(i)->GetName()) == 0) {
            fOut = dynamic_cast<TFile*>(listOfFiles->At(i));
            break;
        }
    if (fOut) {
        if (option.compare("RECREATE") == 0) {
            // if recreating, close current file, which will be overwritten
            BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: File already open. option \"RECREATE\" will now overwrite it!");
            if (fOut->IsWritable())
                fOut->Write(0, TObject::kWriteDelete);
            fOut->Close();
            fOut = NULL;
            closeExistingFile = true; // existing file closed; close newly opened file
        } else if (option.compare("UPDATE") == 0 && !fOut->IsWritable()) {
            BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: File already open but not in readable mode.");
            return;
        }
    } else {
        closeExistingFile = true; // no pre-open file found; close newly opened file
    }

    // else open file
    if (!fOut)
        fOut = TFile::Open(filename.c_str(), option.c_str());

    if (!fOut) {
        BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: Could not open output file.");
        return;
    }

    if (!fOut->IsWritable()) {
        BCLog::OutError("BCEngineMCMC::WriteMarginalizedDistributions: File must be opened in writeable mode.");
        return;
    }

    // write histograms
    for (unsigned i = 0; i < GetNVariables(); ++i) {
        if (MarginalizedHistogramExists(i))
            fOut->WriteTObject(GetMarginalizedHistogram(i));
        for (unsigned j = 0; j < GetNVariables(); ++j)
            if (MarginalizedHistogramExists(i, j))
                fOut->WriteTObject(GetMarginalizedHistogram(i, j));
    }

    // Causes multiple instances of any tree that already exists in the root file.
    // fOut->Write();

    if (closeExistingFile)
        fOut->Close();

    // restore directory
    gDirectory = dir;
}

// --------------------------------------------------------
TH1* BCEngineMCMC::GetMarginalizedHistogram(unsigned index) const
{
    if (index >= fH1Marginalized.size()) {
        BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", index));
        return 0;
    }

    if (fH1Marginalized[index])
        return fH1Marginalized[index];

    // else output warning
    if (index < GetNVariables()) // Marginalization of model parameter
        BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram: marginal distribution not stored for %s %s", GetVariable(index).GetPrefix().data(), GetVariable(index).GetName().data()));

    return NULL;
}

// --------------------------------------------------------
TH2* BCEngineMCMC::GetMarginalizedHistogram(unsigned i, unsigned j) const
{
    if (i == j) {
        BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Called with identical indices %u.", i));
        return NULL;
    }

    if (i >= fH2Marginalized.size()) {
        BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", i));
        return NULL;
    }
    if (j >= fH2Marginalized[i].size()) {
        BCLog::OutError(Form("BCEngineMCMC::GetMarginalizedHistogram. Index %u out of range.", j));
        return NULL;
    }

    if (fH2Marginalized[i][j])
        return fH2Marginalized[i][j];

    if (i < GetNVariables() && j < GetNVariables())
        BCLog::OutWarning(Form("BCEngineMCMC::GetMarginalizedHistogram : marginal distribution not stored for %s %s vs %s %s",
                               GetVariable(i).GetPrefix().data(), GetVariable(i).GetName().data(),
                               GetVariable(j).GetPrefix().data(), GetVariable(j).GetName().data()));
    return NULL;
}

// --------------------------------------------------------
BCH1D BCEngineMCMC::GetMarginalized(unsigned index) const
{
    TH1* h = GetMarginalizedHistogram(index);

    BCH1D bch(h);

    if (bch.Valid() && index < GetBestFitParameters().size())
        // set global mode if available
        bch.SetGlobalMode(GetBestFitParameters()[index]);

    return bch;
}

// --------------------------------------------------------
BCH2D BCEngineMCMC::GetMarginalized(unsigned i, unsigned j) const
{
    TH2* h = NULL;

    if (!MarginalizedHistogramExists(i, j) && MarginalizedHistogramExists(j, i))
        h = BCAux::Transpose(GetMarginalizedHistogram(j, i));
    else
        h = GetMarginalizedHistogram(i, j);

    BCH2D bch(h);

    // set global mode if available
    if (bch.Valid() && i < GetBestFitParameters().size() && j < GetBestFitParameters().size())
        bch.SetGlobalMode(GetBestFitParameters()[i], GetBestFitParameters()[j]);

    return bch;
}

// --------------------------------------------------------
unsigned BCEngineMCMC::GetCurrentChain() const
{
    // serial case is the default
    static const unsigned defaultChain = 0;

    if (fChainIndex.empty())
        return defaultChain;

    int id = 0;
#if THREAD_PARALLELIZATION
    id = omp_get_thread_num();
#endif
    return fChainIndex.at(id);
}

// ---------------------------------------------------------
unsigned BCEngineMCMC::GetNIterationsPreRun() const
{
    if (fMCMCPhase == kUnsetPhase)
        return 0;
    BCLOG_WARNING(Form("global %d min %d max %d", fMCMCNIterationsConvergenceGlobal, int(fMCMCNIterationsPreRunMin), int(fMCMCNIterationsPreRunMax)));
    if (fMCMCNIterationsConvergenceGlobal > -1)
        return fMCMCNIterationsConvergenceGlobal;
    return fMCMCNIterationsPreRunMax;
}

// ---------------------------------------------------------
const std::vector<double>& BCEngineMCMC::GetBestFitParameters() const
{
    return fMCMCStatistics_AllChains.modepar;
}

// ---------------------------------------------------------
const std::vector<double>& BCEngineMCMC::GetBestFitParameterErrors() const
{
    if (fMCMCStatistics_AllChains.stderrpar.size() < GetNParameters())
        BCLOG_ERROR("Standard errors not available. Run Markov chain first");

    return fMCMCStatistics_AllChains.stderrpar;
}

// ---------------------------------------------------------
const std::vector<double>& BCEngineMCMC::GetLocalModes(bool force_recalculation)
{
    if (fLocalModes.empty() || force_recalculation) {
        fLocalModes.clear();
        for (unsigned i = 0; i < GetNVariables(); ++i)
            if (i < GetNParameters() && GetParameter(i).Fixed())
                fLocalModes.push_back(GetParameter(i).GetFixedValue());
            else if (fH1Marginalized[i]) {
                fLocalModes.push_back(fH1Marginalized[i]->GetBinCenter(fH1Marginalized[i]->GetMaximumBin()));
            } else {
                if (i < GetNParameters()) {
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

// --------------------------------------------------------
void BCEngineMCMC::SetInitialPositions(const std::vector<double>& x0)
{
    if (x0.size() != GetNParameters()) {
        BCLOG_ERROR(Form("#initial positions does not match #parameters: %u vs %u", unsigned(x0.size()), GetNParameters()));
        return;
    }
    fMCMCInitialPosition.clear();
    for (unsigned i = 0; i < GetNChains(); ++i)
        fMCMCInitialPosition.push_back(x0);
    SetInitialPositionScheme(BCEngineMCMC::kInitUserDefined);
}

// --------------------------------------------------------
void BCEngineMCMC::SetRandomSeed(unsigned seed)
{
    // set main generator
    fRandom.SetSeed(seed);

    // call once so return value of GetSeed() fixed
    fRandom.Rndm();

    SyncThreadStorage();

    // type conversion to avoid compiler warnings
    if (size_t(fMCMCNChains) != fMCMCThreadLocalStorage.size())
        BCLog::OutError(Form("#chains does not match #(thread local storages): %d vs %u",
                             fMCMCNChains, unsigned(fMCMCThreadLocalStorage.size())));
}

// --------------------------------------------------------
void BCEngineMCMC::InitializeMarkovChainTree(bool replacetree, bool replacefile)
{
    if (replacetree) {
        delete fMCMCTree;
        fMCMCTree = NULL;
        delete fParameterTree;
        fParameterTree = NULL;
    }
    if (replacefile) {
        if (fMCMCOutputFile)
            fMCMCOutputFile->Close();
        delete fMCMCOutputFile;
        fMCMCOutputFile = NULL;
    }

    // don't initialize a 2nd time
    if (fMCMCOutputFile && fMCMCTree && fParameterTree)
        return;

    TDirectory* dir = gDirectory;

    // create file
    if (!fMCMCOutputFile)
        fMCMCOutputFile = TFile::Open(fMCMCOutputFilename.c_str(), fMCMCOutputFileOption.c_str());
    // if failed
    if (!fMCMCOutputFile) {
        BCLog::OutError("BCEngineMCMC::InitializeMarkovChainTree: Could not create new file.");
        WriteMarkovChain(false);
        return;
    }
    // if file mode not writeable
    if (fMCMCOutputFile && !fMCMCOutputFile->IsWritable()) {
        BCLog::OutError("BCEngineMCMC::InitializeMarkovChainTree: File must be opened in a writeable mode.");
        delete fMCMCOutputFile;
        fMCMCOutputFile = 0;
        WriteMarkovChain(false);
        return;
    }

    // change to output file (directory)
    fMCMCOutputFile->cd();

    if (fMCMCTree)
        // Add existing MCMC tree to output file
        fMCMCOutputFile->Add(fMCMCTree);
    else {
        // or create new MCMC tree (under umbrella of output file)
        fMCMCTree = new TTree(TString::Format("%s_mcmc", GetSafeName().data()), TString::Format("%s_mcmc", GetSafeName().data()));
        fMCMCTree->Branch("Chain",          &fMCMCTree_Chain,                 "chain/i");
        fMCMCTree->Branch("Iteration",      &fMCMCTree_State.iteration,       "iteration/i");
        fMCMCTree->Branch("Phase",          &fMCMCPhase,                      "phase/I");
        fMCMCTree->Branch("LogProbability", &fMCMCTree_State.log_probability, "log_probability/D");
        fMCMCTree_State.parameters.assign(GetNParameters(), 0);
        for (unsigned j = 0; j < GetNParameters(); ++j) {
            fMCMCTree->Branch(GetParameter(j).GetSafeName().data(), &fMCMCTree_State.parameters[j], (GetParameter(j).GetSafeName() + "/D").data());
            fMCMCTree->SetAlias(TString::Format("Parameter%i", j), GetParameter(j).GetSafeName().data());
        }
        fMCMCTree_State.observables.assign(GetNObservables(), 0);
        for (unsigned j = 0; j < GetNObservables(); ++j) {
            fMCMCTree->Branch(GetObservable(j).GetSafeName().data(), &fMCMCTree_State.observables[j], (GetObservable(j).GetSafeName() + "/D").data());
            fMCMCTree->SetAlias(TString::Format("Observable%i", j), GetObservable(j).GetSafeName().data());
        }
    }

    if (fParameterTree)
        // add existing parameter tree to output file
        fMCMCOutputFile->Add(fParameterTree);
    else {
        // or create new parameter tree (under umbrella of output file)
        fParameterTree = new TTree(TString::Format("%s_parameters", GetSafeName().data()), TString::Format("%s_parameters", GetSafeName().data()));
        bool p_parameter, p_fill_1d, p_fill_2d, p_fixed;
        unsigned p_index, p_precision, p_nbins;
        char p_name[200], p_safename[200], p_latexname[200], p_unitstring[200];
        double p_lowerlimit, p_upperlimit, p_fixedvalue;
        fParameterTree->Branch("parameter", &p_parameter, "parameter/O");
        fParameterTree->Branch("index", &p_index, "index/i");
        fParameterTree->Branch("name", p_name, "name/C");
        fParameterTree->Branch("safe_name", p_safename, "safe_name/C");
        fParameterTree->Branch("latex_name", p_latexname, "latex_name/C");
        fParameterTree->Branch("unit_string", p_unitstring, "unit_string/C");
        fParameterTree->Branch("lower_limit", &p_lowerlimit, "lower_limit/D");
        fParameterTree->Branch("upper_limit", &p_upperlimit, "upper_limit/D");
        fParameterTree->Branch("precision", &p_precision, "precision/i");
        fParameterTree->Branch("nbins", &p_nbins, "nbins/i");
        fParameterTree->Branch("fill_1d", &p_fill_1d, "fill_1d/O");
        fParameterTree->Branch("fill_2d", &p_fill_2d, "fill_2d/O");
        fParameterTree->Branch("fixed", &p_fixed, "fixed/O");
        fParameterTree->Branch("fixed_value", &p_fixedvalue, "fixed_value/D");
        for (unsigned i = 0; i < GetNVariables(); ++i) {
            p_parameter  = (i < GetNParameters());
            p_index      = (p_parameter) ? i : i - GetNParameters();
            strcpy(p_name,      GetVariable(i).GetName().data());
            strcpy(p_safename,  GetVariable(i).GetSafeName().data());
            strcpy(p_latexname, GetVariable(i).GetLatexName().data());
            strcpy(p_unitstring, GetVariable(i).GetUnitString().data());
            p_lowerlimit = GetVariable(i).GetLowerLimit();
            p_upperlimit = GetVariable(i).GetUpperLimit();
            p_precision  = GetVariable(i).GetPrecision();
            p_nbins      = GetVariable(i).GetNbins();
            p_fill_1d    = GetVariable(i).FillH1();
            p_fill_2d    = GetVariable(i).FillH2();
            p_fixed      = p_parameter && GetParameter(i).Fixed();
            p_fixedvalue = (p_parameter) ? GetParameter(i).GetFixedValue() : 0;
            fParameterTree->Fill();
        }
        fParameterTree->AutoSave("SaveSelf");
        // fParameterTree->ResetBranchAddresses();
    }

    // return to old directory
    gDirectory = dir;
}

// --------------------------------------------------------
void BCEngineMCMC::UpdateParameterTree()
{
    if (!fParameterTree)
        return;

    unsigned nchains = GetNChains();
    std::vector<double> scale(GetNChains(), 0);
    std::vector<double> eff(GetNChains(), 0);

    // check for branch existences
    TBranch* b_nchains = fParameterTree->GetBranch("nchains");
    TBranch* b_scale   = fParameterTree->GetBranch("scale");

    // if nchains branch doesn't exist, create it
    if (b_nchains == 0)
        b_nchains = fParameterTree->Branch("nchains", &nchains, "nchains/i");
    // else set 0, so as not to fill
    else
        b_nchains = 0;

    // if scale branch doesn't exist, create it
    if (b_scale == 0)
        b_scale = fParameterTree->Branch("scale", &(scale.front()), TString::Format("scale[%d]/D", GetNChains()));
    // else set 0, so as not to fill
    else
        b_scale = 0;

    // create next effiency branch
    unsigned i = 0;
    while (fParameterTree->GetBranch(TString::Format("efficiency_%d", i)))
        ++i;
    TBranch* b_eff = fParameterTree->Branch(TString::Format("efficiency_%d", i), &(eff.front()), TString::Format("efficiency_%d[%d]/D", i, GetNChains()));

    for (unsigned n = 0; n < fParameterTree->GetEntries(); ++n) {
        if (b_nchains)
            b_nchains->Fill();

        for (unsigned j = 0; j < nchains; ++j) {
            scale[j] = (n < GetNParameters()) ? (fMCMCProposeMultivariate ? fMCMCProposalFunctionScaleFactor[j][0] : fMCMCProposalFunctionScaleFactor[j][n]) : -1;
            if (!fMCMCProposeMultivariate)
                eff[j]   = (n < GetNParameters()) ? fMCMCStatistics[j].efficiency[n] : -1;
            else if (!fMCMCStatistics[j].efficiency.empty())
                eff[j] = fMCMCStatistics[j].efficiency.front();
            else
                eff[j] = -1;
        }

        if (b_scale)
            b_scale->Fill();

        b_eff->Fill();
    }
    fParameterTree->AutoSave("SaveSelf");
    // fParameterTree->ResetBranchAddresses();
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidMCMCTree(TTree* tree, bool checkObservables) const
{
    if (!tree)
        return false;

    if (!(tree->GetBranch("Chain")))
        return false;

    if (!(tree->GetBranch("Phase")))
        return false;

    unsigned nvar = checkObservables ? GetNObservables() : GetNParameters();
    for (unsigned i = 0; i < nvar; ++i)
        if (!(tree->GetBranch(GetVariable(i).GetSafeName().data())))
            return false;

    return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::ValidParameterTree(TTree* tree) const
{
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

    return true;
}

// --------------------------------------------------------
void BCEngineMCMC::LoadParametersFromTree(TTree* partree, bool loadObservables)
{
    // absolutely necessary branches
    if (!partree->GetBranch("parameter"))
        throw std::runtime_error("BCEngineMCMC::LoadParametersFromTree: tree missing parameter branch");
    if (!partree->GetBranch("index"))
        throw std::runtime_error("BCEngineMCMC::LoadParametersFromTree: tree missing index branch");
    if (!partree->GetBranch("name"))
        throw std::runtime_error("BCEngineMCMC::LoadParametersFromTree: tree missing name branch");
    if (!partree->GetBranch("lower_limit"))
        throw std::runtime_error("BCEngineMCMC::LoadParametersFromTree: tree missing lower_limit branch");
    if (!partree->GetBranch("upper_limit"))
        throw std::runtime_error("BCEngineMCMC::LoadParametersFromTree: tree missing upper_limit branch");

    partree->ResetBranchAddresses();

    char p_name[200];
    double p_lowerlimit;
    double p_upperlimit;

    partree->SetBranchAddress("name", p_name);
    partree->SetBranchAddress("lower_limit", &p_lowerlimit);
    partree->SetBranchAddress("upper_limit", &p_upperlimit);

    // not entirely necessary branches
    char p_latexname[200] = "";
    unsigned p_precision = 6;
    unsigned p_nbins = 100;
    bool p_fill_1d = true;
    bool p_fill_2d = true;
    bool p_fixed = false;
    double p_fixedvalue = 0;

    if (partree->GetBranch("latex_name"))
        partree->SetBranchAddress("latex_name", p_latexname);
    if (partree->GetBranch("precision"))
        partree->SetBranchAddress("precision", &p_precision);
    if (partree->GetBranch("nbins"))
        partree->SetBranchAddress("nbins", &p_nbins);
    if (partree->GetBranch("fill_1d"))
        partree->SetBranchAddress("fill_1d", &p_fill_1d);
    if (partree->GetBranch("fill_2d"))
        partree->SetBranchAddress("fill_2d", &p_fill_2d);
    if (partree->GetBranch("fixed"))
        partree->SetBranchAddress("fixed", &p_fixed);
    if (partree->GetBranch("fixed_value"))
        partree->SetBranchAddress("fixed_value", &p_fixedvalue);

    partree->BuildIndex("parameter", "index");

    // load parameters
    // loop over entries with "parameter" branch = true
    for (unsigned i = 0; partree->GetEntryNumberWithIndex(1, i) >= 0; ++i) {
        partree->GetEntryWithIndex(1, i);
        BCParameter Par(p_name, p_lowerlimit, p_upperlimit, p_latexname);
        if (p_fixed)
            Par.Fix(p_fixedvalue);
        Par.SetPrecision(p_precision);
        Par.FillHistograms(p_fill_1d, p_fill_2d);
        Par.SetNbins(p_nbins);
        AddParameter(Par);
    }

    // load user-defined observables
    if (loadObservables) {
        fObservables = BCObservableSet();
        // loop over entries with "parameter" branch = false
        for (unsigned i = 0; partree->GetEntryNumberWithIndex(0, i) >= 0; ++i) {
            partree->GetEntryWithIndex(0, i);
            BCObservable Obs(p_name, p_lowerlimit, p_upperlimit, p_latexname);
            Obs.SetPrecision(p_precision);
            Obs.FillHistograms(p_fill_1d, p_fill_2d);
            Obs.SetNbins(p_nbins);
            AddObservable(Obs);
        }
    }

    partree->ResetBranchAddresses();
}

// --------------------------------------------------------
void BCEngineMCMC::LoadMCMCParameters(TTree& partree)
{
    if (partree.GetEntries() < 1)
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree is empty");
    if (!partree.GetBranch("parameter"))
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree missing parameter branch");
    if (!partree.GetBranch("index"))
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree missing index branch");
    if (!partree.GetBranch("nchains"))
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree missing nchains branch");
    if (!partree.GetBranch("efficiency_0"))
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree missing efficiency_0 branch");
    if (!partree.GetBranch("scale"))
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: tree missing scale branch");

    partree.ResetBranchAddresses();

    // set number of chains
    unsigned p_nchains;
    partree.SetBranchAddress("nchains", &p_nchains);
    partree.GetEntry(0);
    if (p_nchains == 0)
        throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: no chains in branch");
    SetNChains(p_nchains);

    std::vector<double> p_efficiency(GetNChains(), -1);
    std::vector<double> p_scale(GetNChains(), -1);

    partree.SetBranchAddress("efficiency_0", &p_efficiency[0]);
    partree.SetBranchAddress("scale", &p_scale[0]);

    std::vector<std::vector<double> > scales(GetNChains(), std::vector<double>(GetNParameters(), -1));
    std::vector<std::vector<double> > efficiencies(GetNChains(), std::vector<double>(GetNParameters(), -1));

    partree.BuildIndex("parameter", "index");

    for (unsigned p = 0; p < GetNParameters(); ++p) {
        partree.GetEntryWithIndex(1, p);
        for (unsigned c = 0; c < GetNChains(); ++c) {
            scales[c][p] = p_scale[c];
            efficiencies[c][p] = p_efficiency[c];
        }
    }

    for (unsigned c = 0; c < GetNChains(); ++c)
        for (unsigned p = 0; p < GetNParameters(); ++p)
            if (scales[c][p] < 0 || efficiencies[c][p] == -1)
                throw std::runtime_error("BCEngineMCMC::LoadMCMCParameters: unset scale or efficiency.");

    fMCMCProposalFunctionScaleFactor = scales;

    // if multivariate proposal, then all parameters have same efficiency in each chain
    // (even fixed ones)
    if (GetNParameters() > 1) {
        SetProposeMultivariate(true);
        for (unsigned c = 0; c < GetNChains() and GetProposeMultivariate(); ++c)
            for (unsigned p = 1; p < GetNParameters() and GetProposeMultivariate(); ++p)
                if (efficiencies[c][p] != efficiencies[c][p - 1])
                    SetProposeMultivariate(false);
    }
    // if only one parameter, leave to what user has set.

    partree.ResetBranchAddresses();
}

// --------------------------------------------------------
bool BCEngineMCMC::ParameterTreeMatchesModel(TTree* partree, bool checkObservables)
{
    bool p_fixed;
    char p_name[200];
    double p_lowerlimit, p_upperlimit, p_fixedvalue;
    partree->SetBranchAddress("name", p_name);
    partree->SetBranchAddress("lower_limit", &p_lowerlimit);
    partree->SetBranchAddress("upper_limit", &p_upperlimit);
    partree->BuildIndex("parameter", "index");

    bool has_fixed = true;
    if (partree->GetBranch("fixed"))
        partree->SetBranchAddress("fixed", &p_fixed);
    else
        has_fixed = false;

    bool has_fixed_value = true;
    if (partree->GetBranch("fixed_value"))
        partree->SetBranchAddress("fixed_value", &p_fixedvalue);
    else
        has_fixed_value = false;

    // check parameters
    for (unsigned i = 0; i < GetNParameters(); ++i) {
        if (partree->GetEntryNumberWithIndex(1, i) < 0) {
            BCLog::OutError("BCEngineMCMC::ParameterTreeMatchesModel : Parameter tree contains too few entries.");
            return false;
        }
        partree->GetEntryWithIndex(1, i);
        if (!GetParameter(i).IsNamed(p_name)) {
            BCLog::OutError(Form("BCEngineMCMC::ParameterTreeMatchesModel : Parameter[%d]'s names do not match.", i));
            return false;
        }
        if (GetParameter(i).GetLowerLimit() != p_lowerlimit)
            BCLog::OutDetail(Form("BCEngineMCMC::ParameterTreeMatchesModel : Lower limit of parameter \"%s\" does not match.", p_name));
        if (GetParameter(i).GetUpperLimit() != p_upperlimit)
            BCLog::OutDetail(Form("BCEngineMCMC::ParameterTreeMatchesModel : Upper limit of parameter \"%s\" does not match.", p_name));
        if (has_fixed && GetParameter(i).Fixed() != p_fixed) {
            if (p_fixed) {
                BCLog::OutDetail(Form("BCEngineMCMC::ParameterTreeMatchesModel : Fixed status of parameter \"%s\" does not match. Fixing it to %f.", p_name, p_fixedvalue));
                GetParameter(i).Fix(p_fixedvalue);
            } else {
                BCLog::OutDetail(Form("BCEngineMCMC::ParameterTreeMatchesModel : Fixed status of parameter \"%s\" does not match. Unfixing it.", p_name));
                GetParameter(i).Unfix();
            }
        }
        if (has_fixed && GetParameter(i).Fixed() && has_fixed_value && GetParameter(i).GetFixedValue() != p_fixedvalue) {
            BCLog::OutDetail(Form("BCEngineMCMC::ParameterTreeMatchesModel : Fixed value of parameter \"%s\" does not match. Updating it.", p_name));
            GetParameter(i).Fix(p_fixedvalue);
        }
    }
    if (!checkObservables)
        return true;
    // check observables
    for (unsigned i = 0; i < GetNObservables(); ++i) {
        if (partree->GetEntryNumberWithIndex(0, i) < 0) {
            BCLog::OutError("BCEngineMCMC::ParameterTreeMatchesModel : Parameters tree contains too few entries.");
            return false;
        }
        partree->GetEntryWithIndex(0, i);
        if (!GetObservable(i).IsNamed(p_name)) {
            BCLog::OutError(Form("BCEngineMCMC::ParameterTreeMatchesModel : Observable[%d]'s names do not match.", i));
            return false;
        }
        if (GetObservable(i).GetLowerLimit() != p_lowerlimit)
            BCLog::OutWarning(Form("BCEngineMCMC::ParameterTreeMatchesModel : Lower limit of observable \"%s\" does not match.", p_name));
        if (GetObservable(i).GetUpperLimit() != p_upperlimit)
            BCLog::OutWarning(Form("BCEngineMCMC::ParameterTreeMatchesModel : Upper limit of observable \"%s\" does not match.", p_name));
    }
    return true;
}

// --------------------------------------------------------
void BCEngineMCMC::LoadMCMC(const std::string& filename, std::string mcmcTreeName, std::string parameterTreeName, bool loadObservables)
{
    // save current directory
    TDirectory* dir = gDirectory;

    TFile* inputfile = TFile::Open(filename.data(), "READ");
    if (!inputfile || inputfile->IsZombie()) {
        gDirectory = dir;
        throw std::runtime_error(Form("BCEngineMCMC::LoadMCMC: Could not open file %s.", filename.data()));
    }

    if (mcmcTreeName.empty() && parameterTreeName.empty()) {
        // look through file for trees named according to BAT scheme
        TList* LoK = inputfile->GetListOfKeys();
        std::vector<std::string> mcmc_names;
        std::vector<std::string> parameter_names;
        for (int i = 0; i < LoK->GetEntries(); ++i) {
            TKey* k = (TKey*)(LoK->At(i));
            if (strcmp(k->GetClassName(), "TTree") != 0)
                continue;
            std::string treeName(k->GetName());
            if (treeName.find_last_of("_") == std::string::npos)
                continue;
            if (treeName.substr(treeName.find_last_of("_")) == "_mcmc")
                mcmc_names.push_back(treeName.substr(0, treeName.find_last_of("_")));
            else if (treeName.substr(treeName.find_last_of("_")) == "_parameters")
                parameter_names.push_back(treeName.substr(0, treeName.find_last_of("_")));
        }

        std::vector<std::string> model_names;
        for (unsigned i = 0; i < mcmc_names.size(); ++i)
            for (unsigned j = 0; j < parameter_names.size(); ++j)
                if (mcmc_names[i] == parameter_names[j])
                    model_names.push_back(mcmc_names[i]);

        if (model_names.empty())
            throw std::runtime_error(Form("BCEngineMCMC::LoadMCMC : %s contains no matching MCMC and Parameter trees.", filename.data()));

        if (model_names.size() > 1) {
            std::string model_names_string = model_names[0];
            for (unsigned i = 0; i < model_names.size(); ++i)
                model_names_string += ", " + model_names[i];
            throw std::runtime_error(Form("BCEngineMCMC::LoadMCMC : %s contains more than one model, please select one by providing a model name: %s", filename.data(), model_names_string.data()));
        }

        mcmcTreeName = model_names[0] + "_mcmc";
        parameterTreeName = model_names[0] + "_parameters";

        if (fName.empty())
            SetName(model_names[0]);
    } else if (fName.empty()) {
        // check mcmcTreeName and parameterTreeName for default BAT name scheme [modelname]_mcmc/parameters:
        if (mcmcTreeName.find_last_of("_") != std::string::npos && mcmcTreeName.substr(mcmcTreeName.find_last_of("_")) == "_mcmc"
                && parameterTreeName.find_last_of("_") != std::string::npos && parameterTreeName.substr(parameterTreeName.find_last_of("_")) == "_parameters") {
            fName = mcmcTreeName.substr(0, mcmcTreeName.find_last_of("_"));
        }
    }

    // set tree names if empty
    if (mcmcTreeName.empty())
        mcmcTreeName = Form("%s_mcmc", GetSafeName().data());
    if (parameterTreeName.empty()) // default parameter tree name
        parameterTreeName = Form("%s_parameters", GetSafeName().data());

    TTree* mcmcTree = NULL;
    inputfile->GetObject(mcmcTreeName.data(), mcmcTree);
    if (!mcmcTree)
        throw std::runtime_error(Form("BCEngineMCMC::LoadMCMC : %s does not contain a tree named %s", filename.data(), mcmcTreeName.data()));


    TTree* parTree = NULL;
    inputfile->GetObject(parameterTreeName.data(), parTree);
    if (!parTree)
        throw std::runtime_error(Form("BCEngineMCMC::LoadMCMC : %s does not contain a tree named %s", filename.data(), parameterTreeName.data()));

    gDirectory = dir;
    LoadMCMC(mcmcTree, parTree, loadObservables);
}

// --------------------------------------------------------
void BCEngineMCMC::LoadMCMC(TTree* mcmcTree, TTree* parTree, bool loadObservables)
{
    fMCMCTreeLoaded = false;
    fMCMCTreeReuseObservables = loadObservables;

    if (!mcmcTree || !parTree)
        throw std::runtime_error("BCEngineMCMC::LoadMCMC : empty trees provided");

    // load parameter tree
    if (!ValidParameterTree(parTree))
        throw std::runtime_error("BCEngineMCMC::LoadMCMC : invalid parameter tree");

    delete fParameterTree;
    fParameterTree = parTree;

    // if parameters is empty
    if (GetNParameters() == 0)
        LoadParametersFromTree(fParameterTree, fMCMCTreeReuseObservables);
    // else check parameter tree
    else if (!ParameterTreeMatchesModel(fParameterTree, fMCMCTreeReuseObservables))
        throw std::runtime_error("BCEngineMCMC::LoadMCMC : Parameter tree does not match model.");

    LoadMCMCParameters(*fParameterTree);

    // check mcmc tree
    if (!ValidMCMCTree(mcmcTree, fMCMCTreeReuseObservables))
        throw std::runtime_error("BCEngineMCMC::LoadMCMC : invalid MCMC tree");

    delete fMCMCTree;
    fMCMCTree = mcmcTree;

    fMCMCTreeLoaded = true;
}

// --------------------------------------------------------
void BCEngineMCMC::Remarginalize(bool autorange)
{
    // Check if tree is valid
    if (!ValidMCMCTree(fMCMCTree, fMCMCTreeReuseObservables))
        return;

    // link chain and phase
    fMCMCTree->SetBranchAddress("Chain",          &fMCMCTree_Chain);
    fMCMCTree->SetBranchAddress("Phase",          &fMCMCPhase);

    // link log(probability) if available
    if (fMCMCTree->GetBranch("LogProbability"))
        fMCMCTree->SetBranchAddress("LogProbability", &fMCMCTree_State.log_probability);

    // link iteration if available
    if (fMCMCTree->GetBranch("Iteration"))
        fMCMCTree->SetBranchAddress("Iteration",      &fMCMCTree_State.iteration);
    else
        fMCMCTree_State.iteration = 0;

    // link log(likelihood) if available
    if (fMCMCTree->GetBranch("LogLikelihood"))
        fMCMCTree->SetBranchAddress("LogLikelihood", &fMCMCTree_State.log_likelihood);
    else
        fMCMCTree_State.log_likelihood = -std::numeric_limits<double>::infinity();

    // link log(prior) if available
    if (fMCMCTree->GetBranch("LogPrior"))
        fMCMCTree->SetBranchAddress("LogPrior", &fMCMCTree_State.log_prior);
    else
        fMCMCTree_State.log_prior = -std::numeric_limits<double>::infinity();

    // link parameters
    fMCMCTree_State.parameters.assign(GetNParameters(), 0);
    for (unsigned i = 0; i < GetNParameters(); ++i)
        fMCMCTree->SetBranchAddress(GetParameter(i).GetSafeName().data(), &fMCMCTree_State.parameters[i]);

    // link observables
    if (fMCMCTreeReuseObservables) {
        fMCMCTree_State.observables.assign(GetNObservables(), 0);
        for (unsigned i = 0; i < GetNObservables(); ++i)
            fMCMCTree->SetBranchAddress(GetObservable(i).GetSafeName().data(), &fMCMCTree_State.observables[i]);
    }

    // find out how many chains used to generate tree
    // (this presumes that chains were run in blocks of one iteration at a time!)
    unsigned nchains = 0;
    for (int n = 0; n < fMCMCTree->GetEntries(); ++n) {
        fMCMCTree->GetEntry(n);
        if (nchains > 0 && fMCMCTree_Chain == 0)
            break;
        if (fMCMCTree_Chain + 1 > nchains)
            nchains = fMCMCTree_Chain + 1;
    }
    SetNChains(nchains);

    MCMCInitialize();

    if (autorange) {
        std::vector<double> XMin;
        std::vector<double> XMax;
        if (fMCMCStatistics_AllChains.n_samples > 0) {
            XMin = fMCMCStatistics_AllChains.minimum;
            XMax = fMCMCStatistics_AllChains.maximum;
        } else {
            // find min and max
            XMin.assign(GetNVariables(), +std::numeric_limits<double>::infinity());
            XMax.assign(GetNVariables(), -std::numeric_limits<double>::infinity());
            for (int n = 0; n < fMCMCTree->GetEntries(); ++n) {
                fMCMCTree->GetEntry(n);

                if (fMCMCPhase <= 0)
                    continue;

                for (unsigned i = 0; i < fMCMCTree_State.parameters.size(); ++i) {
                    XMin[i] = std::min(XMin[i], fMCMCTree_State.parameters[i]);
                    XMax[i] = std::max(XMax[i], fMCMCTree_State.parameters[i]);
                }
                if (fMCMCTreeReuseObservables) {
                    for (unsigned i = 0; i < fMCMCTree_State.observables.size(); ++i) {
                        XMin[GetNParameters() + i] = std::min(XMin[GetNParameters() + i], fMCMCTree_State.observables[i]);
                        XMax[GetNParameters() + i] = std::max(XMax[GetNParameters() + i], fMCMCTree_State.observables[i]);
                    }
                } else {
                    UpdateChainIndex(fMCMCTree_Chain);
                    CalculateObservables(fMCMCTree_State.parameters);
                    for (unsigned i = 0; i < GetNObservables(); ++i) {
                        XMin[GetNParameters() + i] = std::min(XMin[GetNParameters() + i], fObservables[i].Value());
                        XMax[GetNParameters() + i] = std::max(XMax[GetNParameters() + i], fObservables[i].Value());
                    }
                }
            }
        }
        if (!XMin.empty() && !XMax.empty()) {
            // store mins and maxes, and change
            std::vector<double> xmin;
            std::vector<double> xmax;
            for (unsigned i = 0; i < GetNVariables(); ++i) {
                xmin.push_back(GetVariable(i).GetLowerLimit());
                xmax.push_back(GetVariable(i).GetUpperLimit());
                GetVariable(i).SetLimits(XMin[i], XMax[i]);
            }
            CreateHistograms();
            // restore mins and maxes
            for (unsigned i = 0; i < GetNVariables(); ++i)
                GetVariable(i).SetLimits(xmin[i], xmax[i]);
        }
    }

    fMCMCStatistics.assign(fMCMCNChains, BCEngineMCMC::Statistics(GetNParameters(), GetNObservables()));
    fMCMCStatistics_AllChains.Init(GetNParameters(), GetNObservables());

    fMCMCTree_State.log_probability = -std::numeric_limits<double>::infinity();

    for (unsigned n = 0; n < fMCMCTree->GetEntries(); ++n) {
        fMCMCTree->GetEntry(n);

        fMCMCStates[fMCMCTree_Chain] = fMCMCTree_State;

        fMCMCCurrentIteration = fMCMCStates[fMCMCTree_Chain].iteration;

        // calculate observables if requested
        if (!fMCMCTreeReuseObservables) {
            CalculateObservables(fMCMCStates[fMCMCTree_Chain].parameters);
            for (unsigned i = 0; i < GetNObservables(); ++i)
                fMCMCStates[fMCMCTree_Chain].observables[i] = GetObservable(i).Value();
        }

        // store iteration of convergance
        if (fMCMCNIterationsConvergenceGlobal < 0 && fMCMCPhase > 0) {
            for (unsigned c = 0; c < fMCMCNChains; ++c)
                fMCMCStatistics[c].Reset(false, true);
            fMCMCNIterationsConvergenceGlobal = fMCMCCurrentIteration;
        }

        fMCMCStatistics[fMCMCTree_Chain].Update(fMCMCStates[fMCMCTree_Chain]);

        if (fMCMCPhase <= 0)
            continue;

        MCMCCurrentPointInterface(fMCMCStates[fMCMCTree_Chain].parameters, fMCMCTree_Chain, true);

        // This should be changed:
        // MCMCUserIterationInterface() should be deleted
        // and the whole block should be replaced with
        // InChainFillHistograms(fMCMCStates[fMCMCTree_Chain])
        if (fMCMCTree_Chain == fMCMCNChains - 1) {
            MCMCUserIterationInterface();
            InChainFillHistograms();
        }
    }

    // set mult. var. proposal function covariance to those of main run if empty
    if (fMCMCProposeMultivariate && fMultivariateProposalFunctionCovariance.empty()) {
        fMultivariateProposalFunctionCovariance.assign(fMCMCNChains, TMatrixDSym(GetNFreeParameters()));
        fMultivariateCovarianceUpdates = 0;
        UpdateMultivariateProposalFunctionCovariances(1.);
    }

    // combine statistics
    for (unsigned c = 0; c < fMCMCNChains; ++c)
        fMCMCStatistics_AllChains += fMCMCStatistics[c];

}

// --------------------------------------------------------
void BCEngineMCMC::PrepareToContinueMarginalization(const std::string& filename, const std::string& mcmcTreeName, const std::string& parameterTreeName, bool loadObservables, bool autorange)
{
    LoadMCMC(filename, mcmcTreeName, parameterTreeName, loadObservables);
    Remarginalize(autorange);
    // delete trees
    delete fMCMCTree;
    fMCMCTree = NULL;
    delete fParameterTree;
    fParameterTree = NULL;
}

// --------------------------------------------------------
bool BCEngineMCMC::UpdateMultivariateProposalFunctionCovariances(double a)
{
    if (fMultivariateProposalFunctionCovariance.size() != fMCMCNChains)
        return false;

    // Update covariance matricies
    unsigned I = 0;
    for (unsigned i = 0; i < GetNParameters(); ++i) {
        if (GetParameter(i).Fixed())
            continue;
        unsigned J = I;
        for (unsigned j = i; j < GetNParameters(); ++j) {
            if (GetParameter(j).Fixed())
                continue;
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                fMultivariateProposalFunctionCovariance[c][I][J] *= (1 - a);
                fMultivariateProposalFunctionCovariance[c][I][J] += a * fMCMCStatistics[c].covariance[i][j];
                fMultivariateProposalFunctionCovariance[c][J][I] = fMultivariateProposalFunctionCovariance[c][I][J];
            }
            ++J;
        }
        ++I;
    }

    return true;
}

// --------------------------------------------------------
bool BCEngineMCMC::UpdateMultivariateProposalFunctionCovariances()
{
    // a = (1+t)^(-lambda)
    double a = pow(1. + fMultivariateCovarianceUpdates + 1, -fMultivariateCovarianceUpdateLambda);

    if (!UpdateMultivariateProposalFunctionCovariances(a))
        return false;

    ++fMultivariateCovarianceUpdates;
    return true;
}
// --------------------------------------------------------
void BCEngineMCMC::CalculateCholeskyDecompositions()
{
    if (fMultivariateProposalFunctionCovariance.size() != fMCMCNChains)
        throw std::runtime_error("BCEngineMCMC::CalculateCholeskyDecompositions: size of fMultivariateProposalFunctionCovariance does not match number of chains.");

    // clear decomps
    fMultivariateProposalFunctionCholeskyDecomposition.assign(fMCMCNChains, TMatrixD(GetNFreeParameters(), GetNFreeParameters()));

    // create decomposer
    TDecompChol CholeskyDecomposer;
    // Update cholesky decompositions
    for (unsigned c = 0; c < fMCMCNChains; ++c) {

        // try cholesky decomposition
        CholeskyDecomposer.SetMatrix(fMultivariateProposalFunctionCovariance[c]*fMCMCProposalFunctionScaleFactor[c][0]);
        if (CholeskyDecomposer.Decompose())
            fMultivariateProposalFunctionCholeskyDecomposition[c].Transpose(CholeskyDecomposer.GetU());

        else {
            // try with covariance + epsilon
            BCLog::OutDetail(Form("BCEngineMCMC::CalculateCholeskyDecompositions : chain %u Cholesky decomposition failed! Adding epsilon*I and trying again.", c));
            TMatrixDSym U(fMultivariateProposalFunctionCovariance[c]*fMCMCProposalFunctionScaleFactor[c][0]);
            for (int i = 0; i < U.GetNrows(); ++i)
                U[i][i] *= (1 + fMultivariateEpsilon);
            CholeskyDecomposer.SetMatrix(U);
            if (CholeskyDecomposer.Decompose())
                fMultivariateProposalFunctionCholeskyDecomposition[c].Transpose(CholeskyDecomposer.GetU());

            else {
                // diagonalize
                BCLog::OutDetail(Form("BCEngineMCMC::CalculateCholeskyDecompositions : chain %u Cholesky decomposition failed! Setting off-diagonal elements of covariance to zero", c));
                TMatrixDSym U(fMultivariateProposalFunctionCovariance[c]*fMCMCProposalFunctionScaleFactor[c][0]);
                for (int i = 0; i < fMultivariateProposalFunctionCholeskyDecomposition[c].GetNrows(); ++i)
                    for (int j = 0; j < fMultivariateProposalFunctionCholeskyDecomposition[c].GetNcols(); ++j)
                        if (i != j)
                            U[i][j] = 0;
                CholeskyDecomposer.SetMatrix(U);
                if (CholeskyDecomposer.Decompose())
                    fMultivariateProposalFunctionCholeskyDecomposition[c].Transpose(CholeskyDecomposer.GetU());
                else {
                    throw std::runtime_error(Form("BCEngineMCMC::CalculateCholeskyDecompositions : chain %u Cholesky decomposition failed! No remedies!", c));
                }
            }
        }
    }
}

// --------------------------------------------------------
double BCEngineMCMC::ProposalFunction(unsigned ichain, unsigned iparameter)
{
    // multiply by 1 (dof <=0, Gauss) or a random variate that scales the Gaussian to a Student's t with dof degrees of freedom
    const double scale = fMCMCThreadLocalStorage[ichain].scale(fMCMCProposalFunctionDof);

    return scale * fMCMCThreadLocalStorage[ichain].rng->Gaus(0, fMCMCProposalFunctionScaleFactor[ichain][iparameter]);
}

// --------------------------------------------------------
bool BCEngineMCMC::GetProposalPointMetropolis(unsigned chain, std::vector<double>& x)
{
    // copy the current point into the new
    x = fMCMCStates[chain].parameters;

    // generate N-Free N(0,1) random values
    TVectorD& y = fMCMCThreadLocalStorage[chain].yLocal;
    for (int i = 0; i < y.GetNrows(); ++i)
        y[i] = fMCMCThreadLocalStorage[chain].rng->Gaus(0, 1);

    // multiply by Cholesky decomposition
    y *= fMultivariateProposalFunctionCholeskyDecomposition[chain];

    // multiply by 1 (dof <=0, Gauss) or a random variate that scales the Gaussian to a Student's t with dof degrees of freedom
    const double scale = fMCMCThreadLocalStorage[chain].scale(fMCMCProposalFunctionDof);

    // add values into x
    int I = 0;
    for (unsigned i = 0; i < GetNParameters() && I < y.GetNrows(); ++i)
        if (!GetParameter(i).Fixed()) {
            x[i] += y[I] * scale;
            ++I;
        }

    // return whether point is within limits, ignoring fixed parameters
    return GetParameters().IsWithinLimits(x);
}

// --------------------------------------------------------
bool BCEngineMCMC::GetProposalPointMetropolis(unsigned ichain, unsigned ipar, std::vector<double>& x)
{
    // copy the current point into the new
    x = fMCMCStates[ichain].parameters;

    // check if parameter is fixed
    if (GetParameter(ipar).Fixed()) {
        x[ipar] = GetParameter(ipar).GetFixedValue();
        return true; // assume that value is inside allowed region
    }

    // get unscaled random point in the dimension of the chosen
    // parameter. this point might not be in the correct volume.
    double proposal = ProposalFunction(ichain, ipar);

    // modify the parameter under study
    x[ipar] += proposal * GetParameter(ipar).GetRangeWidth();

    // check if the point is in the correct volume.
    return GetParameter(ipar).IsWithinLimits(x[ipar]);
}

// --------------------------------------------------------
bool BCEngineMCMC::AcceptOrRejectPoint(unsigned chain, unsigned parameter)
{
    // retrieve current probability
    double p0 = fMCMCStates[chain].log_probability;
    if (!std::isfinite(p0)) p0 = -std::numeric_limits<double>::max();
    // calculate proposed probability
    const double p1 = LogEval(fMCMCThreadLocalStorage[chain].parameters);

    // if the new point is more probable, keep it; or else throw dice
    if (std::isfinite(p1) && (p1 >= p0 || log(fMCMCThreadLocalStorage[chain].rng->Rndm()) < (p1 - p0))) {
        // accept point
        fMCMCStates[chain] = fMCMCThreadLocalStorage[chain];
        // increase efficiency
        fMCMCStatistics[chain].efficiency[parameter] += (1. - fMCMCStatistics[chain].efficiency[parameter]) / (fMCMCStatistics[chain].n_samples_efficiency + 1.);
        // execute user code and return
        MCMCCurrentPointInterface(fMCMCStates[chain].parameters, chain, true);
        return true;
    }

    // else decrease efficiency
    fMCMCStatistics[chain].efficiency[parameter] *= 1.*fMCMCStatistics[chain].n_samples_efficiency / (fMCMCStatistics[chain].n_samples_efficiency + 1.);

    // if log(likelihood) of proposed point was not a finite number
    if (!std::isfinite(p1)) {
        if (fMCMCProposeMultivariate) {
            BCLog::OutDebug(Form("log(probability) evaluated to nan or inf in chain %i while at ", chain));
            PrintParameters(fMCMCThreadLocalStorage[chain].parameters, BCLog::OutDebug);
        } else
            BCLog::OutDebug(Form("log(probability) evaluated to nan or inf in chain %i while varying parameter %s to %.3e", chain, GetParameter(parameter).GetName().data(), fMCMCThreadLocalStorage[chain].parameters[parameter]));
    }

    // execute user code and return
    MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].parameters, chain, false);
    return false;
}

// --------------------------------------------------------
bool BCEngineMCMC::GetNewPointMetropolis(unsigned chain, unsigned parameter)
{
    // increase counter
    ++fMCMCStates[chain].iteration;

    // get proposal point
    if (GetProposalPointMetropolis(chain, parameter, fMCMCThreadLocalStorage[chain].parameters))
        return AcceptOrRejectPoint(chain, parameter);

    // execute user code and return
    MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].parameters, chain, false);
    return false;
}

// --------------------------------------------------------
bool BCEngineMCMC::GetNewPointMetropolis(unsigned chain)
{
    // increase counter
    ++fMCMCStates[chain].iteration;

    // get proposal point
    if (GetProposalPointMetropolis(chain, fMCMCThreadLocalStorage[chain].parameters))
        return AcceptOrRejectPoint(chain, 0);

    // execute user code and return
    MCMCCurrentPointInterface(fMCMCThreadLocalStorage[chain].parameters, chain, false);
    return false;
}

//--------------------------------------------------------
bool BCEngineMCMC::GetNewPointMetropolis()
{
    bool return_value = true;

    if (!fMCMCProposeMultivariate) {
        /* run over parameters one at a time */

        for (unsigned ipar = 0; ipar < GetNParameters(); ++ipar) {
            if (GetParameter(ipar).Fixed())
                continue;

            //loop over chains
            #pragma omp parallel for schedule(static)
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                UpdateChainIndex(c);
                return_value &= GetNewPointMetropolis(c, ipar);
            }
        }

    } else {
        /* run over all pars at once */

        //loop over chains
        #pragma omp parallel for schedule(static)
        for (unsigned c = 0; c < fMCMCNChains; ++c) {
            UpdateChainIndex(c);
            return_value &= GetNewPointMetropolis(c);
        }
    }

    // increase number of iterations used in each chain for calculating efficiencies
    for (unsigned c = 0; c < fMCMCNChains; ++c)
        fMCMCStatistics[c].n_samples_efficiency += 1;

    ++fMCMCCurrentIteration;
    return return_value;
}

// --------------------------------------------------------
void BCEngineMCMC::InChainFillHistograms(const ChainState& cs)
{
    ////////////////////////////////////////
    // fill each 1-dimensional histogram that exists
    for (unsigned i = 0; i < GetNVariables() && i < fH1Marginalized.size(); ++i)
        if (TH1* h = fH1Marginalized[i]) {
            if (i < GetNParameters())
                h->Fill(cs.parameters[i]);
            else if (i - GetNParameters() < GetNObservables())
                h->Fill(cs.observables[i - GetNParameters()]);
        }

    ////////////////////////////////////////
    // fill each 2-dimensional histogram that exists
    for (unsigned j = 0; j < GetNVariables() && j < fH2Marginalized.size(); ++j)
        for (unsigned k = 0; k < GetNVariables() && k < fH2Marginalized[j].size(); ++k)
            if (TH2* h = fH2Marginalized[j][k]) {
                if (j < GetNParameters()) {
                    if (k < GetNParameters())
                        h->Fill(cs.parameters[j], cs.parameters[k]);
                    else if (k - GetNParameters() < GetNObservables())
                        h->Fill(cs.parameters[j], cs.observables[k - GetNParameters()]);
                } else if (j - GetNParameters() < GetNObservables()) {
                    if (k < GetNParameters())
                        h->Fill(cs.observables[j - GetNParameters()], cs.parameters[k]);
                    else if (k - GetNParameters() < GetNObservables())
                        h->Fill(cs.observables[j - GetNParameters()], cs.observables[k - GetNParameters()]);
                }
            }
}

// --------------------------------------------------------
void BCEngineMCMC::InChainFillHistograms()
{
    // loop over chains
    for (unsigned c = 0; c < fMCMCNChains; ++c)
        InChainFillHistograms(fMCMCStates[c]);
}

// --------------------------------------------------------
void BCEngineMCMC::InChainFillTree(const ChainState& cs, unsigned chain_number)
{
    if (!fMCMCTree)
        return;
    fMCMCTree_Chain = chain_number;
    fMCMCTree_State = cs;
    fMCMCTree->Fill();
}

// -------------------------------------------------------
void BCEngineMCMC::InChainFillTree()
{
    // loop over all chains
    for (unsigned i = 0; i < fMCMCNChains; ++i)
        InChainFillTree(fMCMCStates[i], i);
}

//---------------------------------------------------------
void BCEngineMCMC::CloseOutputFile()
{
    if (!fMCMCOutputFile)
        return;

    if (fMCMCOutputFile->IsOpen() && fMCMCOutputFile->IsWritable()) {
        fMCMCOutputFile->Write(0, TObject::kWriteDelete);
        fMCMCOutputFile->Close();
    }
    // ROOT also deletes associated named objects, the trees
    delete fMCMCOutputFile;
    fMCMCOutputFile = NULL;
    fMCMCTree = NULL;
    fParameterTree = NULL;
    fMCMCTreeLoaded = false;
}

// --------------------------------------------------------
bool BCEngineMCMC::MetropolisPreRun()
{
    // print on screen
    BCLog::OutSummary(Form("Pre-run Metropolis MCMC for model \"%s\" ...", GetName().data()));

    // initialize Markov chain
    MCMCInitialize();

    if (fMCMCFlagWritePreRunToFile)
        InitializeMarkovChainTree();

    // perform run
    BCLog::OutSummary(Form(" --> Perform MCMC pre-run with %i chains, each with maximum %i iterations", fMCMCNChains, fMCMCNIterationsPreRunMax));

    const int old_error_ignore_level = gErrorIgnoreLevel;

    fMultivariateProposalFunctionCovariance.clear();
    fMultivariateProposalFunctionCholeskyDecomposition.clear();

    if (fMCMCProposeMultivariate) {
        // multivariate proposal function

        // suppress ROOT errors about Cholesky decomposition
        gErrorIgnoreLevel = kBreak;

        // initialize covariance matrices as diag(var_0, var_1, var_2, ...)
        // initialize cholesky decomposition as diag(std_0, std_1, std_2, ...)
        // for free parameters only
        TMatrixDSym S0(GetNFreeParameters());
        TMatrixD    CD0(GetNFreeParameters(), GetNFreeParameters());
        unsigned I = 0;
        for (unsigned i = 0; i < GetNParameters(); ++i)
            if (!(GetParameter(i).Fixed())) {
                if (GetParameter(i).GetPrior() != NULL) {
                    S0[I][I] = GetParameter(i).GetPriorVariance();
                    if (!std::isfinite(S0[I][I]))
                        S0[I][I] = GetParameter(i).GetRangeWidth() * GetParameter(i).GetRangeWidth() / 12;
                } else
                    S0[I][I] = GetParameter(i).GetRangeWidth() * GetParameter(i).GetRangeWidth() / 12;
                CD0[I][I] = sqrt(fMCMCProposalFunctionScaleFactor[0][0] * S0[I][I]);
                ++I;
            }

        // if chains run a second time on a model, the assign does not update the contents. wtf?
        fMultivariateProposalFunctionCovariance.assign(fMCMCNChains, S0);
        fMultivariateProposalFunctionCholeskyDecomposition.assign(fMCMCNChains, CD0);
    }

    //////////////////////////////////////////////////
    // Adjust scales until all parameters are in correct efficiency range in all chains
    bool allEfficient = false;
    bool inefficientScalesAdjustable = true;
    fMCMCCurrentIteration = 0;
    fMCMCPhase = BCEngineMCMC::kPreRun;

    unsigned nIterationsPreRunCheck = fMCMCNIterationsPreRunCheck;

    // autosave every 10 checks
    if (fMCMCTree) {
        fMCMCTree->AutoSave("SaveSelf");
    }

    fMCMCNIterationsConvergenceGlobal = -1;

    // Cholesky Decomposer for multivariate proposal function
    TDecompChol CholeskyDecomposer;

    unsigned nChecks = 0;

    // While loop criteria---do while:
    //     not yet at maximum number of iterations
    // AND (    not yet above minimum number of iterations
    //       OR an efficiency is out of range, and still adjustable;
    //       OR the chains have not converged (if using more than one chain);
    //       OR the minimum number of tuning steps have been made to the multivariate proposal function, if using it.
    //     )
    while (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax
            && (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin
                || (!allEfficient && inefficientScalesAdjustable)
                || (fMCMCNChains > 1 && fMCMCNIterationsConvergenceGlobal < 0))) {

        // Generate (nIterationsCheckConvergence) new points in each chain
        for (unsigned i = 0; i < nIterationsPreRunCheck && fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMax; ++i) {
            // get new point & calculate observables
            GetNewPointMetropolis();
            EvaluateObservables();

            // update chain statistics
            for (unsigned c = 0; c < fMCMCNChains; ++c)
                fMCMCStatistics[c].Update(fMCMCStates[c]);

            // update output tree
            if (fMCMCFlagWritePreRunToFile)
                InChainFillTree();
        }

        //////////////////////////////////////////
        // Adjust scales until efficiencies within range
        allEfficient = true;
        inefficientScalesAdjustable = false;

        bool scalesAdjusted = false;

        for (unsigned c = 0; c < fMCMCNChains; ++c) {

            if (fMCMCProposeMultivariate) { // multivariate proposal function, one efficiency per chain

                if (fMCMCStatistics[c].efficiency[0] >= fMCMCEfficiencyMin && fMCMCStatistics[c].efficiency[0] <= fMCMCEfficiencyMax)
                    continue; // since chain efficiency is in range,

                if (allEfficient) // print header if encountering first bad efficiency
                    BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within predefined range after %i iterations. Efficiency of ", fMCMCCurrentIteration));
                allEfficient = false;

                double oldScale = fMCMCProposalFunctionScaleFactor[c][0];

                if (fMCMCStatistics[c].efficiency[0] < fMCMCEfficiencyMin) { // efficiency too low... decrease scale factor
                    fMCMCProposalFunctionScaleFactor[c][0] /= fMultivariateScaleMultiplier;

                    if (fMCMCProposalFunctionScaleFactor[c][0] > fMCMCScaleFactorLowerLimit) { // still room to tune
                        if (fMCMCStatistics[c].efficiency[0] == 0)
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% (zero). Scale decreased to %.4g", c, 100 * fMCMCEfficiencyMin, fMCMCProposalFunctionScaleFactor[c][0]));
                        else if (fMCMCStatistics[c].efficiency[0] < 1.e-2)
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% (%.0g %%). Scale decreased to %.4g", c, 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[0], fMCMCProposalFunctionScaleFactor[c][0]));
                        else
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% (%.0f %%). Scale decreased to %.4g", c, 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[0], fMCMCProposalFunctionScaleFactor[c][0]));
                        inefficientScalesAdjustable = true;

                    } else { // no more room to tune
                        fMCMCProposalFunctionScaleFactor[c][0] = fMCMCScaleFactorLowerLimit;
                        if (fMCMCStatistics[c].efficiency[0] == 0)
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% (zero). Scale now at lower limit (%.4g)", c, 100 * fMCMCEfficiencyMin, fMCMCScaleFactorLowerLimit));
                        else if (fMCMCStatistics[c].efficiency[0] < 1.e-2)
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% (%.0g %%). Scale now at lower limit (%.4g)", c, 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[0], fMCMCScaleFactorLowerLimit));
                        else
                            BCLog::OutDetail(Form("         chain %d is below %.0f %% f(%.0f %%). Scale now at lower limit (%.4g)", c, 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[0], fMCMCScaleFactorLowerLimit));

                    }

                } else { // efficiency too high... increase scale factor
                    fMCMCProposalFunctionScaleFactor[c][0] *= fMultivariateScaleMultiplier;

                    if (fMCMCProposalFunctionScaleFactor[c][0] < fMCMCScaleFactorUpperLimit) { // still room to tune
                        BCLog::OutDetail(Form("         chain %d is above %.0f %% (%.0f %%). Scale increased to %.4g", c, 100 * fMCMCEfficiencyMax, 100 * fMCMCStatistics[c].efficiency[0], fMCMCProposalFunctionScaleFactor[c][0]));
                        inefficientScalesAdjustable = true;
                    } else { // no more room to tune
                        fMCMCProposalFunctionScaleFactor[c][0] = fMCMCScaleFactorUpperLimit;
                        BCLog::OutDetail(Form("         chain %d is above %.0f %% (%.0f %%). Scale now at upper limit (%.4g)", c, 100 * fMCMCEfficiencyMax, 100 * fMCMCStatistics[c].efficiency[0], fMCMCScaleFactorUpperLimit));
                    }
                }

                if (oldScale != fMCMCProposalFunctionScaleFactor[c][0])
                    scalesAdjusted = true;

            } else { // factorized proposal function, one efficiency per parameter per chain

                for (unsigned p = 0; p < GetNParameters(); ++p) {

                    if (GetParameter(p).Fixed())
                        continue;

                    if (fMCMCStatistics[c].efficiency[p] >= fMCMCEfficiencyMin && fMCMCStatistics[c].efficiency[p] <= fMCMCEfficiencyMax)
                        continue; // since parameter efficiency is in range for this chain

                    if (allEfficient) // print header if first bad efficiency
                        BCLog::OutDetail(Form("     * Efficiency status: Efficiencies not within pre-defined range after %i iterations. Efficiency of ", fMCMCCurrentIteration));
                    allEfficient = false;

                    double oldScale = fMCMCProposalFunctionScaleFactor[c][p];

                    if (fMCMCStatistics[c].efficiency[p] < fMCMCEfficiencyMin) { // efficiency too low... decrease scale factor
                        fMCMCProposalFunctionScaleFactor[c][p] /= (fMCMCStatistics[c].efficiency[p] < fMCMCEfficiencyMin / 2) ? 4 : 2;

                        if (fMCMCProposalFunctionScaleFactor[c][p] > fMCMCScaleFactorLowerLimit) { // still room to tune
                            if (fMCMCStatistics[c].efficiency[p] == 0)
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (zero) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, c, fMCMCProposalFunctionScaleFactor[c][p]));
                            else if (fMCMCStatistics[c].efficiency[p] < 1.e-2)
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (%.0g %%) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCProposalFunctionScaleFactor[c][p]));
                            else
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (%.0f %%) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCProposalFunctionScaleFactor[c][p]));
                            inefficientScalesAdjustable = true;
                        } else { // no more room to tune
                            if (fMCMCStatistics[c].efficiency[p] == 0)
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (zero) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, c, fMCMCProposalFunctionScaleFactor[c][p]));
                            else if (fMCMCStatistics[c].efficiency[p] < 1.e-2)
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (%.0g %%) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCProposalFunctionScaleFactor[c][p]));
                            else
                                BCLog::OutDetail(Form("         %-*s is below %.0f %% (%.0f %%) in chain %i. Scale decreased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMin, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCProposalFunctionScaleFactor[c][p]));
                            fMCMCProposalFunctionScaleFactor[c][p] = fMCMCScaleFactorLowerLimit;
                        }

                    } else { // if efficiency too high ... increase scale factor
                        fMCMCProposalFunctionScaleFactor[c][p] *= 2;

                        if (fMCMCProposalFunctionScaleFactor[c][p] < fMCMCScaleFactorUpperLimit) { // still room to tune
                            BCLog::OutDetail(Form("         %-*s is above %.0f %% (%.0f %%) in chain %i. Scale increased to %.4g", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMax, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCProposalFunctionScaleFactor[c][p]));
                            inefficientScalesAdjustable = true;
                        } else { // no more room to tune
                            fMCMCProposalFunctionScaleFactor[c][p] = fMCMCScaleFactorUpperLimit;
                            BCLog::OutDetail(Form("         %-*s is above %.0f %% (%.0f %%) in chain %i. Scale now at upper limit (%.4g)", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), 100 * fMCMCEfficiencyMax, 100 * fMCMCStatistics[c].efficiency[p], c, fMCMCScaleFactorUpperLimit));
                        }
                    }

                    if (oldScale != fMCMCProposalFunctionScaleFactor[c][p])
                        scalesAdjusted = true;
                } // end of parameter loop
            }
        } // end of chain loop


        /////////////////////////////////////////////
        // Check convergence
        if (fMCMCNChains > 1) {

            // Calculate & check R values
            fMCMCNIterationsConvergenceGlobal = fMCMCCurrentIteration;
            for (unsigned p = 0; p < GetNParameters(); ++p) {
                if (GetParameter(p).Fixed())
                    continue;
                std::vector<double> means(fMCMCNChains, 0);
                std::vector<double> variances(fMCMCNChains, 0);
                for (unsigned c = 0; c < fMCMCNChains; ++c) {
                    means[c]     = fMCMCStatistics[c].mean[p];
                    variances[c] = fMCMCStatistics[c].variance[p];
                }
                fMCMCRValueParameters[p] = RValue(means, variances, fMCMCStatistics[0].n_samples, fCorrectRValueForSamplingVariability);
                if (fMCMCRValueParameters[p] > fMCMCRValueParametersCriterion)
                    fMCMCNIterationsConvergenceGlobal = -1;
            }

            // output results if not converged
            if (fMCMCNIterationsConvergenceGlobal <= 0) {
                BCLog::OutDetail(Form("     * Convergence status: Set of %i Markov chains did not converge after %i iterations.", fMCMCNChains, fMCMCCurrentIteration));
                BCLog::OutDetail(Form("       - %-*s :  R-Value", fParameters.MaxNameLength(), "Parameter"));

                for (unsigned p = 0; p < GetNParameters(); ++p) {
                    if (GetParameter(p).Fixed())
                        continue;

                    if (fMCMCRValueParameters[p] < fMCMCRValueParametersCriterion)
                        BCLog::OutDetail(Form("         %-*s :  %.03f", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), fMCMCRValueParameters[p]));
                    else if (std::isfinite(fMCMCRValueParameters[p]))
                        BCLog::OutDetail(Form("         %-*s :  %.03f <-- Greater than threshold", fParameters.MaxNameLength(), GetParameter(p).GetName().data(), fMCMCRValueParameters[p]));
                    else
                        BCLog::OutDetail(Form("         %-*s :  error in calculation", fParameters.MaxNameLength(), GetParameter(p).GetName().data()));
                }
            } // end convergence conditional
        } // end chains>1 conditional
        if (fMCMCTree && nChecks % 10 == 0) {
            fMCMCTree->AutoSave("SaveSelf");
        }

        ++nChecks;              // increase number of checks made

        // scales have not been adjusted and convergence has been reached (or only one chain used)
        if (!scalesAdjusted && (fMCMCNChains == 1 || fMCMCNIterationsConvergenceGlobal > 0)) {
            if (fMCMCCurrentIteration < (int)fMCMCNIterationsPreRunMin)
                // still below minimum number of prerun iterations
                BCLog::OutDetail(Form("     * Running until at least %d iterations performed in prerun. Current iteration is %d", fMCMCNIterationsPreRunMin, fMCMCCurrentIteration));
            else
                continue;       // HURRAY!
        }

        // Update multivariate proposal function covariances
        if (fMCMCProposeMultivariate) {
            if (UpdateMultivariateProposalFunctionCovariances())
                CalculateCholeskyDecompositions();
        }

        if (fMCMCCurrentIteration >= (int)fMCMCNIterationsPreRunMax)
            continue;

        // reset statistics
        for (unsigned c = 0; c < fMCMCStatistics.size(); ++c)
            if (fMCMCPreRunCheckClear > 0 && nChecks % fMCMCPreRunCheckClear == 0)
                fMCMCStatistics[c].Reset(false, true); // clear means and (co)variances, preserve mode information, clear efficiency information
            else
                fMCMCStatistics[c].ResetEfficiencies(); // clear only efficiency information

    } // end prerun iteration while loop

    // restore ROOT error ignore level
    gErrorIgnoreLevel = old_error_ignore_level;

    // output results of prerun concerning convergence and scale adjustment
    if (fMCMCNIterationsConvergenceGlobal > 0) {
        if (allEfficient)
            BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations, and all scales are adjusted.", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
        else if (!inefficientScalesAdjustable)
            BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations, but could not adjust all scales (scale limits reached).", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
        else
            BCLog::OutSummary(Form(" --> Set of %i Markov chains converged within %i iterations, but could not adjust all scales (maximum number of iterations reached).", fMCMCNChains, fMCMCNIterationsConvergenceGlobal));
        if (fMCMCProposeMultivariate)
            BCLog::OutSummary(Form(" --> %i updates to multivariate proposal function's covariances were made.", fMultivariateCovarianceUpdates));
    } else if (allEfficient)
        BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, but all scales are adjusted.", fMCMCNChains, fMCMCNIterationsPreRunMax));
    else if (!inefficientScalesAdjustable)
        BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, and could not adjust all scales (scale limits reached).", fMCMCNChains, fMCMCNIterationsPreRunMax));
    else
        BCLog::OutSummary(Form(" --> Set of %i Markov chains did not converge within %i iterations, and could not adjust all scales.", fMCMCNChains, fMCMCNIterationsPreRunMax));


    // combine statistics:
    fMCMCStatistics_AllChains.Reset();
    for (unsigned c = 0; c < fMCMCNChains; ++c)
        fMCMCStatistics_AllChains += fMCMCStatistics[c];

    if (fMCMCProposeMultivariate) {
        BCLog::OutDetail(Form(" --> Scale factors and efficiencies (measured in last %u iterations):", fMCMCStatistics.front().n_samples_efficiency));
        BCLog::OutDetail("       - Chain : Scale factor    Efficiency");
        for (unsigned c = 0; c < fMCMCNChains; ++c)
            BCLog::OutDetail(Form("           %3d :       % 6.4g        %4.1f %%", c, fMCMCProposalFunctionScaleFactor[c][0], 100.*fMCMCStatistics[c].efficiency[0]));
    } else {
        BCLog::OutDetail(Form(" --> Average scale factors and efficiencies (measured in last %u iterations):", fMCMCStatistics.front().n_samples_efficiency));
        BCLog::OutDetail(Form("       - %-*s : Scale factor    Efficiency", fParameters.MaxNameLength(), "Parameter"));
        // print scale factors and efficiencies
        for (unsigned i = 0; i < GetNParameters(); ++i) {
            if (GetParameter(i).Fixed())
                continue;
            double scalefactors = 0;
            for (unsigned j = 0; j < fMCMCNChains; ++j)
                scalefactors += fMCMCProposalFunctionScaleFactor[j][i] / fMCMCNChains;
            BCLog::OutDetail(Form("         %-*s :          % 6.4g %%        %4.1f %%", fParameters.MaxNameLength(), GetParameter(i).GetName().data(), 100.*scalefactors, 100.*fMCMCStatistics_AllChains.efficiency[i]));
        }
    }

    // reset current iteration
    fMCMCCurrentIteration = -1;

    if (fMCMCFlagWritePreRunToFile) {
        // UpdateParameterTree();
        if (fMCMCTree)
            fMCMCTree->AutoSave("SaveSelf");
    }

    // rescale histograms
    if (fRescaleHistogramRangesAfterPreRun)
        CreateHistograms(true);

    // no error
    return true;
}

// --------------------------------------------------------
double BCEngineMCMC::RValue(const std::vector<double>& means, const std::vector<double>& variances, unsigned n, bool correctForSamplingVariability)
{
    // calculate R values according to Brooks & Gelman,
    // "General Methods for Monitoring Convergence of Iterative Simulations, 1998

    if (means.empty() || variances.empty() || means.size() != variances.size() || n == 0)
        return std::numeric_limits<double>::quiet_NaN();

    unsigned m = means.size();

    // means of values
    double mean_of_means     = 0;
    double mean_of_variances = 0;
    for (unsigned c = 0; c < m; ++c) {
        mean_of_means     += means[c];
        mean_of_variances += variances[c];
    }
    mean_of_means     *= 1. / m;
    mean_of_variances *= 1. / m;

    // variances of values
    double variance_of_means     = 0;
    double variance_of_variances = 0;
    for (unsigned c = 0; c < m; ++c) {
        variance_of_means     += (means[c] - mean_of_means) * (means[c] - mean_of_means);
        variance_of_variances += (variances[c] - mean_of_variances) * (variances[c] - mean_of_variances);
    }
    variance_of_means     /= m - 1.;
    variance_of_variances /= m - 1.;

    // variance of all samples:
    double full_variance = m * (n - 1.) / (m * n - 1) * mean_of_variances + n * (m - 1.) / (m * n - 1) * variance_of_means;

    if (mean_of_variances == 0) {
        BCLOG_DEBUG("mean of variances is zero!");
        if (full_variance == 0) {
            BCLOG_DEBUG("variance of all samples is also zero!");
            return 1;
        } else
            return std::numeric_limits<double>::infinity();
    }

    double rvalue = sqrt(full_variance / mean_of_variances); // variance(all samples) / mean(chain variances)

    if (!correctForSamplingVariability)
        return rvalue;

    // else correct for initial sampling variability

    double meansquare_of_means = 0;
    for (unsigned c = 0; c < m; ++c)
        meansquare_of_means += means[c] * means[c];
    meansquare_of_means *= 1. / m;

    double covariance_of_variance_with_mean       = 0;
    double covariance_of_variance_with_meansquare = 0;
    for (unsigned c = 0; c < m; ++c) {
        covariance_of_variance_with_mean       += (variances[c] - mean_of_variances) * (means[c] - mean_of_means);
        covariance_of_variance_with_meansquare += (variances[c] - mean_of_variances) * (means[c] * means[c] - meansquare_of_means);
    }
    covariance_of_variance_with_mean       /= m - 1.;
    covariance_of_variance_with_meansquare /= m - 1.;

    double N = (n - 1.) / n;
    double M = (m + 1.) / m;

    double V = N * full_variance + M * variance_of_means;

    double varV = N * N / m * variance_of_variances + 2 * M / (m - 1) * variance_of_means * variance_of_means
                  + 2 * M * N / m * (covariance_of_variance_with_meansquare - 2 * mean_of_means * covariance_of_variance_with_mean);

    double df = 2 * V * V / varV;

    return rvalue * sqrt((df + 3) / (df + 1));
}

// --------------------------------------------------------
bool BCEngineMCMC::Metropolis()
{
    // check the number of free parameters
    if (GetNFreeParameters() <= 0) {
        BCLog::OutWarning("BCEngineMCMC::MCMCMetropolis. Number of free parameters <= 0. Do not run Metropolis.");
        return false;
    }

    // check if prerun should be performed
    if (fMCMCPhase == BCEngineMCMC::kUnsetPhase) {
        if (!MetropolisPreRun())
            return false;

        // reset statistics
        for (unsigned c = 0; c < fMCMCStatistics.size(); ++c)
            fMCMCStatistics[c].Reset(false, true); // keep mode, reset efficiencies

        // reset iterations
        for (unsigned c = 0; c < fMCMCStates.size(); ++c)
            fMCMCStates[c].iteration = 0;
    }
    if (fMCMCFlagWriteChainToFile)
        InitializeMarkovChainTree(false, false);

    // check that correct objects of correct size have been created
    // these checks will fail if the user has changed the number of chains or parameters
    // or the fixing of parameters between the pre-run and the main run, or between main runs
    if (fMCMCStatistics.size() != fMCMCNChains)
        throw std::runtime_error("BCEngineMCMC::Metropolis: size of fMCMCStatistics does not match number of chains.");
    if (fMCMCStates.size() != fMCMCNChains)
        throw std::runtime_error("BCEngineMCMC::Metropolis: size of fMCMCStates does not match number of chains.");
    if (fMCMCThreadLocalStorage.size() != fMCMCNChains)
        throw std::runtime_error("BCEngineMCMC::Metropolis: size of fMCMCThreadLocalStorage does not match number of chains.");
    for (unsigned c = 0; c < fMCMCThreadLocalStorage.size(); ++c) {
        if (fMCMCThreadLocalStorage[c].parameters.size() != GetNParameters())
            throw std::runtime_error(Form("BCEngineMCMC::Metropolis: size of fMCMCThreadLocalStorage[%u].parameters does not match number of parameters.", c));
        if (!fMCMCThreadLocalStorage[c].rng)
            throw std::runtime_error(Form("BCEngineMCMC::Metropolis: fMCMCThreadLocalStorage[%u] lacks random number generator.", c));
    }
    if (fMCMCProposeMultivariate) {
        if (fMultivariateProposalFunctionCholeskyDecomposition.empty() && fMultivariateProposalFunctionCovariance.empty())
            throw std::runtime_error(Form("BCEngineMCMC::Metropolis: pre-run was run with factorized proposal function."));
        // check if cholesky decompositions must be calculated
        bool calc_cholesky = false;
        if (fMultivariateProposalFunctionCholeskyDecomposition.size() != fMCMCNChains)
            calc_cholesky = true;
        else {
            for (unsigned c = 0; c < fMultivariateProposalFunctionCholeskyDecomposition.size(); ++c)
                if (fMultivariateProposalFunctionCholeskyDecomposition[c].GetNrows() != static_cast<int>(GetNFreeParameters()) or
                        fMultivariateProposalFunctionCholeskyDecomposition[c].GetNcols() != static_cast<int>(GetNFreeParameters()))
                    calc_cholesky = true;
        }
        if (calc_cholesky)
            CalculateCholeskyDecompositions();

        if (fMultivariateProposalFunctionCholeskyDecomposition.size() != fMCMCNChains)
            throw std::runtime_error("BCEngineMCMC::Metropolis: size of fMultivariateProposalFunctionCholeskyDecomposition does not match number of chains.");
        for (unsigned c = 0; c < fMultivariateProposalFunctionCholeskyDecomposition.size(); ++c)
            if (fMultivariateProposalFunctionCholeskyDecomposition[c].GetNrows() != static_cast<int>(GetNFreeParameters()) or
                    fMultivariateProposalFunctionCholeskyDecomposition[c].GetNcols() != static_cast<int>(GetNFreeParameters()))
                throw std::runtime_error(Form("BCEngineMCMC::Metropolis: size of fMultivariateProposalFunctionCholeskyDecomposition[%u] matrix does not match number of free parameters.", c));
        for (unsigned c = 0; c < fMCMCThreadLocalStorage.size(); ++c)
            if (fMCMCThreadLocalStorage[c].yLocal.GetNrows() != static_cast<int>(GetNFreeParameters()))
                throw std::runtime_error(Form("BCEngineMCMC::Metropolis: size of fThreadLocalStorage[%u].yLocal matrix does not match number of free parameters.", c));
    } else {
        // check that pre-run was not run with multivariate proposal
        if (!fMultivariateProposalFunctionCholeskyDecomposition.empty())
            throw std::runtime_error("BCEngineMCMC::Metropolis: pre-run was run with multivariate proposal function.");
        if (fMCMCProposalFunctionScaleFactor.size() != fMCMCNChains)
            throw std::runtime_error("BCEngineMCMC::Metropolis: size of fProposalFunctionScaleFactor does not match number of chains.");
        for (unsigned c = 0; c < fMCMCProposalFunctionScaleFactor.size(); ++c) {
            if (fMCMCProposalFunctionScaleFactor[c].size() != GetNParameters())
                throw std::runtime_error(Form("BCEngineMCMC::Metropolis: size of fMCMCProposalFunctionScaleFactor[%u] does not match number of parameters.", c));
            for (unsigned p = 0; p < fMCMCProposalFunctionScaleFactor[c].size(); ++p) {
                if (!GetParameter(p).Fixed() && fMCMCProposalFunctionScaleFactor[c][p] <= 0)
                    throw std::runtime_error(Form("BCEngineMCMC::Metropolis: fMCMCProposalFunctionScaleFactor[%u][%u] <= 0.", c, p));
                if (GetParameter(p).Fixed() && fMCMCProposalFunctionScaleFactor[c][p] > 0)
                    throw std::runtime_error(Form("BCEngineMCMC::Metropolis: parameter %u has been unfixed without rerunning the pre-run.", p));
            }
        }
    }

    // print to screen
    BCLog::OutSummary(Form("Run Metropolis MCMC for model \"%s\" ...", GetName().data()));

    // set phase and cycle number
    fMCMCPhase = BCEngineMCMC::kMainRun;

    BCLog::OutSummary(Form(" --> Perform MCMC run with %i chains, each with %i iterations.", fMCMCNChains, fMCMCNIterationsRun));

    unsigned nwrite = UpdateFrequency(fMCMCNIterationsRun);

    // start the run
    fMCMCCurrentIteration = 0;
    while (fMCMCCurrentIteration < (int)fMCMCNIterationsRun) {

        GetNewPointMetropolis();
        EvaluateObservables();

        if (fMCMCCurrentIteration % nwrite == 0) {
            BCLog::OutDetail(Form(" --> iteration number %6i (%3.0f %%)", fMCMCCurrentIteration, (double)(fMCMCCurrentIteration) / (double)fMCMCNIterationsRun * 100.));
            if (fMCMCFlagWriteChainToFile && fMCMCTree)
                fMCMCTree->AutoSave("SaveSelf");
        }

        if (fMCMCCurrentIteration % fMCMCNLag != 0) // apply lag
            continue;

        MCMCUserIterationInterface();   // user action (overloadable)

        for (unsigned c = 0; c < fMCMCNChains; ++c)
            fMCMCStatistics[c].Update(fMCMCStates[c]);

        // fill histograms
        if (!fH1Marginalized.empty() || !fH2Marginalized.empty())
            InChainFillHistograms();

        // write chain to file
        if (fMCMCFlagWriteChainToFile)
            InChainFillTree();

    } // end run

    BCLog::OutSummary(Form(" --> Markov chains ran for %i iterations.", fMCMCNIterationsRun));

    // reset total stats
    fMCMCStatistics_AllChains.Reset();
    // add in individual chain stats
    for (unsigned c = 0; c < fMCMCStatistics.size(); ++c)
        fMCMCStatistics_AllChains += fMCMCStatistics[c];

    // print efficiencies
    if (fMCMCProposeMultivariate) {
        BCLog::OutDetail(Form(" --> Efficiencies (measured in %u iterations):", fMCMCStatistics.front().n_samples_efficiency));
        BCLog::OutDetail("       - Chain : Efficiency");
        for (unsigned c = 0; c < fMCMCNChains; ++c)
            BCLog::OutDetail(Form("           %3d :     %4.1f %%", c, 100.*fMCMCStatistics[c].efficiency[0]));
    } else {
        BCLog::OutDetail(Form(" --> Average efficiencies (measured in %u iterations):", fMCMCStatistics_AllChains.n_samples_efficiency / fMCMCNChains));
        BCLog::OutDetail(Form("       - %-*s : Efficiency", fParameters.MaxNameLength(), "Parameter"));
        for (unsigned i = 0; i < GetNParameters(); ++i) {
            if (GetParameter(i).Fixed())
                continue;
            BCLog::OutDetail(Form("         %-*s :     %4.1f %%", fParameters.MaxNameLength(), GetParameter(i).GetName().data(), 100.*fMCMCStatistics_AllChains.efficiency[i]));
        }
    }

    if (fMCMCFlagWriteChainToFile)
        UpdateParameterTree();

    BCLog::OutDetail(" --> Global mode from MCMC:");
    BCLog::OutDebug(Form(" --> Posterior value: %g", fMCMCStatistics_AllChains.probability_at_mode));
    PrintParameters(fMCMCStatistics_AllChains.modepar, BCLog::OutDetail);

    // reset counter
    fMCMCCurrentIteration = -1;

    return true;
}

// --------------------------------------------------------
void BCEngineMCMC::EvaluateObservables()
{
    if (GetNObservables() > 0)
        for (unsigned c = 0; c < fMCMCNChains; ++c)
            EvaluateObservables(c);
}

// --------------------------------------------------------
void BCEngineMCMC::EvaluateObservables(unsigned chain)
{
    if (chain > fMCMCNChains)
        return;

    UpdateChainIndex(chain);
    CalculateObservables(fMCMCStates[chain].parameters);
    for (unsigned j = 0; j < GetNObservables(); ++j)
        fMCMCStates[chain].observables[j] = GetObservable(j).Value();
}

// --------------------------------------------------------
void BCEngineMCMC::ResetResults()
{
    // reset variables
    fMCMCCurrentIteration = -1;
    fMCMCStatistics.clear();
    fMCMCStatistics_AllChains.Clear();
    fMCMCProposalFunctionScaleFactor.clear();
    fMCMCStates.clear();
    fMCMCNIterationsConvergenceGlobal = -1;
    fMCMCRValueParameters.clear();

    for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
        delete fH1Marginalized[i];

    for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
        for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
            delete fH2Marginalized[i][j];

    // clear plots
    fH1Marginalized.clear();
    fH2Marginalized.clear();

    // reset flags
    fMCMCPhase = kUnsetPhase;

    fLocalModes.clear();
}

// --------------------------------------------------------
void BCEngineMCMC::MCMCInitialize()
{
    // reset phase
    fMCMCPhase = BCEngineMCMC::kUnsetPhase;

    // reset convergence
    fMCMCNIterationsConvergenceGlobal = -1;

    // reset statistics counters
    fMCMCStatistics.assign(fMCMCNChains, BCEngineMCMC::Statistics(GetNParameters(), GetNObservables()));
    fMCMCStatistics_AllChains.Init(GetNParameters(), GetNObservables());

    // reset states holder
    fMCMCStates.assign(fMCMCNChains, ChainState(GetNObservables()));

    // rest r value holders
    fMCMCRValueParameters.assign(GetNParameters(), std::numeric_limits<double>::infinity());

    // clear info about local modes
    fLocalModes.clear();

    // clear number of multivariate proposal function covariance updates performed
    fMultivariateCovarianceUpdates = 0;

    SyncThreadStorage();

    CreateHistograms(false);

    // set scale factors
    if (fMCMCProposeMultivariate) {
        // if multivariate

        // if only one scale factor is set, use for all chains
        if (fMCMCInitialScaleFactors.size() == 1) {
            fMCMCProposalFunctionScaleFactor.assign(fMCMCNChains, fMCMCInitialScaleFactors);
        }
        // if one scale factor set per chain:
        else if (fMCMCInitialScaleFactors.size() == fMCMCNChains) {
            fMCMCInitialScaleFactors.reserve(fMCMCNChains);
            for (unsigned i = 0; i < fMCMCNChains; ++i)
                fMCMCProposalFunctionScaleFactor.push_back(std::vector<double>(1, fMCMCInitialScaleFactors[i]));
        }
        // else initialize proposal function scale factors to 2.38^2 / number of dimensions
        else {
            fMCMCProposalFunctionScaleFactor.assign(fMCMCNChains, std::vector<double>(1, 2.38 * 2.38 / GetNFreeParameters()));
        }
    }
    // else factorized proposal
    else {

        // if provided by user
        if (fMCMCInitialScaleFactors.size() == GetNParameters()) {
            fMCMCProposalFunctionScaleFactor.assign(fMCMCNChains, fMCMCInitialScaleFactors);
        }
        // else calculate from priors
        else {
            std::vector<double> temp;
            for (unsigned i = 0; i < GetNParameters(); ++i)
                if (GetParameter(i).Fixed() || GetParameter(i).GetRangeWidth() == 0)
                    temp.push_back(1);
                else {
                    double var = GetParameter(i).GetPriorVariance();
                    if (var > 0 && std::isfinite(var))
                        temp.push_back(sqrt(var) / GetParameter(i).GetRangeWidth());
                    else
                        temp.push_back(1. / sqrt(12));
                }
            fMCMCProposalFunctionScaleFactor.assign(fMCMCNChains, temp);
        }
        // set fixed parameters' scale factors to zero
        for (unsigned i = 0; i < GetNParameters(); ++i)
            if (GetParameter(i).Fixed())
                for (unsigned c = 0; c < GetNChains(); ++c)
                    fMCMCProposalFunctionScaleFactor[c][i] = 0;
    }

    // before we set the initial position and evaluate the likelihood, it's time to let the user initialize the model with #chains etc. fixed
    MCMCUserInitialize();

    /* set initial position */

    // initialize markov chain positions
    switch (fInitialPositionScheme) {

        // use range centers
        case kInitCenter : {
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                fMCMCThreadLocalStorage[c].parameters = GetParameters().GetRangeCenters();
                UpdateChainIndex(c);
                LogEval(fMCMCThreadLocalStorage[c].parameters);
                if (!std::isfinite(fMCMCThreadLocalStorage[c].log_probability))
                    throw std::runtime_error("BCEngineMCMC::MCMCInitialize : Range center as initial point yields invalid probability.");
                fMCMCStates[c] = fMCMCThreadLocalStorage[c];
            }

            break;
        }

        // uniformly distribute all coordinates in provided ranges
        case kInitRandomUniform : {
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                for (unsigned n = 0; n < fInitialPositionAttemptLimit && !std::isfinite(fMCMCThreadLocalStorage[c].log_probability); ++n) {
                    fMCMCThreadLocalStorage[c].parameters = GetParameters().GetUniformRandomValues(fMCMCThreadLocalStorage[c].rng);
                    UpdateChainIndex(c);
                    LogEval(fMCMCThreadLocalStorage[c].parameters);
                }
                if (!std::isfinite(fMCMCThreadLocalStorage[c].log_probability))
                    throw std::runtime_error(Form("BCEngineMCMC::MCMCInitialize : Could not generate uniformly distributed initial point with valid probability in %u tries.", fInitialPositionAttemptLimit));
                fMCMCStates[c] = fMCMCThreadLocalStorage[c];
            }

            break;
        }

        // use user-defined starting points
        case kInitUserDefined : {
            // check initial position vector size
            if (fMCMCInitialPosition.size() < fMCMCNChains)
                throw std::runtime_error("BCEngineMCMC::MCMCInitialize : Too few initial positions provided.");

            // copy positions and set fixed values
            // then check whether initial positions are within bounds
            // (which also checks that initial position vectors are correct size)
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                fMCMCThreadLocalStorage[c].parameters = fMCMCInitialPosition[c];
                GetParameters().ApplyFixedValues(fMCMCThreadLocalStorage[c].parameters);
                if (!GetParameters().IsWithinLimits(fMCMCThreadLocalStorage[c].parameters)) {
                    BCLog::OutDebug(Form("Initial point of chain %d", c));
                    PrintParameters(fMCMCThreadLocalStorage[c].parameters, BCLog::OutDebug);
                    throw std::runtime_error("BCEngineMCMC::MCMCInitialize : User-defined initial point is out of bounds.");
                } else {
                    UpdateChainIndex(c);
                    LogEval(fMCMCThreadLocalStorage[c].parameters);
                    if (!std::isfinite(fMCMCThreadLocalStorage[c].log_probability)) {
                        BCLog::OutDebug(Form("Initial point of chain %d", c));
                        PrintParameters(fMCMCThreadLocalStorage[c].parameters, BCLog::OutDebug);
                        throw std::runtime_error("BCEngineMCMC::MCMCInitialize : User-defined initial point yields invalid probability.");
                    }
                    fMCMCStates[c] = fMCMCThreadLocalStorage[c];
                }
            }

            break;
        }

        // randomly distribute according to factorized priors
        case kInitRandomPrior : {
            if (!GetParameters().ArePriorsSet(true))
                throw std::runtime_error("BCEngineMCMC::MCMCInitialize : Not all unfixed parameters have priors set.");

            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                for (unsigned n = 0; n < fInitialPositionAttemptLimit && !std::isfinite(fMCMCThreadLocalStorage[c].log_probability); ++n) {
                    fMCMCThreadLocalStorage[c].parameters = GetParameters().GetRandomValuesAccordingToPriors(fMCMCThreadLocalStorage[c].rng);
                    // check new point
                    if (!GetParameters().IsWithinLimits(fMCMCThreadLocalStorage[c].parameters))
                        throw std::runtime_error("BCEngineMCMC::MCMCInitialize : Could not generate random point within limits.");

                    UpdateChainIndex(c);
                    LogEval(fMCMCThreadLocalStorage[c].parameters);
                }
                if (!std::isfinite(fMCMCThreadLocalStorage[c].log_probability))
                    throw std::runtime_error(Form("BCEngineMCMC::MCMCInitialize : Could not generate initial point from prior with valid probability in %u tries.", fInitialPositionAttemptLimit));
                fMCMCStates[c] = fMCMCThreadLocalStorage[c];
            }

            break;
        }

        default:
            throw std::runtime_error("BCEngineMCMC::MCMCInitialize : No MCMC position initialization scheme specified.");
    } // (switch)

    if (fMCMCStates.empty())
        throw std::runtime_error("BCEngineMCMC::MCMCInitialize failed.");

    for (unsigned c = 0; c < fMCMCStates.size(); ++c)
        if (fMCMCStates[c].parameters.empty())
            throw std::runtime_error("BCEngineMCMC::MCMCInitialize failed.");
}

// ------------------------------------------------------------
void BCEngineMCMC::SetFillHistogram(int x, int y, bool flag)
{
    // check indices
    if (x >= (int)fParameters.Size() || - x > (int)fObservables.Size())
        return;
    if (y >= (int)fParameters.Size() || - y > (int)fObservables.Size())
        return;

    if (flag) {                 // adding
        // check for combination already in list
        for (unsigned i = 0; i < fRequestedH2.size(); ++i)
            if (fRequestedH2[i].first == x && fRequestedH2[i].second == y)
                return;
        fRequestedH2.push_back(std::make_pair(x, y));
    } else {                    // removing
        for (int i = fRequestedH2.size() - 1; i >= 0; --i)
            if (fRequestedH2[i].first == x && fRequestedH2[i].second == y)
                fRequestedH2.erase(fRequestedH2.begin() + i);
    }
}

// ------------------------------------------------------------
void BCEngineMCMC::CreateHistograms(bool rescale_ranges)
{
    // clear existing histograms
    for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
        delete fH1Marginalized[i];
    fH1Marginalized.assign(GetNVariables(), NULL);

    for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
        for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
            delete fH2Marginalized[i][j];
    fH2Marginalized.assign(GetNVariables(), std::vector<TH2*>(GetNVariables(), NULL));

    // store old bounds, rescale if rescaling:
    std::vector<std::pair<double, double> > original_bounds;
    original_bounds.reserve(GetNVariables());
    for (unsigned i = 0; i < GetNVariables(); ++i) {
        original_bounds.push_back(std::make_pair(GetVariable(i).GetLowerLimit(), GetVariable(i).GetUpperLimit()));
        if (i < GetNParameters() && GetParameter(i).Fixed())
            continue;
        if (rescale_ranges) {
            if (i < fMCMCStatistics_AllChains.minimum.size() && std::isfinite(fMCMCStatistics_AllChains.minimum[i]))
                GetVariable(i).SetLowerLimit(fMCMCStatistics_AllChains.minimum[i]);
            if (i < fMCMCStatistics_AllChains.maximum.size() && std::isfinite(fMCMCStatistics_AllChains.maximum[i]))
                GetVariable(i).SetUpperLimit(fMCMCStatistics_AllChains.maximum[i]);
            if (fHistogramRescalePadding > 0) {
                // calculate enlargement factor
                double range_width_rescale = GetVariable(i).GetRangeWidth() * fHistogramRescalePadding;
                // push bounds out, but not beyond original parameter bounds
                GetVariable(i).SetLowerLimit(std::max<double>(GetVariable(i).GetLowerLimit() - range_width_rescale, original_bounds.back().first));
                GetVariable(i).SetUpperLimit(std::min<double>(GetVariable(i).GetUpperLimit() + range_width_rescale, original_bounds.back().second));
            }
        }
    }

    // define 1-dimensional histograms for marginalization
    int filling = 0;

    for (unsigned i = 0; i < GetNVariables(); ++i)
        if (GetVariable(i).FillH1()) {
            if (i < GetNParameters()) { // parameter
                if (!GetParameter(i).Fixed()) {
                    fH1Marginalized[i] = GetVariable(i).CreateH1(Form("h1_%s_parameter_%s", GetSafeName().data() , GetParameter(i).GetSafeName().data()));
                    ++filling;
                }
            } else {                  // user-defined observable
                fH1Marginalized[i] = GetVariable(i).CreateH1(Form("h1_%s_observable_%s", GetSafeName().data() , GetVariable(i).GetSafeName().data()));
                ++filling;
            }
        }

    if (filling == 0) // if filling no 1D histograms, clear vector
        fH1Marginalized.clear();

    // define 2D histograms for marginalization
    filling = 0;

    // requested 2D histograms in order requested:
    for (std::vector<std::pair<int, int> >::const_iterator h = fRequestedH2.begin(); h != fRequestedH2.end(); ++h) {
        if (h->first >= 0 && GetParameter(h->first).Fixed())
            continue;
        if (h->second >= 0 && GetParameter(h->second).Fixed())
            continue;

        unsigned i = (h->first >= 0)  ? h->first  : -(h->first + 1);
        unsigned j = (h->second >= 0) ? h->second : -(h->second + 1);

        if (h->first >= 0) {
            if (i >= GetNParameters())
                continue;

            if (h->second >= 0) { // par vs par
                if (j >= GetNParameters())
                    continue;
                fH2Marginalized[i][j] = GetParameter(i).CreateH2(Form("h2_%s_par_%s_vs_par_%s", GetSafeName().data(), GetParameter(j).GetSafeName().data(), GetParameter(i).GetSafeName().data()), GetParameter(j));
            } else {              // obs vs par
                if (j >= GetNObservables())
                    continue;
                fH2Marginalized[i][j + GetNParameters()] = GetParameter(i).CreateH2(Form("h2_%s_obs_%s_vs_par_%s", GetSafeName().data(), GetObservable(j).GetSafeName().data(), GetParameter(i).GetSafeName().data()), GetObservable(j));
            }

        } else {
            if (i >= GetNObservables())
                continue;

            if (h->second >= 0) { // par vs obs
                if (j >= GetNParameters())
                    continue;
                fH2Marginalized[i + GetNParameters()][j] = GetObservable(i).CreateH2(Form("h2_%s_par_%s_vs_obs_%s", GetSafeName().data(), GetParameter(j).GetSafeName().data(), GetObservable(i).GetSafeName().data()), GetParameter(j));
            } else { // obs vs obs
                if (j >= GetNObservables())
                    continue;
                fH2Marginalized[i + GetNParameters()][j + GetNParameters()] = GetObservable(i).CreateH2(Form("h2_%s_obs_%s_vs_obs_%s", GetSafeName().data(), GetObservable(j).GetSafeName().data(), GetObservable(i).GetSafeName().data()), GetObservable(j));
            }
        }
        ++filling;
    }

    // automatically produced combinations, using pars' and obvs' FillH2():

    // parameter i as abscissa:
    for (unsigned i = 0; i < GetNParameters(); ++i) {
        if (GetParameter(i).Fixed() || !GetParameter(i).FillH2())
            continue;

        // parameter j as ordinate
        for (unsigned j = i + 1; j < GetNParameters(); ++j)
            if (!GetParameter(j).Fixed() && GetParameter(j).FillH2() && !fH2Marginalized[i][j]) {
                fH2Marginalized[i][j] = GetParameter(i).CreateH2(Form("h2_%s_par_%s_vs_par_%s", GetSafeName().data(), GetParameter(j).GetSafeName().data(), GetParameter(i).GetSafeName().data()), GetParameter(j));
                ++filling;
            }
        // user-defined observable j as ordinate
        for (unsigned j = 0; j < GetNObservables(); ++j)
            if (GetObservable(j).FillH2() && !fH2Marginalized[i][j + GetNParameters()]) {
                fH2Marginalized[i][j + GetNParameters()] = GetParameter(i).CreateH2(Form("h2_%s_obs_%s_vs_par_%s", GetSafeName().data(), GetObservable(j).GetSafeName().data(), GetParameter(i).GetSafeName().data()), GetObservable(j));
                ++filling;
            }
    }

    // user-defined observable i as abscissa
    for (unsigned i = 0; i < GetNObservables(); ++i) {
        if (!GetObservable(i).FillH2())
            continue;

        // user-defined observable j as ordinate
        for (unsigned j = i + 1; j < GetNObservables(); ++j)
            if (GetObservable(j).FillH2() && !fH2Marginalized[i + GetNParameters()][j + GetNParameters()]) {
                fH2Marginalized[i + GetNParameters()][j + GetNParameters()] = GetObservable(i).CreateH2(Form("h2_%s_obs_%s_vs_obs_%s", GetSafeName().data(), GetObservable(j).GetSafeName().data(), GetObservable(i).GetSafeName().data()), GetObservable(j));
                ++filling;
            }
    }

    if (filling == 0) // if filling no 2D histograms, clear vector
        fH2Marginalized.clear();

    // restore bounds
    for (unsigned i = 0; i < GetNVariables(); ++i)
        GetVariable(i).SetLimits(original_bounds[i].first, original_bounds[i].second);
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintSummary() const
{
    PrintModelSummary();
    BCLog::OutSummary("");
    PrintBestFitSummary();
    BCLog::OutSummary("");
    PrintMarginalizationSummary();
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintModelSummary() const
{
    BCLog::OutSummary("");
    BCLog::OutSummary(" -----------------------------------------------------");
    BCLog::OutSummary(" Summary");
    BCLog::OutSummary(" -----------------------------------------------------");
    BCLog::OutSummary("");

    BCLog::OutSummary(" Model summary");
    BCLog::OutSummary(" =============");
    BCLog::OutSummary(" Model: " + GetName());
    BCLog::OutSummary(Form(" Number of parameters: %u", GetNParameters()));

    if (!fParameters.Empty()) {
        BCLog::OutSummary("");
        BCLog::OutSummary(" List of parameters and ranges:");
        fParameters.PrintSummary();
    }

    if (!fObservables.Empty()) {
        BCLog::OutSummary("");
        BCLog::OutSummary(" List of observables and ranges:");
        fObservables.PrintSummary();
    }
    BCLog::OutSummary("");
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintBestFitSummary() const
{
    if (GetBestFitParameters().size() != GetNParameters() && GetBestFitParameters().size() != GetNVariables()) {
        BCLog::OutSummary("No best fit information available.");
        return;
    }

    BCLog::OutSummary(" Best Fit Results");
    BCLog::OutSummary(" ===========================");
    BCLog::OutSummary(Form(" Log of the maximum posterior: %f", GetLogMaximum()));
    BCLog::OutSummary("");
    BCLog::OutSummary(" Global mode:");

    for (unsigned i = 0; i < GetBestFitParameters().size(); ++i)
        BCLog::OutSummary(GetBestFitSummary(i));
}

// ---------------------------------------------------------
std::string BCEngineMCMC::GetBestFitSummary(unsigned i) const
{
    if (i >= GetNVariables())
        return std::string("");

    unsigned n = (int)log10(GetNVariables()) + 1;
    std::string par_summary = Form(" %*d) %10s \"%s\"%*s : %.*g", n, i, GetVariable(i).GetPrefix().data(),
                                   GetVariable(i).GetName().data(), (int)(GetMaximumParameterNameLength() - GetVariable(i).GetName().length()), "",
                                   GetVariable(i).GetPrecision(), GetBestFitParameters()[i]);

    if (i < GetNParameters() && GetParameter(i).Fixed())
        par_summary += " (fixed)";

    return par_summary;
}

// ---------------------------------------------------------
void BCEngineMCMC::PrintMarginalizationSummary() const
{
    BCLog::OutSummary(" Results of the marginalization");
    BCLog::OutSummary(" ==============================");

    if (fMCMCNIterationsConvergenceGlobal < 0) {
        // give warning if MCMC did not converge
        BCLog::OutSummary(" WARNING: the Markov Chain did not converge!");
        BCLog::OutSummary(" Be cautious using the following results!");
        BCLog::OutSummary("");
    }

    BCLog::OutSummary(" List of variables and properties of the marginalized distributions:");
    BCLog::OutSummary("");

    for (unsigned i = 0; i < GetNVariables(); ++i) {
        std::string par_summary = Form("  (%u) ", i) + GetVariable(i).GetPrefix() + " \"" + GetVariable(i).GetName() + "\" :";

        if (i < GetNParameters() && GetParameter(i).Fixed()) {
            par_summary += Form(" fixed at %.*g", GetVariable(i).GetPrecision(), GetParameter(i).GetFixedValue());
            BCLog::OutSummary(par_summary);

        } else if (!MarginalizedHistogramExists(i)) {
            par_summary += " histogram does not exist.";
            BCLog::OutSummary(par_summary);

        } else {
            BCLog::OutSummary(par_summary);
            GetMarginalized(i).PrintSummary("      ", GetVariable(i).GetPrecision(), std::vector<double>(1, 0.68));
        }

        BCLog::OutSummary("");
    }

    if (fMCMCPhase == kMainRun) {

        BCLog::OutSummary(" Status of the MCMC");
        BCLog::OutSummary(" ==================");

        if (GetNIterationsConvergenceGlobal() > 0) {
            BCLog::OutSummary(" Convergence reached:                    yes");
            BCLog::OutSummary(Form(" Number of iterations until convergence: %d", fMCMCNIterationsConvergenceGlobal));
        } else
            BCLog::OutSummary(" Convergence reached:                    no");

        BCLog::OutSummary(Form(" Number of chains:                       %u", fMCMCNChains));
        BCLog::OutSummary(Form(" Number of iterations per chain:         %u", fMCMCNIterationsRun));

        if (fMCMCProposeMultivariate) {
            BCLog::OutDetail(Form(" Scale factors and efficiencies (measured in last %u iterations):", fMCMCStatistics.front().n_samples_efficiency));
            BCLog::OutDetail(" Chain : Scale factor    Efficiency");
            for (unsigned c = 0; c < fMCMCNChains; ++c) {
                BCLog::OutDetail(Form("   %3d :       %-6.4g        %4.1f %%", c, fMCMCProposalFunctionScaleFactor[c][0], 100.*fMCMCStatistics[c].efficiency[0]));
            }

        } else {
            BCLog::OutDetail(" Average scale factors and efficiencies:");
            BCLog::OutDetail(Form(" %-*s : Scale factor    Efficiency", fParameters.MaxNameLength(), "Parameter"));
            for (unsigned i = 0; i < GetNParameters(); ++i) {
                if (GetParameter(i).Fixed())
                    continue;
                double scalefactor = 0;
                for (unsigned j = 0; j < fMCMCNChains; ++j)
                    scalefactor += fMCMCProposalFunctionScaleFactor[j][i] / fMCMCNChains;
                BCLog::OutDetail(Form(" %-*s :          % 6.4g %%        %4.1f %%", fParameters.MaxNameLength(), GetParameter(i).GetName().data(), 100.*scalefactor, 100.*fMCMCStatistics_AllChains.efficiency[i]));
            }
        }
    }
}


// ---------------------------------------------------------
void BCEngineMCMC::PrintParameters(const std::vector<double>& P, void (*output)(const std::string&)) const
{
    if (P.size() != GetNParameters() && P.size() != GetNVariables())
        return;

    for (unsigned i = 0; i < P.size(); ++i)
        output(Form("          %-*s :   % .*g", GetMaximumParameterNameLength(P.size() > GetNParameters()), GetVariable(i).GetName().data(), GetVariable(i).GetPrecision(), P[i]));
}

// ---------------------------------------------------------
std::vector<unsigned> BCEngineMCMC::GetH1DPrintOrder() const
{
    std::vector<unsigned> H1Indices;
    for (unsigned i = 0; i < GetNVariables(); ++i)
        if (MarginalizedHistogramExists(i))
            H1Indices.push_back(i);

    return H1Indices;
}

// ---------------------------------------------------------
std::vector<std::pair<unsigned, unsigned> > BCEngineMCMC::GetH2DPrintOrder() const
{
    // create vector for ordering 2D hists properly
    std::vector<std::pair<unsigned, unsigned> > H2Coords;
    H2Coords.reserve(GetNVariables()*GetNVariables() - 1);
    // par vs par
    for (unsigned i = 0; i < GetNParameters(); ++i)
        for (unsigned j = 0; j < GetNParameters(); ++j)
            if (MarginalizedHistogramExists(i, j))
                H2Coords.push_back(std::make_pair(i, j));
    // obs vs par
    for (unsigned i = 0; i < GetNParameters(); ++i)
        for (unsigned j = GetNParameters(); j < GetNVariables(); ++j)
            if (MarginalizedHistogramExists(i, j))
                H2Coords.push_back(std::make_pair(i, j));
    // par vs obs
    for (unsigned i = GetNParameters(); i < GetNVariables(); ++i)
        for (unsigned j = 0; j < GetNParameters(); ++j)
            if (MarginalizedHistogramExists(i, j))
                H2Coords.push_back(std::make_pair(i, j));
    // obs vs obs
    for (unsigned i = GetNParameters(); i < GetNVariables(); ++i)
        for (unsigned j = GetNParameters(); j < GetNVariables(); ++j)
            if (MarginalizedHistogramExists(i, j))
                H2Coords.push_back(std::make_pair(i, j));

    return H2Coords;
}

// ---------------------------------------------------------
unsigned BCEngineMCMC::PrintAllMarginalized(const std::string& filename, unsigned hdiv, unsigned vdiv) const
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : No variables defined!");
        return 0;
    }

    if (fH1Marginalized.empty() && fH2Marginalized.empty()) {
        BCLog::OutError("BCEngineMCMC::PrintAllMarginalized : Marginalized distributions not stored.");
        return 0;
    }
    std::string newFilename(filename);
    BCAux::DefaultToPDF(newFilename);
    if (newFilename.empty())
        return 0;

    // Find nonempty H1's
    std::vector<unsigned> H1Indices = GetH1DPrintOrder();
    std::vector<BCH1D> h1;
    h1.reserve(H1Indices.size());
    for (unsigned i = 0; i < H1Indices.size(); ++i) {
        if (GetMarginalizedHistogram(H1Indices[i])->Integral() == 0) { // histogram was never filled in range
            BCLog::OutWarning(Form("BCEngineMCMC::PrintAllMarginalized : 1D Marginalized histogram for \"%s\" is empty; printing is skipped.", GetVariable(H1Indices[i]).GetName().data()));
            continue;
        }
        h1.push_back(GetMarginalized(H1Indices[i]));
        if (h1.back().Valid())
            h1.back().CopyOptions(fBCH1DdrawingOptions);
        else
            h1.pop_back();
    }

    // Find nonempty H2's
    std::vector<std::pair<unsigned, unsigned> > H2Coords = GetH2DPrintOrder();
    std::vector<BCH2D> h2;
    h2.reserve(H2Coords.size());
    for (unsigned k = 0; k < H2Coords.size(); ++k) {
        unsigned i = H2Coords[k].first;
        unsigned j = H2Coords[k].second;
        if (fH2Marginalized[i][j]->Integral() == 0) { // histogram was never filled in range
            BCLog::OutWarning(Form("BCEngineMCMC::PrintAllMarginalized : 2D Marginalized histogram for \"%s\":\"%s\" is empty; printing is skipped.", GetVariable(i).GetName().data(), GetVariable(i).GetName().data()));
            continue;
        }
        h2.push_back(GetMarginalized(i, j));
        if (h2.back().Valid())
            h2.back().CopyOptions(fBCH2DdrawingOptions);
        else
            h2.pop_back();
    }

    return BCAux::PrintPlots(h1, h2, newFilename, hdiv, vdiv);
}

// ---------------------------------------------------------
unsigned BCEngineMCMC::PrintParameterPlot(const std::string& filename, int npar, double interval_content, std::vector<double> quantiles, bool rescale_ranges) const
{
    std::string newFilename(filename);
    BCAux::DefaultToPDF(newFilename);
    if (newFilename.empty())
        return 0;

    TCanvas c_par("c_parplot_init");
    c_par.Print(Form("%s[", newFilename.data()));
    c_par.cd();
    c_par.SetTicky(1);
    c_par.SetFrameLineWidth(0);
    c_par.SetFrameLineColor(0);

    if (npar <= 0) // all parameters on one page, all user-defined observables on the next
        npar = std::max<int> (GetNParameters(), GetNObservables());

    unsigned pages_printed = 0;

    // parameters first
    for (unsigned i = 0; i < GetNParameters(); i += npar)
        if (DrawParameterPlot(i, std::min<int>(npar, GetNParameters() - i), interval_content, quantiles, rescale_ranges)) {
            c_par.Print(newFilename.data());
            c_par.Clear();
            ++pages_printed;
        }

    // then user-defined observables
    for (unsigned i = GetNParameters(); i < GetNVariables(); i += npar)
        if (DrawParameterPlot(i, std::min<int>(npar, GetNVariables() - i), interval_content, quantiles, rescale_ranges)) {
            c_par.Print(newFilename.data());
            c_par.Clear();
            ++pages_printed;
        }

    c_par.Print(Form("%s]", newFilename.data()));
    return pages_printed > 0;
}

// ---------------------------------------------------------
bool BCEngineMCMC::DrawParameterPlot(unsigned i0, unsigned npar, double interval_content, std::vector<double> quantiles, bool rescale_ranges) const
{

    // if npar==0, print all remaining observables
    unsigned i1 = (npar > 0 && i0 + npar <= GetNVariables()) ? i0 + npar : GetNVariables();

    if (i1 <= i0) {
        BCLog::OutError(Form("BCSummaryTool::PrintParameterPlot : invalid parameter range [%d, %d)", i0, i1));
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
    original_min.reserve(i1 - i0);
    original_max.reserve(i1 - i0);
    if (rescale_ranges) {
        for (unsigned i = i0; i < i1; ++i) {
            BCVariable& var = const_cast<BCVariable&>(GetVariable(i));
            original_min.push_back(GetVariable(i).GetLowerLimit());
            original_max.push_back(GetVariable(i).GetUpperLimit());
            if (i < fMCMCStatistics_AllChains.minimum.size() && std::isfinite(fMCMCStatistics_AllChains.minimum[i]))
                var.SetLowerLimit(fMCMCStatistics_AllChains.minimum[i]);
            if (i < fMCMCStatistics_AllChains.maximum.size() && std::isfinite(fMCMCStatistics_AllChains.maximum[i]))
                var.SetUpperLimit(fMCMCStatistics_AllChains.maximum[i]);
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

        if (i < GetNParameters() && GetParameter(i).Fixed())
            continue;

        // Global Mode
        x_i_bf.push_back(i);
        global_mode.push_back(GetVariable(i).PositionInRange(GetBestFitParameters()[i]));
        mean.push_back(std::numeric_limits<double>::infinity());
        rms.push_back(0);

        if (i < fH1Marginalized.size() && fH1Marginalized[i]) {
            BCH1D bch1d_temp = GetMarginalized(i);
            if (bch1d_temp.Valid()) {

                x_i.push_back(i);

                // quantiles
                x_quantiles.insert(x_quantiles.end(), quantiles.size(), i);
                std::vector<double> q(quantiles.size(), 0);
                bch1d_temp.GetHistogram()->GetQuantiles(quantiles.size(), &q[0], &quantiles[0]);
                for (unsigned j = 0; j < quantiles.size(); ++j)
                    quantile_vals.push_back(GetVariable(i).PositionInRange(q[j]));

                local_mode.push_back(GetVariable(i).PositionInRange(bch1d_temp.GetLocalMode(0)));
                mean.back() = GetVariable(i).PositionInRange(bch1d_temp.GetHistogram()->GetMean());
                rms.back() = bch1d_temp.GetHistogram()->GetRMS() / GetVariable(i).GetRangeWidth();

                // smallest interval
                BCH1D::BCH1DSmallestInterval SI = bch1d_temp.GetSmallestIntervals(interval_content);
                if (SI.intervals.empty()) {
                    interval_lo.push_back(0);
                    interval_hi.push_back(0);
                } else {
                    interval_lo.push_back(fabs(SI.intervals.front().mode - SI.intervals.front().xmin) / GetVariable(i).GetRangeWidth());
                    interval_hi.push_back(fabs(SI.intervals.front().mode - SI.intervals.front().xmax) / GetVariable(i).GetRangeWidth());
                }
            }
        }

        // use chain statistics if they exist:
        if (i < fMCMCStatistics_AllChains.mean.size() && std::isfinite(fMCMCStatistics_AllChains.mean[i]))
            mean.back() = GetVariable(i).PositionInRange(fMCMCStatistics_AllChains.mean[i]);
        if (i < fMCMCStatistics_AllChains.variance.size() && std::isfinite(fMCMCStatistics_AllChains.variance[i]))
            rms.back() = sqrt(fMCMCStatistics_AllChains.variance[i]) / GetVariable(i).GetRangeWidth();
    }

    if (x_i.empty() && x_i_bf.empty())
        return false;

    /////////////////////////
    // Draw it all

    // axes have a name, the other objects don't, so they are not registered
    // with ROOT and we don't need a guard

    // Create, label, and draw axes
    TH2D* hist_axes;
    {
        BCAux::RootSideEffectGuard g;
        hist_axes = new TH2D(Form("h2_axes_parplot_%s_%d_%d", GetSafeName().data(), i0, i1), "",  //";;Scaled range [a.u.]",
                             i1 - i0, i0 - 0.5, i1 - 0.5, 10, -0.05 + 1e-3, 1.05 - 1e-3);
    }
    fObjectTrash.Put(hist_axes);
    hist_axes->SetStats(kFALSE);
    hist_axes->GetXaxis()->SetAxisColor(0);
    hist_axes->GetXaxis()->SetLabelOffset(0.015);
    hist_axes->GetXaxis()->SetLabelSize(std::max<double>(0.01, 0.05 * std::min<double>(1, 4. / hist_axes->GetNbinsX())));
    hist_axes->GetXaxis()->SetTickLength(0.0);
    hist_axes->GetYaxis()->SetLabelSize(0);

    // set bin labels
    for (int i = 0; i < hist_axes->GetNbinsX(); ++i)
        hist_axes->GetXaxis()->SetBinLabel(i + 1, GetVariable(i0 + i).GetLatexNameWithUnits().data());
    hist_axes->Draw("");

    // Draw lines
    TLine* line = new TLine();
    fObjectTrash.Put(line);
    line->SetLineColor(kBlack);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->DrawLine(hist_axes->GetXaxis()->GetXmin(), 0.0, hist_axes->GetXaxis()->GetXmax(), 0.0);
    line->DrawLine(hist_axes->GetXaxis()->GetXmin(), 1.0, hist_axes->GetXaxis()->GetXmax(), 1.0);

    // Mark parameter ranges
    TLatex* latex = new TLatex();
    fObjectTrash.Put(latex);
    latex->SetTextSize(0.02);
    for (unsigned i = i0; i < i1; ++i)
        if (i < GetNParameters() && GetParameter(i).Fixed()) {
            latex->SetTextAlign(22);
            latex->DrawLatex((double)i, 0.52, "Fixed at");
            latex->DrawLatex((double)i, 0.47, Form("%.*g", GetVariable(i).GetPrecision(), GetParameter(i).GetFixedValue()));
        } else {
            latex->SetTextAlign(21);
            latex->DrawLatex((double)i,  1.015, Form("%+.*g", GetVariable(i).GetPrecision(), GetVariable(i).GetUpperLimit()));
            latex->SetTextAlign(23);
            latex->DrawLatex((double)i, -0.015, Form("%+.*g", GetVariable(i).GetPrecision(), GetVariable(i).GetLowerLimit()));
        }

    // create legend
    TLegend* legend = new TLegend();
    fObjectTrash.Put(legend);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetNColumns(2);
    legend->SetTextAlign(12);
    legend->SetTextSize(0.02 * gPad->GetWNDC());

    if (!x_i.empty()) {

        // Smallest Interval
        std::vector<double> x_i_err(x_i.size(), 0.5);
        TGraphAsymmErrors* graph_intervals = new TGraphAsymmErrors(x_i.size(), x_i.data(), local_mode.data(), x_i_err.data(), x_i_err.data(), interval_lo.data(), interval_hi.data());
        fObjectTrash.Put(graph_intervals);
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
            std::vector<double> quantiles_err(x_quantiles.size(), 0.5);
            TGraphErrors* graph_quantiles = new TGraphErrors(x_quantiles.size(), x_quantiles.data(), quantile_vals.data(), quantiles_err.data(), 0);
            fObjectTrash.Put(graph_quantiles);
            graph_quantiles->SetMarkerSize(0);
            graph_quantiles->SetLineColor(38);
            graph_quantiles->SetLineStyle(2);
            graph_quantiles->Draw("SAMEZ");
            std::string quantiles_text = "Quantiles (";
            for (unsigned i = 0; i < quantiles.size() - 1; ++i)
                quantiles_text += Form("%.0f%%, ", quantiles[i] * 100);
            quantiles_text += (quantiles.size() > 0) ? Form("%.0f%%)", quantiles.back() * 100) : "none)";
            legend->AddEntry(graph_quantiles, quantiles_text.c_str(), "L");
        }

        // Means & RMSs
        TGraphErrors* graph_mean = new TGraphErrors(x_i.size(), x_i.data(), mean.data(), 0, rms.data());
        fObjectTrash.Put(graph_mean);
        graph_mean->SetMarkerColor(kBlack);
        graph_mean->SetMarkerStyle(21);
        graph_mean->SetMarkerSize(1 * gPad->GetWNDC());
        graph_mean->Draw("SAMEP");

        legend->AddEntry(graph_mean, "Mean and RMS", "LEP");
        legend->AddEntry(graph_intervals, Form("Smallest %.0f%% interval and local mode", 100.*interval_content), "FL");
    }

    // Global Modes
    if (!x_i_bf.empty()) {
        TGraph* graph_mode = new TGraph(x_i_bf.size(), x_i_bf.data(), global_mode.data());
        fObjectTrash.Put(graph_mode);
        graph_mode->SetMarkerColor(kRed);
        graph_mode->SetMarkerStyle(20);
        graph_mode->SetMarkerSize(1 * gPad->GetWNDC());
        graph_mode->Draw("SAMEP");
        legend->AddEntry(graph_mode, "Global mode", "P");
    }

    gPad->SetTopMargin(0.02);

    // place legend on top of histogram
#if ROOTVERSION >= 6000000
    legend->SetX1(gPad->GetLeftMargin());
    legend->SetX2(1. - gPad->GetRightMargin());
    double y1 = gPad->GetTopMargin() + legend->GetTextSize() * legend->GetNRows();
    legend->SetY1(1. - y1);
    legend->SetY2(1. - gPad->GetTopMargin());
#else
    legend->SetX1NDC(gPad->GetLeftMargin());
    legend->SetX2NDC(1. - gPad->GetRightMargin());
    double y1 = gPad->GetTopMargin() + legend->GetTextSize() * legend->GetNRows();
    legend->SetY1NDC(1. - y1);
    legend->SetY2NDC(1. - gPad->GetTopMargin());
#endif
    legend->Draw("SAME");

    gPad->SetTopMargin(y1 + 0.01);

    gPad->RedrawAxis();
    gPad->Update();

    // restore ranges
    for (unsigned i = 0; i < original_min.size(); ++i)
        const_cast<BCVariable&>(GetVariable(i0 + i)).SetLimits(original_min[i], original_max[i]);

    // no error
    return true;
}

// ---------------------------------------------------------
bool BCEngineMCMC::PrintCorrelationMatrix(const std::string& filename) const
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCEngineMCMC::PrintCorrelationMatrix : No variables defined!");
        return 0;
    }

    // create histogram
    TH2D hist_corr(Form("hist_correlation_matrix_%s", GetSafeName().data()), ";;", GetNVariables(), -0.5, GetNVariables() - 0.5, GetNVariables(), -0.5, GetNVariables() - 0.5);
    hist_corr.SetStats(false);
    hist_corr.GetXaxis()->SetTickLength(0.0);
    hist_corr.GetYaxis()->SetTickLength(0.0);
    hist_corr.GetXaxis()->SetLabelSize(0);
    hist_corr.GetYaxis()->SetLabelSize(0);

    // vector of unfilled values:
    std::vector<std::pair<unsigned, unsigned> > unfilled;

    // fill histogram
    for (unsigned i = 0; i < GetNVariables(); ++i) {
        hist_corr.SetBinContent(i + 1, GetNVariables() - i, 1);

        double var_i = (i < fMCMCStatistics_AllChains.variance.size()) ? fMCMCStatistics_AllChains.variance[i] : std::numeric_limits<double>::infinity();
        for (unsigned j = i + 1; j < GetNVariables(); ++j) {
            double var_j = (j < fMCMCStatistics_AllChains.variance.size()) ? fMCMCStatistics_AllChains.variance[j] : std::numeric_limits<double>::infinity();
            double covar_ij = (i < fMCMCStatistics_AllChains.covariance.size() && j < fMCMCStatistics_AllChains.covariance[i].size()) ? fMCMCStatistics_AllChains.covariance[i][j] : std::numeric_limits<double>::infinity();
            double corr_ij = std::numeric_limits<double>::infinity();
            if (std::isfinite(covar_ij) && std::isfinite(var_i) && std::isfinite(var_j))
                corr_ij = covar_ij / sqrt(var_i * var_j);
            else if (i < fH2Marginalized.size() && j < fH2Marginalized[i].size() && fH2Marginalized[i][j])
                corr_ij = fH2Marginalized[i][j]->GetCorrelationFactor();
            else if (j < fH2Marginalized.size() && i < fH2Marginalized[j].size() && fH2Marginalized[j][i])
                corr_ij = fH2Marginalized[j][i]->GetCorrelationFactor();
            if (std::isfinite(corr_ij)) {
                hist_corr.SetBinContent(i + 1, GetNVariables() - j, corr_ij);
                hist_corr.SetBinContent(j + 1, GetNVariables() - i, corr_ij);
            } else {
                unfilled.push_back(std::make_pair(i, j));
            }
        }
    }

    // print to file
    TCanvas c_corr("c_corr_matrix");
    c_corr.cd();

    double text_size = std::max<double>(0.005, 0.02 * std::min<double>(1., 5. / GetNVariables()));

    TLatex xlabel;
    xlabel.SetTextSize(text_size);
    xlabel.SetTextAlign(22);

    TLatex ylabel;
    ylabel.SetTextSize(text_size);
    ylabel.SetTextAlign(22);
    ylabel.SetTextAngle(90);

    TLatex corr_number;
    corr_number.SetTextSize(text_size);
    corr_number.SetTextAlign(22);

    gStyle->SetPalette(54);
    gStyle->SetPaintTextFormat("+.2g");
    hist_corr.GetZaxis()->SetRangeUser(-1, 1);
    hist_corr.GetZaxis()->SetDecimals(true);
    hist_corr.GetZaxis()->SetLabelSize(text_size);
    hist_corr.Draw("colz");

    // Draw labels and correlations
    for (int i = 1; i <= hist_corr.GetNbinsX(); ++i) {
        // labels
        xlabel.DrawLatex(hist_corr.GetXaxis()->GetBinCenter(i),
                         hist_corr.GetYaxis()->GetXmax() + 12e-2,
                         GetVariable(i - 1).GetLatexNameWithUnits().data());

        ylabel.DrawLatex(hist_corr.GetXaxis()->GetXmin() - 12e-2,
                         hist_corr.GetYaxis()->GetBinCenter(GetNVariables() - i + 1),
                         GetVariable(i - 1).GetLatexNameWithUnits().data());
        for (int j = 1; j <= hist_corr.GetNbinsY(); ++j) {
            if (hist_corr.GetBinContent(i, j) >= 0)
                corr_number.SetTextColor(kBlack);
            else
                corr_number.SetTextColor(kWhite);
            corr_number.DrawLatex(hist_corr.GetXaxis()->GetBinCenter(i),
                                  hist_corr.GetYaxis()->GetBinCenter(j),
                                  Form("%+.2g", hist_corr.GetBinContent(i, j)));
        }
    }

    // Blank out empty squares
    TBox bcorr;
    bcorr.SetFillColor(kWhite);
    for (unsigned i = 0; i < unfilled.size(); ++i) {
        bcorr.DrawBox(hist_corr.GetXaxis()->GetBinLowEdge(unfilled[i].first + 1), hist_corr.GetYaxis()->GetBinLowEdge(GetNVariables() - unfilled[i].second),
                      hist_corr.GetXaxis()->GetBinUpEdge (unfilled[i].first + 1), hist_corr.GetYaxis()->GetBinUpEdge (GetNVariables() - unfilled[i].second));
        bcorr.DrawBox(hist_corr.GetXaxis()->GetBinLowEdge(unfilled[i].second + 1), hist_corr.GetYaxis()->GetBinLowEdge(GetNVariables() - unfilled[i].first),
                      hist_corr.GetXaxis()->GetBinUpEdge (unfilled[i].second + 1), hist_corr.GetYaxis()->GetBinUpEdge (GetNVariables() - unfilled[i].first));
    }

    // redraw top and right lines
    TLine lA;
    lA.DrawLine(hist_corr.GetXaxis()->GetXmin(), hist_corr.GetYaxis()->GetXmax(), hist_corr.GetXaxis()->GetXmax(), hist_corr.GetYaxis()->GetXmax());
    lA.DrawLine(hist_corr.GetXaxis()->GetXmax(), hist_corr.GetYaxis()->GetXmin(), hist_corr.GetXaxis()->GetXmax(), hist_corr.GetYaxis()->GetXmax());
    // draw line between parameters and user-defined observables
    if (GetNObservables() > 0) {
        lA.DrawLine(hist_corr.GetXaxis()->GetXmin() - 0.40, hist_corr.GetYaxis()->GetBinLowEdge(hist_corr.GetNbinsY() - GetNParameters() + 1),
                    hist_corr.GetXaxis()->GetBinUpEdge(GetNParameters()), hist_corr.GetYaxis()->GetBinLowEdge(hist_corr.GetNbinsY() - GetNParameters() + 1));
        lA.DrawLine(hist_corr.GetXaxis()->GetBinUpEdge(GetNParameters()), hist_corr.GetYaxis()->GetBinLowEdge(hist_corr.GetNbinsY() - GetNParameters() + 1),
                    hist_corr.GetXaxis()->GetBinUpEdge(GetNParameters()), hist_corr.GetYaxis()->GetXmax() + 0.45);
    }

    gPad->RedrawAxis();
    c_corr.Print(filename.data());

    // no error
    return true;
}

// ---------------------------------------------------------
bool BCEngineMCMC::PrintCorrelationPlot(const std::string& filename, bool include_observables) const
{
    if (GetNVariables() == 0) {
        BCLog::OutError("BCEngineMCMC::PrintCorrelationPlot : No variables defined!");
        return 0;
    }

    // Array of indices for which any marginalizations were stored
    std::vector<unsigned> I;
    unsigned n = (include_observables) ? GetNVariables() : GetNParameters();
    for (unsigned i = 0; i < n && i < fH1Marginalized.size(); ++i) {
        if (MarginalizedHistogramExists(i))
            I.push_back(i);
        else {
            for (unsigned j = 0; j < n && j < fH2Marginalized[i].size(); ++j)
                if (i != j && MarginalizedHistogramExists(i, j)) {
                    I.push_back(i);
                    break;
                }
        }
    }

    if (I.empty()) {
        BCLog::OutError("BCEngineMCMC::PrintCorrelationPlot : Marginalized distributions not stored. Cannot print correlation plots.");
        return false;
    }

    TCanvas c("c_correlation_plot");
    c.cd();

    double margin = 0.1;
    double padsize = (1 - 2 * margin) / I.size();

    // array with pads holding the histograms
    std::vector<std::vector<TPad*> > pad (I.size(), std::vector<TPad*>(I.size(), NULL));

    // position of pads
    double xlow, xup, ylow, yup;
    double margintop    = 0.01;
    double marginbottom = margintop;
    double marginleft   = margintop;
    double marginright  = marginleft;

    TLatex ylabel;
    ylabel.SetTextSize(5e-2 / I.size());
    ylabel.SetTextAlign(22);     // set to 32, if latex names too long
    ylabel.SetNDC();
    ylabel.SetTextAngle(90);     // set to 80, if latex names too long

    TLatex xlabel;
    xlabel.SetTextSize(5e-2 / I.size());
    xlabel.SetTextAlign(22);     // set to 12, if latex names too long
    xlabel.SetNDC();
    xlabel.SetTextAngle(0);      // set to 350, if latex names too long

    // Box + Text for empty squares:
    TText text_na;
    text_na.SetTextAlign(22);
    text_na.SetTextSize(8e-1 / I.size());
    text_na.SetTextColor(kGray);

    // store histograms for persistency for ROOT until drawing is saved
    std::vector<BCHistogramBase*> bh;

    // drawing all histograms
    for (unsigned i = 0; i < I.size(); ++i) {
        xlow = i * padsize + margin;
        xup = xlow + padsize;

        for (unsigned j = 0; j < I.size(); ++j) {
            yup = 1. - j * padsize - margin;
            ylow = yup - padsize;

            // preparing the pad
            {
                BCAux::RootSideEffectGuard g;
                pad[i][j] =  new TPad(Form("pad_correlation_plots_%d_%d", i, j), "", xlow, ylow, xup, yup);
                // despite the guard, the pads are deleted before the trash
                // cleans up, which leads to a double-delete problem
                // fObjectTrash.Put(pad[i][j]);
            }
            pad[i][j]->SetMargin(marginleft, marginright, marginbottom, margintop);
            pad[i][j]->SetFillColor(kWhite);
            pad[i][j]->Draw();
            pad[i][j]->cd();

            if (i == j)
                bh.push_back(new BCH1D(GetMarginalized(I[i])));
            else {
                if (MarginalizedHistogramExists(I[i], I[j]))
                    bh.push_back(new BCH2D(GetMarginalized(I[i], I[j])));
                else
                    bh.push_back(NULL);
            }
            // bh = MarginalizedHistogramExists(I[i], I[j]) ? GetMarginalized(I[i], I[j]) : NULL;

            if (bh.back()) {

                bh.back()->GetHistogram()->GetXaxis()->SetLabelSize(0);
                bh.back()->GetHistogram()->GetYaxis()->SetLabelSize(0);
                bh.back()->GetHistogram()->GetXaxis()->SetTitleSize(0);
                bh.back()->GetHistogram()->GetYaxis()->SetTitleSize(0);

                if (bh.back()->GetHistogram()->GetDimension() == 1)
                    bh.back()->CopyOptions(fBCH1DdrawingOptions);
                else if (bh.back()->GetHistogram()->GetDimension() == 2)
                    bh.back()->CopyOptions(fBCH2DdrawingOptions);
                bh.back()->SetDrawLegend(false);
                bh.back()->SetStats(false);
                bh.back()->Draw();

            } else if (!(MarginalizedHistogramExists(I[j], I[i])) && I[i] >= I[j]) { // if the histogram is not available, draw "N/A"
                // pad[i][j]->SetFillColor(kWhite);
                text_na.DrawText(.5, .5, "N/A");
            }

            c.cd();

            if (i == 0)               // y-axis labels
                ylabel.DrawLatex(margin * (1 - 8 * ylabel.GetTextSize()), yup - padsize / 2., GetVariable(I[j]).GetLatexNameWithUnits().data());
            if (j == I.size() - 1)        // x-axis labels
                xlabel.DrawLatex(xlow + padsize / 2., margin * (1 - 8 * xlabel.GetTextSize()), GetVariable(I[i]).GetLatexNameWithUnits().data());
        }
    }

    // gPad->RedrawAxis();
    c.Print(filename.data());

    // delete BCHistogramBase objects
    for (unsigned i = 0; i < bh.size(); ++i)
        delete bh[i];

    return true;
}

// ---------------------------------------------------------
bool BCEngineMCMC::PrintParameterLatex(const std::string& filename) const
{
    // open file
    std::ofstream ofi(filename.data());
    ofi.precision(3);

    // check if file is open
    if (!ofi.is_open()) {
        BCLog::OutError("BCEngineMCMC::PrintParameterLatex : Couldn't open file " + filename + ".");
        return false;
    }

    const char* blank = "---";
    unsigned texwidth = 15;

    // print table
    ofi << "\\documentclass[11pt, a4paper]{article}\n\n"
        << "\\usepackage[landscape]{geometry}\n\n"
        << "\\begin{document}\n\n"
        << "  \\begin{table}[ht!]\n\n"
        << "    \\begin{center}\n\n"
        << "      \\begin{tabular}{llllllll}\n\n"
        << "        Parameter &\n"
        << Form("        %*s & %*s & %*s & %*s & %*s & %*s & %*s\\\\\n",
                texwidth, "Mean",
                texwidth, "RMS",
                texwidth, "Gl. Mode",
                texwidth, "Mode",
                texwidth, "Median",
                texwidth, "16\\% quant.",
                texwidth, "84\\% quant.")
        << "        \\hline\n" << std::endl;

    for (unsigned i = 0; i < GetNVariables(); ++i) {

        if (i == GetNParameters())    // first user-defined observable
            ofi << "        &\\\\\n"
                << "        User-defined Observable &\n"
                << Form("        %*s & %*s & %*s & %*s & %*s & %*s & %*s\\\\\n",
                        texwidth, "Mean",
                        texwidth, "RMS",
                        texwidth, "Gl. Mode",
                        texwidth, "Mode",
                        texwidth, "Median",
                        texwidth, "16\\% quant.",
                        texwidth, "84\\% quant.")
                << "        \\hline\n" << std::endl;

        // formate LaTeX name for LaTeX ('#'->'\')
        std::string latexName(GetVariable(i).GetLatexNameWithUnits());
        std::replace(latexName.begin(), latexName.end(), '#', '\\');

        ofi << "        \\ensuremath{" << latexName << "} &" << std::endl;

        unsigned prec = GetVariable(i).GetPrecision();

        // fixed parameter
        if (i < GetNParameters() && GetParameter(i).Fixed())
            ofi << Form("        \\multicolumn{7}{c}{--- fixed to %*.*g ---}\\\\\n", texwidth, prec, GetParameter(i).GetFixedValue()) << std::endl;

        // not fixed
        else {
            BCH1D bch1d = GetMarginalized(i);

            // marginalization exists
            if (bch1d.Valid())
                ofi << Form("        %*.*g & %*.*g & %*.*g & %*.*g & %*.*g & %*.*g & %*.*g\\\\\n",
                            texwidth, prec, bch1d.GetHistogram()->GetMean(),
                            texwidth, prec, bch1d.GetHistogram()->GetRMS(),
                            texwidth, prec, GetBestFitParameters()[i],
                            texwidth, prec, bch1d.GetLocalMode(0),
                            texwidth, prec, bch1d.GetMedian(),
                            texwidth, prec, bch1d.GetQuantile(0.16),
                            texwidth, prec, bch1d.GetQuantile(0.84))
                    << std::endl;

            // marginalization does not exist
            else
                ofi << Form("        %*s & %*s & %*.*g & %*s & %*s & %*s & %*s\\\\\n",
                            texwidth, blank,
                            texwidth, blank,
                            texwidth, prec, GetBestFitParameters()[i],
                            texwidth, blank,
                            texwidth, blank,
                            texwidth, blank,
                            texwidth, blank)
                    << std::endl;
        }
    }

    ofi << "        \\hline\n" << std::endl
        << "      \\end{tabular}\n" << std::endl
        << "      \\caption{Summary of the parameter estimates.}\n" << std::endl
        << "    \\end{center}\n" << std::endl
        << "  \\end{table}" << std::endl
        << std::endl
        << "\\end{document}" << std::endl;

    // close file
    ofi.close();

    // no error
    return true;
}

// ---------------------------------------------------------
unsigned BCEngineMCMC::UpdateFrequency(unsigned N) const
{
    int n = N / 10;
    if (n < 100)
        return 100;
    if (n < 10000)
        return 1000;
    if (n < 100000)
        return 10000;
    return 100000;
}

// ---------------------------------------------------------
BCEngineMCMC::ThreadLocalStorage::ThreadLocalStorage(unsigned dim) :
    parameters(dim, 0.0),
    log_prior(-std::numeric_limits<double>::infinity()),
    log_likelihood(-std::numeric_limits<double>::infinity()),
    log_probability(-std::numeric_limits<double>::infinity()),
    rng(new TRandom3(0)),
    yLocal(dim)
{
}

// ---------------------------------------------------------
BCEngineMCMC::ThreadLocalStorage::ThreadLocalStorage(const ThreadLocalStorage& other) :
    parameters(other.parameters),
    log_prior(other.log_prior),
    log_likelihood(other.log_likelihood),
    log_probability(other.log_probability),
    rng(new TRandom3(*other.rng)),
    yLocal(other.yLocal)
{
}

// ---------------------------------------------------------
BCEngineMCMC::ThreadLocalStorage& BCEngineMCMC::ThreadLocalStorage::operator = (ThreadLocalStorage other)
{
    swap(*this, other);
    return *this;
}

void BCEngineMCMC::ThreadLocalStorage::swap(BCEngineMCMC::ThreadLocalStorage& A, BCEngineMCMC::ThreadLocalStorage& B)
{
    std::swap(A.parameters, B.parameters);
    std::swap(A.log_prior, B.log_prior);
    std::swap(A.log_likelihood, B.log_likelihood);
    std::swap(A.log_probability, B.log_probability);
    std::swap(A.rng, B.rng);
    std::swap(A.yLocal, B.yLocal);
}

// ---------------------------------------------------------
BCEngineMCMC::ThreadLocalStorage::~ThreadLocalStorage()
{
    delete rng;
}

// ---------------------------------------------------------
double BCEngineMCMC::ThreadLocalStorage::scale(double dof)
{
    // when Z is normally distributed with expected value 0 and std deviation sigma
    // and  V is chi-squared distributed with dof degrees of freedom
    // and  Z and V are independent
    // then Z*sqrt(dof/V) is t-distributed with dof degrees of freedom and std deviation sigma
    if (dof <= 0)
        return 1;
    const double chi2 = BCMath::Random::Chi2(rng, dof);
    return sqrt(dof / chi2);
}

// ---------------------------------------------------------
void BCEngineMCMC::SyncThreadStorage()
{
    if (fMCMCNChains > fMCMCThreadLocalStorage.size())
        fRandom.Rndm();         // fix return value of GetSeed()

    // add storage until equal to number of chains
    fMCMCThreadLocalStorage.reserve(fMCMCNChains);
    while (fMCMCThreadLocalStorage.size() < fMCMCNChains)
        fMCMCThreadLocalStorage.push_back(ThreadLocalStorage(GetNParameters()));

    // remove storage until equal to number of chains
    while (fMCMCThreadLocalStorage.size() > fMCMCNChains)
        fMCMCThreadLocalStorage.pop_back();

    // reserve memory to map each thread to a chain. Initially all zero to
    // handle calling LogEval in a serial context.
    int nthreads = 1;
#if THREAD_PARALLELIZATION
    nthreads = omp_get_max_threads();
#endif
    fChainIndex.assign(nthreads, 0);

    // update for each chain
    for (unsigned i = 0; i < fMCMCThreadLocalStorage.size(); ++i) {
        // need full number of parameters, this is passed into user function
        fMCMCThreadLocalStorage[i].parameters.assign(GetNParameters(), 0.0);
        // need only free parameters, these ones are transformed by Cholesky
        fMCMCThreadLocalStorage[i].yLocal.ResizeTo(GetNFreeParameters());

        fMCMCThreadLocalStorage[i].log_prior = -std::numeric_limits<double>::infinity();
        fMCMCThreadLocalStorage[i].log_likelihood = -std::numeric_limits<double>::infinity();
        fMCMCThreadLocalStorage[i].log_probability = -std::numeric_limits<double>::infinity();

        // each chains gets a different seed. fRandom always returns same seed after the fixing done above
        fMCMCThreadLocalStorage[i].rng->SetSeed(fRandom.GetSeed() + i);
        fMCMCThreadLocalStorage[i].rng->Rndm();
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::UpdateChainIndex(int chain)
{
    int id = 0;
#if THREAD_PARALLELIZATION
    id = omp_get_thread_num();
#endif
    fChainIndex.at(id) = chain;
}

// ---------------------------------------------------------
BCEngineMCMC::Statistics::Statistics(unsigned n_par, unsigned n_obs) :
    n_samples(0),
    mean(n_par + n_obs, 0),
    variance(mean.size(), 0),
    stderrpar(n_par, 0),
    stderrobs(n_obs, 0),
    covariance(mean.size(), std::vector<double>(mean.size(), 0)),
    minimum(mean.size(), +std::numeric_limits<double>::infinity()),
    maximum(mean.size(), -std::numeric_limits<double>::infinity()),
    probability_mean(0),
    probability_variance(0),
    modepar(n_par, 0),
    modeobs(n_obs, 0),
    probability_at_mode(-std::numeric_limits<double>::infinity()),
    n_samples_efficiency(0),
    efficiency(n_par, 0.)
{
}

// ---------------------------------------------------------
void BCEngineMCMC::Statistics::Clear(bool clear_mode, bool clear_efficiency)
{
    n_samples = 0;
    mean.clear();
    variance.clear();
    stderrpar.clear();
    stderrobs.clear();
    covariance.clear();
    minimum.clear();
    maximum.clear();
    probability_mean = 0;
    probability_variance = 0;
    if (clear_mode) {
        probability_at_mode = -std::numeric_limits<double>::infinity();
        modepar.clear();
        modeobs.clear();
    }
    if (clear_efficiency) {
        n_samples_efficiency = 0;
        efficiency.clear();
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::Statistics::Init(unsigned n_par, unsigned n_obs)
{
    n_samples = 0;
    mean.assign(n_par + n_obs, 0);
    variance.assign(mean.size(), 0);
    stderrpar.assign(n_par, 0);
    stderrobs.assign(n_obs, 0);
    covariance.assign(mean.size(), std::vector<double>(mean.size(), 0));
    minimum.assign(mean.size(), +std::numeric_limits<double>::infinity());
    maximum.assign(mean.size(), -std::numeric_limits<double>::infinity());
    probability_mean = 0;
    probability_variance = 0;
    probability_at_mode = -std::numeric_limits<double>::infinity();
    modepar.assign(n_par, 0);
    modeobs.assign(n_obs, 0);
    n_samples_efficiency = 0;
    efficiency.assign(n_par, 0.);
}

// ---------------------------------------------------------
void BCEngineMCMC::Statistics::Reset(bool reset_mode, bool reset_efficiency)
{
    n_samples = 0;
    mean.assign(mean.size(), 0);
    variance.assign(variance.size(), 0);
    stderrpar.assign(stderrpar.size(), 0);
    stderrobs.assign(stderrobs.size(), 0);
    covariance.assign(covariance.size(), std::vector<double>(covariance.front().size(), 0));
    minimum.assign(minimum.size(), +std::numeric_limits<double>::infinity());
    maximum.assign(maximum.size(), -std::numeric_limits<double>::infinity());
    probability_mean = 0;
    probability_variance = 0;
    if (reset_mode) {
        probability_at_mode = -std::numeric_limits<double>::infinity();
        modepar.assign(modepar.size(), 0);
        modeobs.assign(modeobs.size(), 0);
    }
    if (reset_efficiency) {
        efficiency.assign(efficiency.size(), 0);
        n_samples_efficiency = 0;
    }
}

// ---------------------------------------------------------
void BCEngineMCMC::Statistics::ResetEfficiencies()
{
    efficiency.assign(efficiency.size(), 0);
    n_samples_efficiency = 0;
}

// ---------------------------------------------------------
void BCEngineMCMC::Statistics::Update(const ChainState& cs)
{
    if (mean.size() != cs.parameters.size() + cs.observables.size())
        return;

    // increment number of samples
    ++n_samples;

    // check mode
    if (cs.log_probability > probability_at_mode) {
        modepar = cs.parameters;
        modeobs = cs.observables;
        probability_at_mode = cs.log_probability;
    }

    // update probability mean and variance
    double prob_delta = cs.log_probability - probability_mean;
    probability_mean += prob_delta / n_samples;
    probability_variance += (n_samples > 1) ? prob_delta * prob_delta / n_samples - probability_variance / (n_samples - 1) : 0;

    // update parameter means and (co)variances, and maxima and minima:

    // vector to store difference from current mean
    std::vector<double> delta(mean.size(), 0);

    // loop over values
    for (unsigned i = 0; i < mean.size(); ++i) {
        // get value
        double x = (i < cs.parameters.size()) ? cs.parameters[i] : cs.observables[i - cs.parameters.size()];
        // store difference to current mean
        delta[i] = x - mean[i];
        // update mean
        mean[i] += delta[i] / n_samples;
        // update variance
        variance[i] += (n_samples > 1) ? delta[i] * delta[i] / n_samples - variance[i] / (n_samples - 1.) : 0;
        // update minimum & maximum
        minimum[i] = std::min(minimum[i], x);
        maximum[i] = std::max(maximum[i], x);
    }
    // update covariances
    if (n_samples > 1) {
        for (unsigned i = 0; i < mean.size(); ++i)
            for (unsigned j = i; j < mean.size(); ++j)
                covariance[i][j] += delta[i] * delta[j] / n_samples - covariance[i][j] / (n_samples - 1);
    }
    // update stderrors
    for (unsigned i = 0; i < modepar.size(); ++i)
        stderrpar[i] = std::sqrt(variance[i]);
    for (unsigned i = 0; i < modeobs.size(); ++i)
        stderrobs[i] = std::sqrt(variance[stderrpar.size() + i]);
}

// ---------------------------------------------------------
BCEngineMCMC::Statistics& BCEngineMCMC::Statistics::operator+=(const BCEngineMCMC::Statistics& rhs)
{
    // if rhs is empty
    if (rhs.n_samples == 0)
        return *this;

    // if this is empty
    if (n_samples == 0)
        return operator=(rhs);

    if (mean.size() != rhs.mean.size())
        return *this;

    // check mode:
    if (rhs.probability_at_mode > probability_at_mode) {
        probability_at_mode = rhs.probability_at_mode;
        modepar = rhs.modepar;
        modeobs = rhs.modeobs;
    }

    double n = n_samples + rhs.n_samples;
    double N = (n > 0) ? n_samples * rhs.n_samples / n : 0;

    probability_variance = (probability_variance * (n_samples - 1) + rhs.probability_variance * (rhs.n_samples - 1) + (rhs.probability_mean - probability_mean) * (rhs.probability_mean - probability_mean) * N) / (n - 1);
    probability_mean = (n_samples * probability_mean + rhs.n_samples * rhs.probability_mean) / n;

    // parameter variables:
    for (unsigned i = 0; i < mean.size(); ++i) {
        // check minimum
        if (rhs.minimum[i] < minimum[i])
            minimum[i] = rhs.minimum[i];
        // check maximum
        if (rhs.maximum[i] > maximum[i])
            maximum[i] = rhs.maximum[i];
        // combine variances
        variance[i] = (variance[i] * (n_samples - 1) + rhs.variance[i] * (rhs.n_samples - 1) + (rhs.mean[i] - mean[i]) * (rhs.mean[i] - mean[i]) * N) / (n - 1);
        // combine covariances
        for (unsigned j = 0; j < covariance[i].size(); ++j)
            covariance[i][j] = (covariance[i][j] * (n_samples - 1) + rhs.covariance[i][j] * (rhs.n_samples - 1) + (rhs.mean[i] - mean[i]) * (rhs.mean[j] - mean[j]) * N) / (n - 1);
        // combine means
        mean[i] = (n_samples * mean[i] + rhs.n_samples * rhs.mean[i]) / n;
    }
    // update stderrors
    for (unsigned i = 0; i < stderrpar.size(); ++i)
        stderrpar[i] = std::sqrt(variance[i]);
    for (unsigned i = 0; i < stderrobs.size(); ++i)
        stderrobs[i] = std::sqrt(variance[stderrpar.size() + i]);

    // combine n_samples
    n_samples = n;

    // combine efficiencies
    double n_eff = n_samples_efficiency + rhs.n_samples_efficiency;
    if (n_eff > 0)
        for (unsigned i = 0; i < efficiency.size(); ++i)
            efficiency[i] = (n_samples_efficiency * efficiency[i] + rhs.n_samples_efficiency * rhs.efficiency[i]) / (n_eff);

    // combine efficiency samples
    n_samples_efficiency = n_eff;

    return *this;
}
