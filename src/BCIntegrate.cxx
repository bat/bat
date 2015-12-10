/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------
#include <config.h>

#include "BCIntegrate.h"
#include "BCLog.h"
#include "BCMath.h"
#include "BCParameter.h"
#include "BCH2D.h"
#include "BCH1D.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#ifdef HAVE_CUBA_H
#include <cuba.h>
#endif

#include <algorithm>
#include <limits>
#include <math.h>

namespace
{
/**
 * Hold an instance of BCIntegrate to emulate a global variable for use with Minuit
 */
class BCIntegrateHolder
{
private:
    BCIntegrate* global_this;

public:

    BCIntegrateHolder() :
        global_this(NULL)
    {
    }

    /**
     * Set and/or retrieve the static BCIntegrate object
     */
    static BCIntegrate* instance(BCIntegrate* obj = NULL)
    {
        static BCIntegrateHolder result;
        if (obj)
            result.global_this = obj;

        return result.global_this;
    }
};
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(const std::string& name) :
    BCEngineMCMC(name),
    fMinuit(0),
    fMinuitErrorFlag(0),
    fFlagIgnorePrevOptimization(false),
    fSAT0(100),
    fSATmin(0.1),
    fSATree(0),
    fFlagWriteSAToFile(false),
    fSANIterations(0),
    fSATemperature(0),
    fSALogProb(0),
    fFlagMarginalized(false),
    fSAOutputFile(0),
    fSAOutputFilename(""),
    fSAOutputFileOption(""),
    fSAOutputFileAutoclose(false),
    fOptimizationMethodCurrent(BCIntegrate::kOptDefault),
    fOptimizationMethodUsed(BCIntegrate::kOptEmpty),
    fIntegrationMethodCurrent(BCIntegrate::kIntDefault),
    fIntegrationMethodUsed(BCIntegrate::kIntEmpty),
    fMarginalizationMethodCurrent(BCIntegrate::kMargDefault),
    fMarginalizationMethodUsed(BCIntegrate::kMargEmpty),
    fSASchedule(BCIntegrate::kSACauchy),
    fNIterationsMin(0),
    fNIterationsMax(1000000),
    fNIterationsPrecisionCheck(1000),
    fNIterations(0),
    fLogMaximum(-std::numeric_limits<double>::infinity()),
    fIntegral(-1),
    fRelativePrecision(1e-2),
    fAbsolutePrecision(1e-6),
    fError(-999.),
    fCubaIntegrationMethod(BCIntegrate::kCubaDefault)
{
    fMinuitArglist[0] = 20000;
    fMinuitArglist[1] = 0.01;
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(const std::string& filename, const std::string& name, bool loadObservables) :
    BCEngineMCMC(filename, name, loadObservables),
    fMinuit(0),
    fMinuitErrorFlag(0),
    fFlagIgnorePrevOptimization(false),
    fSAT0(100),
    fSATmin(0.1),
    fSATree(0),
    fFlagWriteSAToFile(false),
    fSANIterations(0),
    fSATemperature(0),
    fSALogProb(0),
    fFlagMarginalized(false),
    fSAOutputFile(0),
    fSAOutputFilename(""),
    fSAOutputFileOption(""),
    fSAOutputFileAutoclose(false),
    fOptimizationMethodCurrent(BCIntegrate::kOptDefault),
    fOptimizationMethodUsed(BCIntegrate::kOptEmpty),
    fIntegrationMethodCurrent(BCIntegrate::kIntDefault),
    fIntegrationMethodUsed(BCIntegrate::kIntEmpty),
    fMarginalizationMethodCurrent(BCIntegrate::kMargDefault),
    fMarginalizationMethodUsed(BCIntegrate::kMargEmpty),
    fSASchedule(BCIntegrate::kSACauchy),
    fNIterationsMin(0),
    fNIterationsMax(1000000),
    fNIterationsPrecisionCheck(1000),
    fNIterations(0),
    fLogMaximum(-std::numeric_limits<double>::infinity()),
    fIntegral(-1),
    fRelativePrecision(1e-2),
    fAbsolutePrecision(1e-6),
    fError(-999.),
    fCubaIntegrationMethod(BCIntegrate::kCubaDefault)
{
    fMinuitArglist[0] = 20000;
    fMinuitArglist[1] = 0.01;
}

// ---------------------------------------------------------
BCIntegrate::BCIntegrate(const BCIntegrate& other)
    : BCEngineMCMC(other),
      fMinuit(new TMinuit()),
      fMinuitErrorFlag(other.fMinuitErrorFlag),
      fFlagIgnorePrevOptimization(other.fFlagIgnorePrevOptimization),
      fSAT0(other.fSAT0),
      fSATmin(other.fSATmin),
      fSATree(0),
      fFlagWriteSAToFile(other.fFlagWriteSAToFile),
      fSANIterations(other.fSANIterations),
      fSATemperature(other.fSATemperature),
      fSALogProb(other.fSALogProb),
      fSAx(other.fSAx),
      fFlagMarginalized(other.fFlagMarginalized),
      fSAOutputFile(0),
      fSAOutputFilename(other.fSAOutputFilename),
      fSAOutputFileOption(other.fSAOutputFileOption),
      fSAOutputFileAutoclose(other.fSAOutputFileAutoclose),
      fOptimizationMethodCurrent(other.fOptimizationMethodCurrent),
      fOptimizationMethodUsed(other.fOptimizationMethodUsed),
      fIntegrationMethodCurrent(other.fIntegrationMethodCurrent),
      fIntegrationMethodUsed(other.fIntegrationMethodUsed),
      fMarginalizationMethodCurrent(other.fMarginalizationMethodCurrent),
      fMarginalizationMethodUsed(other.fMarginalizationMethodUsed),
      fSASchedule(other.fSASchedule),
      fNIterationsMin(other.fNIterationsMin),
      fNIterationsMax(other.fNIterationsMax),
      fNIterationsPrecisionCheck(other.fNIterationsPrecisionCheck),
      fNIterations(other.fNIterations),
      fBestFitParameters(other.fBestFitParameters),
      fBestFitParameterErrors(other.fBestFitParameterErrors),
      fLogMaximum(other.fLogMaximum),
      fIntegral(other.fIntegral),
      fRelativePrecision(other.fRelativePrecision),
      fAbsolutePrecision(other.fAbsolutePrecision),
      fError(other.fError),
      fCubaIntegrationMethod(other.fCubaIntegrationMethod),
      fCubaVegasOptions(other.fCubaVegasOptions),
      fCubaSuaveOptions(other.fCubaSuaveOptions),
      fCubaDivonneOptions(other.fCubaDivonneOptions),
      fCubaCuhreOptions(other.fCubaCuhreOptions)
{
    fMinuitArglist[0] = other.fMinuitArglist[0];
    fMinuitArglist[1] = other.fMinuitArglist[1];
}

// ---------------------------------------------------------
BCIntegrate::~BCIntegrate()
{
    delete fMinuit;
}

// ---------------------------------------------------------
void swap(BCIntegrate& A, BCIntegrate& B)
{
    swap(static_cast<BCEngineMCMC&>(A), static_cast<BCEngineMCMC&>(B));
    std::swap(A.fMinuit, B.fMinuit);
    std::swap(A.fMinuitArglist, B.fMinuitArglist);
    std::swap(A.fMinuitErrorFlag, B.fMinuitErrorFlag);
    std::swap(A.fFlagIgnorePrevOptimization, B.fFlagIgnorePrevOptimization);
    std::swap(A.fSAT0, B.fSAT0);
    std::swap(A.fSATmin, B.fSATmin);
    std::swap(A.fSATree, B.fSATree);
    std::swap(A.fFlagWriteSAToFile, B.fFlagWriteSAToFile);
    std::swap(A.fSANIterations, B.fSANIterations);
    std::swap(A.fSATemperature, B.fSATemperature);
    std::swap(A.fSALogProb, B.fSALogProb);
    std::swap(A.fSAx, B.fSAx);
    std::swap(A.fFlagMarginalized, B.fFlagMarginalized);
    std::swap(A.fSAOutputFile, A.fSAOutputFile);
    std::swap(A.fSAOutputFilename, B.fSAOutputFilename);
    std::swap(A.fSAOutputFileOption, B.fSAOutputFileOption);
    std::swap(A.fSAOutputFileAutoclose, B.fSAOutputFileAutoclose);
    std::swap(A.fOptimizationMethodCurrent, B.fOptimizationMethodCurrent);
    std::swap(A.fOptimizationMethodUsed, B.fOptimizationMethodUsed);
    std::swap(A.fIntegrationMethodCurrent, B.fIntegrationMethodCurrent);
    std::swap(A.fIntegrationMethodUsed, B.fIntegrationMethodUsed);
    std::swap(A.fMarginalizationMethodCurrent, B.fMarginalizationMethodCurrent);
    std::swap(A.fMarginalizationMethodUsed, B.fMarginalizationMethodUsed);
    std::swap(A.fSASchedule, B.fSASchedule);
    std::swap(A.fNIterationsMin, B.fNIterationsMin);
    std::swap(A.fNIterationsMax, B.fNIterationsMax);
    std::swap(A.fNIterationsPrecisionCheck, B.fNIterationsPrecisionCheck);
    std::swap(A.fNIterations, B.fNIterations);
    std::swap(A.fBestFitParameters, B.fBestFitParameters);
    std::swap(A.fBestFitParameterErrors, B.fBestFitParameterErrors);
    std::swap(A.fLogMaximum, B.fLogMaximum);
    std::swap(A.fIntegral, B.fIntegral);
    std::swap(A.fRelativePrecision, B.fRelativePrecision);
    std::swap(A.fAbsolutePrecision, B.fAbsolutePrecision);
    std::swap(A.fError, B.fError);
    std::swap(A.fCubaIntegrationMethod, B.fCubaIntegrationMethod);
    std::swap(A.fCubaVegasOptions, B.fCubaVegasOptions);
    std::swap(A.fCubaSuaveOptions, B.fCubaSuaveOptions);
    std::swap(A.fCubaDivonneOptions, B.fCubaDivonneOptions);
    std::swap(A.fCubaCuhreOptions, B.fCubaCuhreOptions);
}


// ---------------------------------------------------------
const std::vector<double>& BCIntegrate::GetBestFitParameters() const
{
    if (fBestFitParameters.empty() and !BCEngineMCMC::GetBestFitParameters().empty())
        return BCEngineMCMC::GetBestFitParameters();
    return fBestFitParameters;
}

// ---------------------------------------------------------
void BCIntegrate::ResetResults()
{
    BCEngineMCMC::ResetResults();

    fBestFitParameterErrors.clear();
    fBestFitParameters.clear();
    fLogMaximum = -std::numeric_limits<double>::infinity();

    // remove marginalized histograms
    // set marginalization flag
    fFlagMarginalized = false;
}

// ---------------------------------------------------------
double BCIntegrate::Integrate()
{
    // check if parameters are defined
    if (fParameters.Size() < 1) {
        BCLog::OutError("BCIntegrate::Integrate : No parameters defined. Aborting.");
        return -1.;
    }

    // output
    if (!(fIntegrationMethodCurrent == BCIntegrate::kIntDefault) && !(fIntegrationMethodCurrent == BCIntegrate::kIntEmpty))
        BCLog::OutSummary("Integrate using " + DumpCurrentIntegrationMethod());

    switch (fIntegrationMethodCurrent) {

        // Empty
        case BCIntegrate::kIntEmpty: {
            BCLog::OutWarning("BCIntegrate::Integrate : No integration method chosen.");
            return 0;
        }

        // Monte Carlo Integration
        case BCIntegrate::kIntMonteCarlo: {
            std::vector<double> sums (2, 0.0);
            sums.push_back(GetParameters().Volume());
            fIntegral = Integrate(kIntMonteCarlo,
                                  &BCIntegrate::GetRandomVectorInParameterSpace,
                                  &BCIntegrate::EvaluatorMC,
                                  &IntegralUpdaterMC,
                                  sums);

            // set used integration method
            fIntegrationMethodUsed = BCIntegrate::kIntMonteCarlo;

            // return integral
            return fIntegral;
        }

        // CUBA library
        case BCIntegrate::kIntCuba:
            fIntegral = IntegrateCuba();

            // set used integration method
            fIntegrationMethodUsed = BCIntegrate::kIntCuba;

            // return integral
            return fIntegral;

        // CUBA library
        case BCIntegrate::kIntGrid:
            fIntegral = IntegrateSlice();

            // set used integration method
            fIntegrationMethodUsed = BCIntegrate::kIntGrid;

            // return integral
            return fIntegral;

        // default
        case BCIntegrate::kIntDefault: {
#ifdef HAVE_CUBA_H
            SetIntegrationMethod(BCIntegrate::kIntCuba);
#else
            if (GetNFreeParameters() <= 3)
                SetIntegrationMethod(BCIntegrate::kIntGrid);
            else
                SetIntegrationMethod(BCIntegrate::kIntMonteCarlo);
#endif
            return Integrate();
        }

        default:
            BCLog::OutError(Form("BCIntegrate::Integrate : Invalid integration method: %d", fIntegrationMethodCurrent));
            break;
    }

    return 0;
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod intmethod)
{
    // remember original method
    BCIntegrationMethod method_temp = fIntegrationMethodCurrent;

    // set method
    SetIntegrationMethod(intmethod);

    // run algorithm
    double integral = Integrate();

    // re-set original method
    SetIntegrationMethod(method_temp);

    // return integral
    return integral;

}

// ---------------------------------------------------------
void BCIntegrate::SetBestFitParameters(const std::vector<double>& x, const double& new_value, double& old_value)
{
    if (new_value < old_value)
        return;
    old_value = new_value;
    SetBestFitParameters(x);
}

// ---------------------------------------------------------
void BCIntegrate::SetBestFitParameters(const std::vector<double>& x)
{
    fBestFitParameters = x;
    if (fBestFitParameters.size() == GetNParameters() and GetNObservables() > 0) {
        Eval(fBestFitParameters); // in case user uses likelihood to set values needed for observable calculation
        CalculateObservables(fBestFitParameters);
        for (unsigned i = 0; i < GetNObservables(); ++i)
            fBestFitParameters.push_back(GetObservable(i).Value());
    }
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtStartOfIntegration(BCIntegrationMethod type, BCCubaMethod cubatype)
{

    const unsigned NVarNow = GetParameters().GetNFreeParameters();

    BCLog::LogLevel level = BCLog::summary;

    if (GetNParameters() != NVarNow) {

        level = BCLog::detail;
        bool printed = false;

        if (type == kIntCuba) {
            BCLog::OutDetail(Form("Running %s (%s) integration over %i dimensions out of %i.",
                                  DumpIntegrationMethod(type).c_str(),
                                  DumpCubaIntegrationMethod(cubatype).c_str(),
                                  NVarNow, GetNParameters()));
            printed = true;
        }

        if (not printed)
            BCLog::OutDetail(Form("Running %s integration over %i dimensions out of %i.",
                                  DumpIntegrationMethod(type).c_str(),
                                  NVarNow, GetNParameters()));

        BCLog::OutDetail(" --> Fixed parameters:");
        for (unsigned i = 0; i < GetNParameters(); i++)
            if (GetParameter(i).Fixed())
                BCLog::OutDetail(Form("      %3i :  %g", i, GetParameter(i).GetFixedValue()));
    } else {
        bool printed = false;

        if (type == kIntCuba) {
            BCLog::OutDetail(Form("Running %s (%s) integration over %i dimensions.",
                                  DumpIntegrationMethod(type).c_str(),
                                  DumpCubaIntegrationMethod(cubatype).c_str(),
                                  NVarNow));
            printed = true;
        }

        if (not printed)
            BCLog::OutDetail(Form("Running %s integration over %i dimensions.",
                                  DumpIntegrationMethod(type).c_str(),
                                  NVarNow));
    }

    if (GetNIterationsMin() > 0 && GetNIterationsMax() > 0 ) {
        BCLog::Out(level, Form(" --> Minimum number of iterations: %i", GetNIterationsMin()));
        BCLog::Out(level, Form(" --> Maximum number of iterations: %i", GetNIterationsMax()));
    }
    BCLog::Out(level, Form(" --> Target relative precision:    %e", GetRelativePrecision()));
    BCLog::Out(level, Form(" --> Target absolute precision:    %e", GetAbsolutePrecision()));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtEndOfIntegration(double integral, double absprecision, double relprecision, int nIterations)
{
    BCLog::OutSummary(Form(" --> Result of integration:        %e +- %e", integral, absprecision));
    BCLog::OutSummary(Form(" --> Obtained relative precision:  %e. ", relprecision));
    if (nIterations >= 0)
        BCLog::OutSummary(Form(" --> Number of iterations:         %i", nIterations));
}

// ---------------------------------------------------------
void BCIntegrate::LogOutputAtIntegrationStatusUpdate(BCIntegrationMethod type, double integral, double absprecision, int nIterations)
{
    BCLog::OutDetail(Form("%s. Iteration %i, integral: %e +- %e.", DumpIntegrationMethod(type).c_str(), nIterations, integral, absprecision));
}

// ---------------------------------------------------------
double BCIntegrate::Integrate(BCIntegrationMethod type, tRandomizer randomizer, tEvaluator evaluator, tIntegralUpdater updater,
                              std::vector<double>& sums)
{
    LogOutputAtStartOfIntegration(type, NCubaMethods);

    // reset variables
    double pmax = 0.;
    double integral = 0.;
    double absprecision = 2.*fAbsolutePrecision;
    double relprecision = 2.*fRelativePrecision;

    std::vector<double> randx (GetNParameters(), 0.);

    // how often to print out the info line to screen
    int nwrite = UpdateFrequency(fNIterationsMax);

    // reset number of iterations
    fNIterations = 0;
    bool accepted;

    // iterate while number of iterations is lower than minimum number of iterations
    // or precision is not reached and the number of iterations is lower than maximum number of iterations
    while ((GetRelativePrecision() < relprecision and GetAbsolutePrecision() < absprecision and GetNIterations() < GetNIterationsMax())
            or GetNIterations() < GetNIterationsMin()) {

        // get random numbers
        (this->*randomizer)(randx);

        // evaluate function at sampled point
        // updating sums & checking for maximum probability

        SetBestFitParameters(randx, (this->*evaluator)(sums, randx, accepted), pmax);

        // increase number of iterations
        if (accepted)
            ++fNIterations;

        // update precisions
        if (fNIterations % fNIterationsPrecisionCheck == 0) {
            (*updater)(sums, fNIterations, integral, absprecision);
            relprecision = absprecision / integral;
        }

        // write status
        if (fNIterations % nwrite == 0) {
            double temp_integral;
            double temp_absprecision;
            (*updater)(sums, fNIterations, temp_integral, temp_absprecision);
            LogOutputAtIntegrationStatusUpdate(type, temp_integral, temp_absprecision, fNIterations);
        }
    }

    // calculate integral
    (*updater)(sums, fNIterations, integral, absprecision);
    relprecision = absprecision / integral;

    if (fNIterations >= GetNIterationsMax())
        BCLog::OutWarning("BCIntegrate::Integrate: Did not converge within maximum number of iterations");

    // print to log
    LogOutputAtEndOfIntegration(integral, absprecision, relprecision, fNIterations);

    fError = absprecision;
    return integral;
}

// ---------------------------------------------------------
double BCIntegrate::EvaluatorMC(std::vector<double>& sums, const std::vector<double>& point, bool& accepted)
{
    const double value = Eval(point);

    // all samples accepted in sample-mean integration
    accepted = true;

    // add value to sum and sum of squares
    sums[0] += value;
    sums[1] += value * value;

    return value;
}

// ---------------------------------------------------------
void BCIntegrate::IntegralUpdaterMC(const std::vector<double>& sums, const int& nIterations, double& integral, double& absprecision)
{
    // sample mean including the volume of the parameter space
    integral = sums[2] * sums[0] / nIterations;

    // unbiased estimate
    absprecision = sqrt((1.0 / (nIterations - 1)) * (sums[2] * sums[2] * sums[1] / double(nIterations) - integral * integral));
}

// ---------------------------------------------------------
bool BCIntegrate::CheckMarginalizationAvailability(BCMarginalizationMethod type)
{
    switch (type) {
        case BCIntegrate::kMargMonteCarlo:
            return false;
        case BCIntegrate::kMargMetropolis:
            return true;
        case BCIntegrate::kMargGrid:
            return true;
        case BCIntegrate::kMargDefault:
            return true;
        default:
            BCLog::OutError(Form("BCIntegrate::CheckMarginalizationAvailability. Invalid marginalization method: %d.", type));
            break;
    }
    return false;
}

// ---------------------------------------------------------
bool BCIntegrate::CheckMarginalizationIndices(TH1* hist, const std::vector<unsigned>& index)
{
    if (index.empty()) {
        BCLog::OutError("BCIntegrate::Marginalize : No marginalization parameters chosen.");
        return false;
    }

    if (index.size() >= 4 or index.size() > GetNParameters()) {
        BCLog::OutError("BCIntegrate::Marginalize : Too many marginalization parameters.");
        return false;
    }

    if ((int)index.size() < hist->GetDimension()) {
        BCLog::OutError(Form("BCIntegrate::Marginalize : Too few (%d) indices supplied for histogram dimension (%d)", (int)index.size(), hist->GetDimension()));
        return false;
    }

    for (unsigned i = 0; i < index.size(); i++) {
        // check if indices are in bounds
        if (index[i] >= GetNParameters()) {
            BCLog::OutError(Form("BCIntegrate::Marginalize : Parameter index (%d) out of bound.", index[i]));
            return false;
        }
        // check for duplicate indices
        for (unsigned j = 0; j < index.size(); j++)
            if (i != j and index[i] == index[j]) {
                BCLog::OutError(Form("BCIntegrate::Marginalize : Parameter index (%d) appears more than once", index[i]));
                return false;
            }
    }
    return true;
}

// ---------------------------------------------------------
int BCIntegrate::MarginalizeAll()
{

    // check if parameters are defined
    if (GetNParameters() < 1) {
        BCLog::OutError("BCIntegrate::MarginalizeAll : No parameters defined. Aborting.");
        return 0;
    }

    // check if marginalization method is defined
    if (!CheckMarginalizationAvailability(GetMarginalizationMethod())) {
        BCLog::OutError(Form("BCIntegrate::MarginalizeAll : Marginalization method not implemented \'%s\'. Aborting.", DumpCurrentMarginalizationMethod().c_str()));
        return 0;
    }

    // output
    if (!(fMarginalizationMethodCurrent == BCIntegrate::kMargDefault) && !(fMarginalizationMethodCurrent == BCIntegrate::kMargEmpty))
        BCLog::OutSummary(Form("Marginalize using %s", DumpCurrentMarginalizationMethod().c_str()));


    switch (GetMarginalizationMethod()) {

        // Empty
        case BCIntegrate::kMargEmpty: {
            BCLog::OutWarning("BCIntegrate::MarginalizeAll : No marginalization method chosen.");
            return 0;
        }

        // Markov Chain Monte Carlo
        case BCIntegrate::kMargMetropolis: {
            // start preprocess
            MarginalizePreprocess();

            // run the Markov chains
            Metropolis();

            // start postprocess
            MarginalizePostprocess();

            // set used marginalization method
            fMarginalizationMethodUsed = BCIntegrate::kMargMetropolis;

            // check if mode is better than previous one
            if ( (!fFlagIgnorePrevOptimization) && (fLogMaximum < BCEngineMCMC::GetLogMaximum()) ) {
                fBestFitParameters      = BCEngineMCMC::GetBestFitParameters();
                fBestFitParameterErrors.assign(fBestFitParameters.size(), std::numeric_limits<double>::infinity());
                fLogMaximum             = BCEngineMCMC::GetLogMaximum();
                fOptimizationMethodUsed = BCIntegrate::kOptMetropolis;
            }

            break;
        }

        // Sample Mean
        case BCIntegrate::kMargMonteCarlo: {
            return 0;
        }

        // Slice
        case BCIntegrate::kMargGrid: {
            if (GetParameters().GetNFreeParameters() > 2) {
                BCLog::OutWarning("BCIntegrate::MarginalizeAll : Grid marginalization only works for 1D and 2D problems.");
                return 0;
            }
            if (GetNObservables() > 0)
                BCLog::OutWarning("BCIntegrate::MarginalizeAll : Grid marginalization will not store resutls for user-defined observables.");

            // start preprocess
            MarginalizePreprocess();

            std::vector<double> fixpoint = GetParameters().GetFixedValues();

            // delete pre-existing marginalizations
            for (unsigned i = 0; i < fH1Marginalized.size(); ++i)
                if (fH1Marginalized[i])
                    delete fH1Marginalized[i];
            for (unsigned i = 0; i < fH2Marginalized.size(); ++i)
                for (unsigned j = 0; j < fH2Marginalized[i].size(); ++j)
                    if (fH2Marginalized[i][j])
                        delete fH2Marginalized[i][j];
            // set all pointers to zero
            fH1Marginalized.assign(GetNParameters(), NULL);
            fH2Marginalized.assign(GetNParameters(), std::vector<TH2*>(GetNParameters(), NULL));

            // count how often posterior is evaluated
            unsigned nIterations = 0;

            // store highest probability
            double log_max_val = -std::numeric_limits<double>::infinity();
            std::vector<double> bestfit_parameters = fixpoint;
            std::vector<double> bestfit_errors(GetNVariables(), std::numeric_limits<double>::infinity());

            if (GetParameters().GetNFreeParameters() == 1) { // Marginalize the free parameter to a 1D histogram
                for (unsigned i = 0; i < GetNParameters(); ++i) {
                    if (GetParameter(i).Fixed())
                        continue;

                    // calculate slice
                    fH1Marginalized[i] = GetSlice(i, nIterations, log_max_val, fixpoint, 0, false);

                    if ( fH1Marginalized[i]->Integral() > 0 ) { // histogram is nonempty
                        const int index = fH1Marginalized[i]->GetMaximumBin();
                        bestfit_parameters[i] = fH1Marginalized[i]->GetBinCenter(index);
                        // approximate error by flat distribution in bin
                        bestfit_errors[i] = fH1Marginalized[i]->GetBinWidth(index) / sqrt(12);
                    }
                }
            }

            else if (GetNFreeParameters() == 2) { // marginalize the two free parameters to a 2D histogram
                for (unsigned i = 0; i < GetNParameters(); ++i)
                    for (unsigned j = i + 1; j < GetNParameters(); ++j) {
                        if (GetParameter(i).Fixed() or GetParameter(j).Fixed())
                            continue;

                        // calculate slice
                        fH2Marginalized[i][j] = GetSlice(i, j, nIterations, log_max_val, fixpoint, 0, false);

                        if ( fH2Marginalized[i][j]->Integral() > 0 ) { // histogram is nonempty
                            int index1, index2, useless_index;
                            fH2Marginalized[i][j]->GetMaximumBin(index1, index2, useless_index);
                            bestfit_parameters[i] = fH2Marginalized[i][j]->GetXaxis()->GetBinCenter(index1);
                            bestfit_parameters[j] = fH2Marginalized[i][j]->GetYaxis()->GetBinCenter(index2);
                            // approximate errors by flat distribution in bin
                            bestfit_errors[i] = fH2Marginalized[i][j]->GetXaxis()->GetBinWidth(index1) / sqrt(12);
                            bestfit_errors[j] = fH2Marginalized[i][j]->GetYaxis()->GetBinWidth(index2) / sqrt(12);

                            // 1D marginalize by projecting
                            fH1Marginalized[i] = fH2Marginalized[i][j]->ProjectionX(Form("h1_%s_parameter_%i", GetSafeName().data() , i));
                            fH1Marginalized[j] = fH2Marginalized[i][j]->ProjectionY(Form("h1_%s_parameter_%i", GetSafeName().data() , j));
                        }
                    }
            }

            // save if improved the log posterior

            if (fBestFitParameters.empty() or log_max_val > GetLogMaximum()) {
                fLogMaximum = log_max_val;
                SetBestFitParameters(bestfit_parameters);
                fBestFitParameterErrors = bestfit_errors;
            }

            BCLog::OutDetail(" --> Global mode from Grid:");
            BCLog::OutDebug(Form(" --> Posterior value: %g", log_max_val));
            PrintParameters(GetBestFitParameters(), BCLog::OutDetail);

            // start postprocess
            MarginalizePostprocess();

            // set used marginalization method
            fMarginalizationMethodUsed = BCIntegrate::kMargGrid;

            break;
        }

        // default
        case BCIntegrate::kMargDefault: {
            if ( GetNFreeParameters() <= 2 and GetNObservables() == 0)
                SetMarginalizationMethod(BCIntegrate::kMargGrid);
            else
                SetMarginalizationMethod(BCIntegrate::kMargMetropolis);

            // call again
            return MarginalizeAll();

            break;
        }

        default: {
            BCLog::OutError(Form("BCIntegrate::MarginalizeAll : Invalid marginalization method: %d", GetMarginalizationMethod()));
            return 0;
            break;
        }

    }	// end switch

    // set flag
    fFlagMarginalized = true;

    return 1;
}

// ---------------------------------------------------------
int BCIntegrate::MarginalizeAll(BCIntegrate::BCMarginalizationMethod margmethod)
{
    // remember original method
    BCMarginalizationMethod method_temp = fMarginalizationMethodCurrent;

    // set method
    SetMarginalizationMethod(margmethod);

    // run algorithm
    int result = MarginalizeAll();

    // re-set original method
    SetMarginalizationMethod(method_temp);

    // return result
    return result;
}

// ---------------------------------------------------------
TH1* BCIntegrate::GetSlice(std::vector<unsigned> indices, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters, int nbins, bool normalize)
{
    if (indices.empty()) {
        BCLog::OutError("BCIntegrate::GetSlice : No parameter indices provided.");
        return NULL;
    }

    if (indices.size() > 3) {
        BCLog::OutError("BCIntegrate::GetSlice : Too many parameter indices provided. Max is three.");
        return NULL;
    }

    for (unsigned i = 0; i < indices.size(); ++i)
        if (indices[i] >= GetNParameters()) {
            BCLog::OutError("BCIntegrate::GetSlice : Parameter index out of range.");
            return NULL;
        } else // check for duplicates
            for (unsigned j = i + 1; j < indices.size(); ++j)
                if (indices[i] == indices[j]) {
                    BCLog::OutError("BCIntegrate::GetSlice : duplicate parameter indices provided.");
                    return NULL;
                }

    // project (N+n)D slice for (N)D slice
    if (parameters.empty()  and indices.size() < GetNFreeParameters() and GetNFreeParameters() <= 3) {
        unsigned N = indices.size();

        // find unfixed parameters and add them to indices
        for (unsigned i = 0; i < GetNParameters(); ++i)
            if (std::find(indices.begin(), indices.end(), i) == indices.end() and !GetParameter(i).Fixed())
                indices.push_back(i);

        // slice in full free dimensions
        TH1* slice_NnD = GetSlice(indices, nIterations, log_max_val, parameters, nbins, normalize);
        if (slice_NnD == NULL)
            return NULL;

        if (N == 1) {
            TH1* h = NULL;
            // 2->1
            if (dynamic_cast<TH2*>(slice_NnD) != NULL)
                h = dynamic_cast<TH2*>(slice_NnD)->ProjectionX(Form("h1_slice_%s_%d", GetSafeName().data(), indices[0]));
            // 3->1
            if (dynamic_cast<TH3*>(slice_NnD) != NULL)
                h = dynamic_cast<TH3*>(slice_NnD)->ProjectionX(Form("h1_slice_%s_%d", GetSafeName().data(), indices[0]));
            if (h)
                h->SetYTitle(Form("p(%s | data)", GetParameter(indices[0]).GetLatexName().data()));
            return h;
        }
        // 3->2
        if (N == 2 and dynamic_cast<TH3*>(slice_NnD) != NULL) {
            TH1* h = dynamic_cast<TH3*>(slice_NnD)->Project3D("xy_temp_slice_projection");
            if (h) {
                h->SetName(Form("h1_slice_%s_%d_%d", GetSafeName().data(), indices[0], indices[1]));
                h->SetZTitle(Form("p(%s, %s | data)", GetParameter(indices[0]).GetLatexName().data(), GetParameter(indices[1]).GetLatexName().data()));
            }
            return h;
        }
    }

    // create local copy of parameters
    std::vector<double> parameters_temp = parameters;
    if (parameters_temp.empty()) {
        // check that remaining parameters are fixed
        for (unsigned i = 0; i < GetNParameters(); ++i)
            if (std::find(indices.begin(), indices.end(), i) == indices.end() and !GetParameter(i).Fixed()) {
                BCLog::OutError("BCIntegrate::GetSlice : All non-sliced parameters must be fixed or provided values in function call.");
                return NULL;
            }
        // set to fixed values
        parameters_temp = GetParameters().GetFixedValues();
    }

    // store previous binning preferences
    std::vector<unsigned> nbins_temp;
    for (unsigned i = 0; i < indices.size(); ++i) {
        nbins_temp.push_back(GetParameter(indices[i]).GetNbins());
        if (nbins > 0)								// set new binning
            GetParameter(indices[i]).SetNbins(nbins);
    }

    // create histogram
    TH1* h = NULL;
    if (indices.size() == 1) {
        h = GetParameter(indices[0]).CreateH1(Form("h1_slice_%s_%d", GetSafeName().data(), indices[0]));
        // set y-axis label
        if (GetNParameters() == 1)
            h->SetYTitle(Form("p(%s|data)", GetParameter(indices[0]).GetLatexName().data()));
        else
            h->SetYTitle(Form("p(%s|data, all other parameters fixed)", GetParameter(indices[0]).GetLatexName().data()));
    } else if (indices.size() == 2) {
        h = GetParameter(indices[0]).CreateH2(Form("h2_slice_%s_%d_%d", GetSafeName().data(), indices[0], indices[1]), GetParameter(indices[1]));
        // set z-axis label
        if (GetNParameters() == 2)
            h->SetZTitle(Form("p(%s, %s | data)", GetParameter(indices[0]).GetLatexName().data(), GetParameter(indices[1]).GetLatexName().data()));
        else
            h->SetZTitle(Form("p(%s, %s | data, all other parameters fixed)", GetParameter(indices[0]).GetLatexName().data(), GetParameter(indices[1]).GetLatexName().data()));
    } else if (indices.size() == 3) {
        h = GetParameter(indices[0]).CreateH3(Form("h3_slice_%s_%d_%d_%d", GetSafeName().data(), indices[0], indices[1], indices[2]), GetParameter(indices[1]), GetParameter(indices[2]));
        // set z-axis label
        if (GetNParameters() == 3)
            h->SetZTitle(Form("p(%s, %s, %s | data)", GetParameter(indices[0]).GetLatexName().data(), GetParameter(indices[1]).GetLatexName().data(), GetParameter(indices[2]).GetLatexName().data()));
        else
            h->SetZTitle(Form("p(%s, %s, %s | data, all other parameters fixed)", GetParameter(indices[0]).GetLatexName().data(), GetParameter(indices[1]).GetLatexName().data(), GetParameter(indices[2]).GetLatexName().data()));
    }

    // reset log_max_val
    log_max_val = -std::numeric_limits<double>::infinity();

    double log_min_val = std::numeric_limits<double>::infinity();

    // fill histogram
    // calculate number of bins
    int N = h->GetBin(h->GetNbinsX(), h->GetNbinsY(), h->GetNbinsZ());
    int bx, by, bz;
    // loop over all bins
    for  (int b = 1; b <= N; ++b) {
        // skip if bin is underflow or overflow
        if (h->IsBinUnderflow(b) or h->IsBinOverflow(b))
            continue;
        // get local bin numbers for each axis
        h->GetBinXYZ(b, bx, by, bz);
        // update x axis value
        parameters_temp[indices[0]] = h->GetXaxis()->GetBinCenter(bx);
        // update y axis value if 2D
        if (by > 0 and indices.size() > 1)
            parameters_temp[indices[1]] = h->GetYaxis()->GetBinCenter(by);
        // update z axis value if 3D
        if (bz > 0 and indices.size() > 2)
            parameters_temp[indices[2]] = h->GetZaxis()->GetBinCenter(bz);
        // calculate log of function value at parameters
        double log_eval = LogEval(parameters_temp);
        ++nIterations;

        // check max val
        log_max_val = std::max<double>(log_max_val, log_eval);
        // check min val
        log_min_val = std::min<double>(log_min_val, log_eval);
        // set bin content by global bin number
        h->SetBinContent(b, log_eval);
    }

    // remove log pedestal and exponentiate resulting value
    for (int b = 1; b <= N; ++b) {
        if (h->IsBinUnderflow(b) or h->IsBinOverflow(b))
            continue;
        h->SetBinContent(b, exp(h->GetBinContent(b) - log_max_val));
    }

    // normalize
    if (normalize) {
        double integral = h->Integral("width");
        if (integral != 0)
            h->Scale(1. / integral);
    }

    // reset binning
    for (unsigned i = 0; i < indices.size(); ++i)
        GetParameter(indices[i]).SetNbins(nbins_temp[i]);

    return h;
}

// ---------------------------------------------------------
TH2* BCIntegrate::GetSlice(unsigned index1, unsigned index2, unsigned& nIterations, double& log_max_val, const std::vector<double> parameters, int nbins, bool normalize)
{
    std::vector<unsigned> indices(1, index1);
    indices.push_back(index2);
    return (TH2*) GetSlice(indices, nIterations, log_max_val, parameters, nbins, normalize);
}

// ---------------------------------------------------------
double BCIntegrate::GetRandomPoint(std::vector<double>& x)
{
    GetRandomVectorInParameterSpace(x);
    return Eval(x);
}

// ---------------------------------------------------------
void BCIntegrate::GetRandomVectorInParameterSpace(std::vector<double>& x) const
{
    // get random vector in unit hypercube
    const_cast<TRandom3*>(&fRandom)->RndmArray(x.size(), &x.front());
    // translate it into values in ranges of the parameters
    fParameters.ValueFromPositionInRange(x);
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindMode(std::vector<double> start)
{
    if (GetNParameters() < 1) {
        BCLog::OutError("FindMode : No parameters defined. Aborting.");
        return std::vector<double>();
    }

    if (start.empty() and GetBestFitParameters().size() >= GetNParameters())
        start = GetBestFitParameters();

    if (start.size() > GetNParameters())
        start.resize(GetNParameters());

    std::vector<double> mode_temp(GetNParameters());
    std::vector<double> errors_temp(GetNParameters());
    BCIntegrate::BCOptimizationMethod method_temp = fOptimizationMethodCurrent;

    // output
    if (!(fOptimizationMethodCurrent == BCIntegrate::kOptDefault) && !(fOptimizationMethodCurrent == BCIntegrate::kOptEmpty))
        BCLog::OutSummary(Form("Finding mode using %s", DumpCurrentOptimizationMethod().c_str()));

    switch (fOptimizationMethodCurrent) {

        case BCIntegrate::kOptEmpty: {
            BCLog::OutWarning("BCIntegrate::FindMode : No optimization method chosen.");
            return std::vector<double>();
        }

        case BCIntegrate::kOptSimAnn: {
            FindModeSA(mode_temp, errors_temp, start);
            break;
        }

        case BCIntegrate::kOptMetropolis: {
            FindModeMCMC(mode_temp, errors_temp);
            break;
        }

        case BCIntegrate::kOptDefault:
            SetOptimizationMethod(BCIntegrate::kOptMinuit);
        case BCIntegrate::kOptMinuit: {
            int printlevel = -1;
            if (BCLog::GetLogLevelScreen() <= BCLog::detail)
                printlevel = 0;
            if (BCLog::GetLogLevelScreen() <= BCLog::debug)
                printlevel = 1;

            BCIntegrate::FindModeMinuit(mode_temp, errors_temp, start, printlevel);
            break;
        }

        default:
            BCLog::OutError(Form("BCIntegrate::FindMode : Invalid mode finding method: %d", GetOptimizationMethod()));
            return std::vector<double>();
    }

    // calculate function at new mode
    double fcnatmode_temp = Eval(mode_temp);

    // replace previous mode
    if (fFlagIgnorePrevOptimization or fcnatmode_temp > fLogMaximum) {
        SetBestFitParameters(mode_temp);
        fBestFitParameterErrors = errors_temp;
        fBestFitParameterErrors.resize(GetBestFitParameters().size(), std::numeric_limits<double>::infinity());
        fOptimizationMethodUsed = method_temp;
        fLogMaximum = fcnatmode_temp;
    }

    // return the new mode
    return fBestFitParameters;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindMode(BCIntegrate::BCOptimizationMethod optmethod, std::vector<double> start)
{
    // remember original method
    BCOptimizationMethod method_temp = fOptimizationMethodCurrent;

    // set method
    SetOptimizationMethod(optmethod);

    // run algorithm
    std::vector<double> mode = FindMode(start);

    // re-set original method
    SetOptimizationMethod(method_temp);

    // return mode
    return mode;
}

// ---------------------------------------------------------
TMinuit* BCIntegrate::GetMinuit()
{
    if (!fMinuit)
        fMinuit = new TMinuit();

    return fMinuit;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindModeMinuit(std::vector<double>& mode, std::vector<double>& errors, std::vector<double> start, int printlevel)
{
    if (fParameters.Size() < 1) {
        BCLog::OutError("BCIntegrate::FindModeMinuit : No parameters defined. Aborting.");
        return std::vector<double>();
    }

    // check start values
    if (!start.empty() and start.size() != fParameters.Size()) {
        BCLog::OutWarning("BCIntegrate::FindModeMinuit : Start point not valid (mismatch of dimensions), set to center.");
        start.clear();
    }


    // check if point is allowed
    if (!start.empty() and !GetParameters().IsWithinLimits(start)) {
        BCLog::OutWarning("BCIntegrate::FindModeMinuit : Start point not valid (parameter not inside valid range), set to center.");
        start.clear();
    }

    // check fixed values and issue warning before forcing to fixed
    if (!start.empty() and !GetParameters().IsAtFixedValues(start)) {
        BCLog::OutWarning("BCIntegrate::FindModeMinuit : Start point fixed values not properly set. Forcing to fixed values.");
        GetParameters().ApplyFixedValues(start);
    }

    if (start.empty())
        // if empty, set to center, with fixed values fixed
        start = GetParameters().GetRangeCenters();

    // set global this
    ::BCIntegrateHolder::instance(this);

    // define minuit
    delete fMinuit;
    fMinuit = new TMinuit(fParameters.Size());

    // set print level
    fMinuit->SetPrintLevel(printlevel);

    // set function
    fMinuit->SetFCN(&BCIntegrate::FCNLikelihood);

    // set UP for likelihood
    fMinuit->SetErrorDef(0.5);

    // set parameters
    int flag;
    for (unsigned i = 0; i < fParameters.Size(); i++)
        fMinuit->mnparm(i, GetParameter(i).GetName().data(), start[i],
                        GetParameter(i).GetRangeWidth() / 100.,
                        GetParameter(i).GetLowerLimit(),
                        GetParameter(i).GetUpperLimit(),
                        flag);

    for (unsigned i = 0; i < fParameters.Size(); i++)
        if (GetParameter(i).Fixed())
            fMinuit->FixParameter(i);

    // do mcmc minimization
    //   fMinuit->mnseek();

    // do minimization
    fMinuit->mnexcm("MIGRAD", fMinuitArglist, 2, flag);

    // improve search for local minimum
    //   fMinuit->mnimpr();

    // copy flag and result
    fMinuitErrorFlag = flag;
    //   std::vector<double> localMode(fParameters.Size(), 0);
    //   std::vector<double> errors(fParameters.Size(), 0);
    for (unsigned i = 0; i < fParameters.Size(); i++) {
        fMinuit->GetParameter(i, mode[i], errors[i]);
    }

    // delete minuit
    delete fMinuit;
    fMinuit = 0;

    return mode;
}

// --------------------------------------------------------
void BCIntegrate::WriteSAToFile(bool flag)
{
    if (flag)
        BCLog::OutError("BCIntegrate::WriteSAToFile: To turn on output use WriteSAToFile(filename,option).");
    fFlagWriteSAToFile = false;
}

// --------------------------------------------------------
void BCIntegrate::WriteSAToFile(const std::string& filename, const std::string& option, bool autoclose)
{
    if (filename.empty()) {
        BCLog::OutError("BCIntegrate::WriteSAToFile: You must specify the filename when turning on simlated annealing output.");
        return WriteSAToFile(false);
    }
    fSAOutputFilename = filename;
    fSAOutputFileOption = option;
    fSAOutputFileAutoclose = autoclose;
    fFlagWriteSAToFile = true;
}

// ---------------------------------------------------------
void BCIntegrate::InitializeSATree(bool replacetree, bool replacefile)
{
    if (fSATree and replacetree) {
        delete fSATree;
        fSATree = 0;
    }
    if (fSAOutputFile and replacefile) {
        fSAOutputFile->Close();
        delete fSAOutputFile;
        fSAOutputFile = 0;
    }

    TDirectory* dir = gDirectory;

    // create file
    if (!fSAOutputFile)
        fSAOutputFile = TFile::Open(fSAOutputFilename.c_str(), fSAOutputFileOption.c_str());
    // if failed
    if (!fSAOutputFile) {
        BCLog::OutError("BCIntegrate::InitializeSATree: Could not create new file.");
        WriteSAToFile(false);
        return;
    }
    // if file mode not writeable
    if (fSAOutputFile and !fSAOutputFile->IsWritable()) {
        BCLog::OutError("BCIntegrate::InitializeSATree: File must be opened in a writeable mode.");
        delete fSAOutputFile;
        fSAOutputFile = 0;
        WriteSAToFile(false);
        return;
    }

    if (!fSATree) {
        fSATree = new TTree(Form("%s_sa", GetSafeName().data()), Form("%s_sa", GetSafeName().data()));

        fSATree->Branch("Iteration",      &fSANIterations,   "Iteration/I");
        fSATree->Branch("Temperature",    &fSATemperature,   "Temperature/D");
        fSATree->Branch("LogProbability", &fSALogProb,       "LogProbability/D");

        fSAx.assign(fParameters.Size(), 0.0);
        for (unsigned i = 0; i < fParameters.Size(); ++i) {
            fSATree->Branch(GetParameter(i).GetSafeName().data(), &fSAx[i], (GetParameter(i).GetSafeName() + "/D").data());
            fSATree->SetAlias(Form("Parameter%i", i), GetParameter(i).GetSafeName().data());
        }
    }
    gDirectory = dir;
}

// ---------------------------------------------------------
void BCIntegrate::CloseSAOutputFile()
{
    if (!fSAOutputFile or !fSAOutputFile->IsOpen())
        return;
    fSAOutputFile->Write();
    fSAOutputFile->Close();
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindModeSA(std::vector<double>& mode, std::vector<double>& errors, std::vector<double> start)
{
    // note: if f(x) is the function to be minimized, then
    // f(x) := - LogEval(parameters)

    if (fFlagWriteSAToFile)
        InitializeSATree();

    std::vector<double> x, y; // vectors for current state, new proposed state and best fit up to now
    double fval_x, fval_y, fval_mode; // function values at points x, y and mode (we save them rather than to re-calculate them every time)
    int t = 1; // time iterator

    // check start values
    if (!start.empty() and start.size() != fParameters.Size()) {
        BCLog::OutWarning("BCIntegrate::FindModeSA : Start point not valid (mismatch of dimensions), set to center.");
        start.clear();
    }

    // check if point is allowed
    if (!start.empty() and !GetParameters().IsWithinLimits(start)) {
        BCLog::OutWarning("BCIntegrate::FindModeSA : Start point not valid (parameter not inside valid range), set to center.");
        start.clear();
    }

    // check fixed values and issue warning before forcing to fixed
    if (!start.empty() and !GetParameters().IsAtFixedValues(start)) {
        BCLog::OutWarning("BCIntegrate::FindModeSA : Start point fixed values not properly set. Forcing to fixed values.");
        GetParameters().ApplyFixedValues(start);
    }

    if (start.empty())
        // if empty, set to center, with fixed values fixed
        start = GetParameters().GetRangeCenters();

    // set current state and best fit to starting point
    x = start;
    mode = start;

    // calculate function value at starting point
    fval_x = fval_mode = LogEval(x);

    // run while still "hot enough"
    while ( SATemperature(t) > fSATmin ) {
        // generate new state
        y = GetProposalPointSA(x, t);

        // check if the proposed point is inside the phase space (ignoring fixed parameters)
        if (GetParameters().IsWithinLimits(y)) {
            // calculate function value at new point
            fval_y = LogEval(y);

            // is it better than the last one?
            // if so, update state and chef if it is the new best fit...
            if (fval_y >= fval_x) {
                x = y;

                fval_x = fval_y;

                if (fval_y > fval_mode) {
                    mode = y;
                    fval_mode = fval_y;
                }
            }
            // ...else, only accept new state w/ certain probability
            else {
                if (fRandom.Rndm() <= exp( (fval_y - fval_x) / SATemperature(t) )) {
                    x = y;
                    fval_x = fval_y;
                }
            }
        }

        // update tree variables
        fSANIterations = t;
        fSATemperature = SATemperature(t);
        fSALogProb = fval_x;
        fSAx = x;

        // fill tree
        if (fFlagWriteSAToFile && fSATree)
            fSATree->Fill();

        // increate t
        t++;
    }

    // calculate uncertainty
    errors.assign(fParameters.Size(), -1);

    if (fFlagWriteSAToFile and fSAOutputFileAutoclose)
        CloseSAOutputFile();

    return mode;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperature(double t) const
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
double BCIntegrate::SATemperatureBoltzmann(double t) const
{
    return fSAT0 / log(t + 1);
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCauchy(double t) const
{
    return fSAT0 / t;
}

// ---------------------------------------------------------
double BCIntegrate::SATemperatureCustom(double /*t*/) const
{
    BCLog::OutError("BCIntegrate::SATemperatureCustom : No custom temperature schedule defined");
    return 0.;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSA(const std::vector<double>& x, int t) const
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
std::vector<double> BCIntegrate::GetProposalPointSABoltzmann(const std::vector<double>& x, int t) const
{
    std::vector<double> y;
    y.clear();
    double new_val, norm;

    for (unsigned i = 0; i < fParameters.Size(); i++) {
        if (GetParameter(i).Fixed()) {
            y.push_back(GetParameter(i).GetFixedValue());
        } else {
            norm = GetParameter(i).GetRangeWidth() * SATemperature(t) / 2.;
            new_val = x[i] + norm * const_cast<TRandom3*>(&fRandom)->Gaus();
            y.push_back(new_val);
        }
    }

    return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACauchy(const std::vector<double>& x, int t) const
{
    std::vector<double> y;
    y.clear();

    if (GetNParameters() == 1) {
        double cauchy, new_val, norm;

        if (GetParameter(0).Fixed()) {
            y.push_back(GetParameter(0).GetFixedValue());
        } else {
            norm = GetParameter(0).GetRangeWidth() * SATemperature(t) / 2.;
            cauchy = tan(3.14159 * (const_cast<TRandom3*>(&fRandom)->Rndm() - 0.5));
            new_val = x[0] + norm * cauchy;
            y.push_back(new_val);
        }
    } else {
        // use sampling to get radial n-dim Cauchy distribution

        // first generate a random point uniformly distributed on a
        // fNvar-dimensional hypersphere
        y = SAHelperGetRandomPointOnHypersphere();

        // scale the vector by a random factor determined by the radial
        // part of the fNvar-dimensional Cauchy distribution
        double radial = SATemperature(t) * SAHelperGetRadialCauchy();

        // scale y by radial part and the size of dimension i in phase space
        // afterwards, move by x
        for (unsigned i = 0; i < GetNParameters(); i++) {
            if (GetParameter(i).Fixed()) {
                y[i] = GetParameter(i).GetFixedValue();
            } else {
                y[i] = GetParameter(i).GetRangeWidth() * y[i] * radial / 2. + x[i];
            }
        }
    }

    return y;
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::GetProposalPointSACustom(const std::vector<double>& /*x*/, int /*t*/) const
{
    BCLog::OutError("BCIntegrate::GetProposalPointSACustom : No custom proposal function defined");
    return std::vector<double>(fParameters.Size());
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::SAHelperGetRandomPointOnHypersphere() const
{
    std::vector<double> rand_point(fParameters.Size());

    // This method can only be called with fNvar >= 2 since the 1-dim case
    // is already hard wired into the Cauchy annealing proposal function.
    // To speed things up, hard-code fast method for 2 and dimensions.
    // The algorithm for 2D can be found at
    // http://mathworld.wolfram.com/CirclePointPicking.html
    // For 3D just using ROOT's algorithm.
    if (fParameters.Size() == 2) {
        double x1, x2, s;
        do {
            x1 = const_cast<TRandom3*>(&fRandom)->Rndm() * 2. - 1.;
            x2 = const_cast<TRandom3*>(&fRandom)->Rndm() * 2. - 1.;
            s = x1 * x1 + x2 * x2;
        } while (s >= 1);

        rand_point[0] = (x1 * x1 - x2 * x2) / s;
        rand_point[1] = (2.*x1 * x2) / s;
    } else if (fParameters.Size() == 3) {
        const_cast<TRandom3*>(&fRandom)->Sphere(rand_point[0], rand_point[1], rand_point[2], 1.0);
    } else {
        double s = 0.,
               gauss_num;

        for (unsigned i = 0; i < fParameters.Size(); i++) {
            gauss_num = const_cast<TRandom3*>(&fRandom)->Gaus();
            rand_point[i] = gauss_num;
            s += gauss_num * gauss_num;
        }
        s = sqrt(s);

        for (unsigned i = 0; i < fParameters.Size(); i++)
            rand_point[i] = rand_point[i] / s;
    }

    return rand_point;
}

// ---------------------------------------------------------
double BCIntegrate::SAHelperGetRadialCauchy() const
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
    static unsigned map_dimension = 0;

    // is the lookup-table already initialized? if not, do it!
    if (!initialized or map_dimension != fParameters.Size()) {
        double init_theta;
        double beta = SAHelperSinusToNIntegral(fParameters.Size() - 1, 1.57079632679);

        for (int i = 0; i <= 10000; i++) {
            init_theta = 3.14159265 * (double)i / 5000.;
            map_theta[i] = init_theta;

            map_u[i] = SAHelperSinusToNIntegral(fParameters.Size() - 1, init_theta) / beta;
        }

        map_dimension = fParameters.Size();
        initialized = true;
    } // initializing is done.

    // generate uniform random number for sampling
    double u = const_cast<TRandom3*>(&fRandom)->Rndm();

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
double BCIntegrate::SAHelperSinusToNIntegral(int dim, double theta) const
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
void BCIntegrate::FCNLikelihood(int& /*npar*/, double* /*grad*/, double& fval, double* par, int /*flag*/)
{
    // copy parameters
    static std::vector<double> parameters;

    // calculate number of active + fixed parameters
    // remember: npar is just the number of _active_ parameters while
    // par is a vector of _all_ parameters
    int nparameters = ::BCIntegrateHolder::instance()->GetNParameters();

    // adjust size if needed
    parameters.resize(nparameters, 0.0);

    // copy values
    std::copy(par, par + nparameters, parameters.begin());

    // evaluate, for efficiency don't check if npar matches
    fval = - ::BCIntegrateHolder::instance()->LogEval(parameters);
}

// ---------------------------------------------------------
std::vector<double> BCIntegrate::FindModeMCMC(std::vector<double>& mode, std::vector<double>& errors)
{
    // call PreRun
    MetropolisPreRun();

    mode = fMCMCStatistics_AllChains.mode;
    errors.assign(fParameters.Size(), -1.);

    return mode;
}

// ---------------------------------------------------------
void BCIntegrate::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
    if (method >= BCIntegrate::NIntMethods) {
        BCLog::OutError(Form("BCIntegrate::SetIntegrationMethod: Invalid method '%d' ", method));
        return;
    }
    if (method == BCIntegrate::kIntCuba) {
#ifndef HAVE_CUBA_H
        BCLog::OutError("BCIntegrate::SetIntegrationMethod: Cuba not enabled during configure");
#endif
    }
    fIntegrationMethodCurrent = method;
}

// ---------------------------------------------------------
void BCIntegrate::SetCubaIntegrationMethod(BCIntegrate::BCCubaMethod type)
{
#ifdef HAVE_CUBA_H
    switch (type) {
        case BCIntegrate::kCubaVegas:
        case BCIntegrate::kCubaSuave:
        case BCIntegrate::kCubaDivonne:
        case BCIntegrate::kCubaCuhre:
        case BCIntegrate::kCubaDefault:
            fCubaIntegrationMethod = type;
            return;
        default:
            BCLog::OutError(Form("Integration method of type %d is not defined for Cuba", type));
            return;
    }
#else
    (void) type; // suppress compiler warning about unused parameters
    BCLog::OutError("SetCubaIntegrationMethod: Cuba not enabled during configure");
#endif
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateCuba(BCCubaMethod cubatype)
{
#if HAVE_CUBA_H
    // integrand has only one component
    static const int ncomp = 1;

    // don't store integration in a slot for reuse
    static const int gridno = -1;

    // we don't want dumps of internal state
    static const char* statefile = "";

    // only evaluate at one position in parameter space at a time
    static const int nvec = 1;

#if CUBAVERSION > 33
    // we have no spinning cores
    int* spin = 0;
#endif

    // cuba returns info in these variables
    int fail = 0;
    int nregions = 0;
    std::vector<double> integral(ncomp, -1);
    std::vector<double> error(ncomp, -1);
    std::vector<double> prob(ncomp, -1);

    // reset number of iterations
    fNIterations = 0;

    // Cuba needs int variable
    int nIntegrationVariables = GetParameters().GetNFreeParameters();

    if (cubatype == BCIntegrate::kCubaDefault) {
        switch (nIntegrationVariables) {
            case 1:
                cubatype = BCIntegrate::kCubaVegas;
                break;
            case 2:
            case 3:
                cubatype = BCIntegrate::kCubaCuhre;
                break;
            default:
                cubatype = BCIntegrate::kCubaDivonne;
                break;
        }
        if (nIntegrationVariables > 33)
            cubatype = BCIntegrate::kCubaVegas;
    }
    LogOutputAtStartOfIntegration(kIntCuba, cubatype);

    switch (cubatype) {

        case BCIntegrate::kCubaVegas:
            Vegas(nIntegrationVariables, ncomp,
                  &BCIntegrate::CubaIntegrand, static_cast<void*>(this),
                  nvec,
                  fRelativePrecision, fAbsolutePrecision,
                  fCubaVegasOptions.flags, fRandom.GetSeed(),
                  fNIterationsMin, fNIterationsMax,
                  fCubaVegasOptions.nstart, fCubaVegasOptions.nincrease, fCubaVegasOptions.nbatch,
                  gridno, statefile,
#if CUBAVERSION > 33
                  spin,
#endif
                  &fNIterations, &fail,
                  &integral[0], &error[0], &prob[0]);
            break;

        case BCIntegrate::kCubaSuave:
            Suave(nIntegrationVariables, ncomp,
                  &BCIntegrate::CubaIntegrand, static_cast<void*>(this),
                  nvec,
                  fRelativePrecision, fAbsolutePrecision,
                  fCubaSuaveOptions.flags, fRandom.GetSeed(),
                  fNIterationsMin, fNIterationsMax, fCubaSuaveOptions.nnew,
#if CUBAVERSION > 40
                  fCubaSuaveOptions.nmin,
#endif
                  fCubaSuaveOptions.flatness, statefile,
#if CUBAVERSION > 33
                  spin,
#endif
                  &nregions, &fNIterations, &fail,
                  &integral[0], &error[0], &prob[0]);
            break;

        case BCIntegrate::kCubaDivonne:
            if (nIntegrationVariables < 2 || nIntegrationVariables > 33)
                BCLOG_WARNING("Divonne only works in 2 <= d <= 33 dimensions");
            else {
                // no extra info supported
                static const int ngiven = 0;
                static const int nextra = ngiven;
                Divonne(nIntegrationVariables, ncomp,
                        &BCIntegrate::CubaIntegrand, static_cast<void*>(this),
                        nvec,
                        fRelativePrecision, fAbsolutePrecision,
                        fCubaDivonneOptions.flags, fRandom.GetSeed(),
                        fNIterationsMin, fNIterationsMax,
                        fCubaDivonneOptions.key1, fCubaDivonneOptions.key2, fCubaDivonneOptions.key3,
                        fCubaDivonneOptions.maxpass, fCubaDivonneOptions.border,
                        fCubaDivonneOptions.maxchisq, fCubaDivonneOptions.mindeviation,
                        ngiven, nIntegrationVariables /*ldxgiven*/, NULL, nextra, NULL,
                        statefile,
#if CUBAVERSION > 33
                        spin,
#endif
                        &nregions, &fNIterations, &fail,
                        &integral[0], &error[0], &prob[0]);
            }
            break;

        case BCIntegrate::kCubaCuhre:
            if (nIntegrationVariables < 2)
                BCLog::OutError("BCIntegrate::IntegrateCuba(Cuhre): Cuhre(cubature) only works in d > 1");

            Cuhre(nIntegrationVariables, ncomp,
                  &BCIntegrate::CubaIntegrand, static_cast<void*>(this),
                  nvec,
                  fRelativePrecision, fAbsolutePrecision,
                  fCubaCuhreOptions.flags, fNIterationsMin, fNIterationsMax,
                  fCubaCuhreOptions.key,
                  statefile,
#if CUBAVERSION > 33
                  spin,
#endif
                  &nregions, &fNIterations, &fail,
                  &integral[0], &error[0], &prob[0]);
            break;

        case BCIntegrate::NCubaMethods:
        default:
            BCLog::OutError("Cuba integration method not set.");
            error[0] = -1;
            integral[0] = -1;
            break;
    }

    fError = error[0];
    double result = integral[0];

    if (fail != 0) {
        BCLog::OutWarning(" Warning, integral did not converge with the given set of parameters. ");
        BCLog::OutWarning(Form(" neval    = %d", fNIterations));
        BCLog::OutWarning(Form(" fail     = %d", fail));
        BCLog::OutWarning(Form(" integral = %e", result));
        BCLog::OutWarning(Form(" error    = %e", fError));
        BCLog::OutWarning(Form(" prob     = %e", prob[0]));

        // handle cases in which cuba does not alter integral
        if (fNIterations == 0) {
            fError = -1;
            result = -1;
        }

    } else
        LogOutputAtEndOfIntegration(result, fError, fError / result, fNIterations);

    return result;
#else
    (void) cubatype; // suppress compiler warning about unused parameters
    BCLog::OutError("IntegrateCuba: Cuba not enabled during configure");
    return -1;
#endif
}

// ---------------------------------------------------------
int BCIntegrate::CubaIntegrand(const int* ndim, const double xx[],
                               const int* /*ncomp*/, double ff[], void* userdata)
{
    BCIntegrate* local_this = static_cast<BCIntegrate*>(userdata);

    // scale variables
    double jacobian = 1.0;

    // create local parameter vector
    // important for thread safety, though not super efficient
    std::vector<double> scaled_parameters(local_this->fParameters.Size());

    // stay in sync with the possible lower number of parameters
    // that cuba sees due to fixing in BAT
    unsigned cubaIndex = 0;
    unsigned batIndex = 0;
    for (batIndex = 0; batIndex < local_this->fParameters.Size(); ++batIndex) {
        const BCParameter& p = local_this->GetParameter(batIndex);

        // get the scaled parameter value
        if (p.Fixed())
            scaled_parameters[batIndex] = p.GetFixedValue();
        else {
            // convert from unit hypercube to actual parameter hyperrectangle
            scaled_parameters[batIndex] = p.GetLowerLimit() + xx[cubaIndex] * p.GetRangeWidth();

            // multiply range to jacobian
            jacobian *= p.GetRangeWidth();

            // one more parameter that cuba varies
            ++cubaIndex;
        }
    }

    if (cubaIndex != unsigned(*ndim))
        BCLog::OutError(Form("BCIntegrate::CubaIntegrand: mismatch between variable parameters"
                             "in BAT (%d) and Cuba(%d)", batIndex, cubaIndex));

    // call function to integrate
    ff[0] = local_this->Eval(scaled_parameters);

    // multiply jacobian
    ff[0] *= jacobian;

    return 0;
}

// ---------------------------------------------------------
double BCIntegrate::IntegrateSlice()
{
    // print to log
    LogOutputAtStartOfIntegration(fIntegrationMethodCurrent, NCubaMethods);

    double integral = -1;
    fError = -1;

    // precision isn't estimated here
    double absprecision  = -1;
    double relprecision  = -1;

    // get vector of fixed values
    std::vector<double> fixpoint = GetParameters().GetFixedValues();

    // get vector of indices of unfixed parameters
    std::vector<unsigned> indices;
    for (unsigned i = 0; i < GetNParameters(); ++i)
        if (!GetParameter(i).Fixed())
            indices.push_back(i);

    // check size of vector of indices
    if (indices.size() > 3) {
        BCLog::OutWarning("BCIntegrate::IntegrateSlice: Method only implemented for 1D, 2D, and 3D problems. Return -1.");
        integral = -1;
    }

    unsigned nIterations = 0;
    double log_max_val = -std::numeric_limits<double>::infinity();

    // get slice, without normalizing
    TH1* h = GetSlice(indices, nIterations, log_max_val, fixpoint, 0, false);
    integral = h->Integral("width") * exp(log_max_val);

    // print to log
    LogOutputAtEndOfIntegration(integral, absprecision, relprecision, nIterations);

    return integral;
}


// ---------------------------------------------------------
std::string BCIntegrate::DumpIntegrationMethod(BCIntegrate::BCIntegrationMethod type) const
{
    switch (type) {
        case BCIntegrate::kIntEmpty:
            return "Empty";
        case BCIntegrate::kIntMonteCarlo:
            return "Sample Mean Monte Carlo";
        case BCIntegrate::kIntCuba:
            return "Cuba";
        case BCIntegrate::kIntGrid:
            return "Grid";
        default:
            return "Undefined";
    }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpMarginalizationMethod(BCIntegrate::BCMarginalizationMethod type) const
{
    switch (type) {
        case BCIntegrate::kMargEmpty:
            return "Empty";
        case BCIntegrate::kMargMonteCarlo:
            return "Sample Mean Monte Carlo";
        case BCIntegrate::kMargMetropolis:
            return "Metropolis";
        case BCIntegrate::kMargGrid:
            return "Grid";
        case BCIntegrate::kMargDefault:
            return "Default";
        default:
            return "Undefined";
    }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpOptimizationMethod(BCIntegrate::BCOptimizationMethod type) const
{
    switch (type) {
        case BCIntegrate::kOptEmpty:
            return "Empty";
        case BCIntegrate::kOptSimAnn:
            return "Simulated Annealing";
        case BCIntegrate::kOptMetropolis:
            return "Metropolis MCMC";
        case BCIntegrate::kOptMinuit:
            return "Minuit";
        case BCIntegrate::kOptDefault:
            return "Default";
        default:
            return "Undefined";
    }
}

// ---------------------------------------------------------
std::string BCIntegrate::DumpCubaIntegrationMethod(BCIntegrate::BCCubaMethod type) const
{
    switch (type) {
        case BCIntegrate::kCubaVegas:
            return "Vegas";
        case BCIntegrate::kCubaSuave:
            return "Suave";
        case BCIntegrate::kCubaDivonne:
            return "Divonne";
        case BCIntegrate::kCubaCuhre:
            return "Cuhre";
        case BCIntegrate::kCubaDefault:
            return "Default";
        default:
            return "Undefined";
    }
}

// ---------------------------------------------------------
BCCubaOptions::General::General():
    ncomp(1),
    flags(0),
    nregions(0),
    neval(0),
    fail(0),
    error(0),
    prob(0)
{}

// ---------------------------------------------------------
BCCubaOptions::General::~General()
{}

/*
 * Copy default values from sec. 7.3.1 of the cuba manual for version 4.2.
 */

// ---------------------------------------------------------
BCCubaOptions::Vegas::Vegas():
    General(),
    nstart(1000),
    nincrease(500),
    nbatch(1000),
    gridno(0)
{}

// ---------------------------------------------------------
BCCubaOptions::Suave::Suave():
    General(),
    nmin(2),
    nnew(1000),
    flatness(50)
{}

// ---------------------------------------------------------
BCCubaOptions::Divonne::Divonne():
    General(),
    key1(47),
    key2(1),
    key3(1),
    maxpass(5),
    border(0),
    maxchisq(10),
    mindeviation(0.25)
{}

// ---------------------------------------------------------
BCCubaOptions::Cuhre::Cuhre():
    General(),
    key(0) // let cuba choose default cubature rule
{}

// ---------------------------------------------------------
void BCIntegrate::PrintMarginalizationSummary() const
{
    if (GetIntegral() >= 0) {
        BCLog::OutSummary(" Results of the integration");
        BCLog::OutSummary(" ============================");
        BCLog::OutSummary(" Integration method used: " + DumpUsedIntegrationMethod());
        if (GetError() >= 0)
            BCLog::OutSummary(Form(" Evidence: %f +- %f", GetIntegral(), GetError()));
        else
            BCLog::OutSummary(Form(" Evidence: %f (no error estimate available)", GetIntegral()));
    }

    if (fFlagMarginalized) {
        BCLog::OutSummary(" Marginalization algorithm used: " + DumpUsedMarginalizationMethod());
        BCLog::OutSummary("");
        BCEngineMCMC::PrintMarginalizationSummary();
    }
}

// ---------------------------------------------------------
void BCIntegrate::PrintBestFitSummary() const
{
    if (GetBestFitParameters().empty()) {
        BCLog::OutSummary("No best fit information available.");
        return;
    }

    BCLog::OutSummary(" Results of the optimization");
    BCLog::OutSummary(" ===========================");
    BCLog::OutSummary(" Optimization algorithm used: " + DumpUsedOptimizationMethod());

    BCEngineMCMC::PrintBestFitSummary();
}

// ---------------------------------------------------------
std::string BCIntegrate::GetBestFitSummary(unsigned i) const
{
    if (i >= GetNVariables())
        return std::string("");
    if (i < GetBestFitParameterErrors().size() and GetBestFitParameterErrors()[i] != std::numeric_limits<double>::infinity())
        return BCEngineMCMC::GetBestFitSummary(i) + Form(" +- %.*g", GetVariable(i).GetPrecision(), GetBestFitParameterErrors()[i]);
    else if (!GetBestFitParameterErrors().empty())
        return BCEngineMCMC::GetBestFitSummary(i) + " (no error estimate available)";
    else
        return BCEngineMCMC::GetBestFitSummary(i);
}
