/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCModelManager.h"

#include "BCDataPoint.h"
#include "BCLog.h"

#include <TString.h>

// ---------------------------------------------------------
BCModelManager::BCModelManager()
    : fDataSet(NULL)
{
}

// ---------------------------------------------------------
BCModelManager::BCModelManager(const BCModelManager& other)
    : fModels(other.fModels),
      fAPrioriProbability(other.fAPrioriProbability),
      fAPosterioriProbability(other.fAPosterioriProbability),
      fDataSet(other.fDataSet)
{
}

// ---------------------------------------------------------
BCModelManager::~BCModelManager()
{
}

// ---------------------------------------------------------
BCModelManager& BCModelManager::operator=(BCModelManager rhs)
{
    swap(*this, rhs);
    return *this;
}

// ---------------------------------------------------------
void swap(BCModelManager& A, BCModelManager& B)
{
    std::swap(A.fModels,                 B.fModels);
    std::swap(A.fAPrioriProbability,     B.fAPrioriProbability);
    std::swap(A.fAPosterioriProbability, B.fAPosterioriProbability);
    std::swap(A.fDataSet,                B.fDataSet);
}

// ---------------------------------------------------------
void BCModelManager::SetDataSet(BCDataSet* dataset)
{
    // set data set
    fDataSet = dataset;

    // set data set of all models in the manager
    for (unsigned int i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetDataSet(fDataSet);
}

// ---------------------------------------------------------
void BCModelManager::AddModel(BCModel* model, double probability)
{
    // set data set
    if (fDataSet)
        model->SetDataSet(fDataSet);

    // fill model into container
    fModels.push_back(model);
    // set probabilities
    fAPrioriProbability.push_back(probability);
    fAPosterioriProbability.push_back(-1);
}

// ---------------------------------------------------------
void BCModelManager::SetPrecision(BCEngineMCMC::Precision precision)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetPrecision(precision);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsMax(int niterations)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetNIterationsMax(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsMin(int niterations)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetNIterationsMin(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetNIterationsPrecisionCheck(int niterations)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetNIterationsPrecisionCheck(niterations);
}

// ---------------------------------------------------------
void BCModelManager::SetIntegrationMethod(BCIntegrate::BCIntegrationMethod method)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetIntegrationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetMarginalizationMethod(BCIntegrate::BCMarginalizationMethod method)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetMarginalizationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetOptimizationMethod(BCIntegrate::BCOptimizationMethod method)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetOptimizationMethod(method);
}

// ---------------------------------------------------------
void BCModelManager::SetRelativePrecision(double relprecision)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetRelativePrecision(relprecision);
}

// ---------------------------------------------------------
void BCModelManager::SetAbsolutePrecision(double absprecision)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetAbsolutePrecision(absprecision);
}

// ---------------------------------------------------------
void BCModelManager::SetNbins(unsigned int n)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetNbins(n);
}

// ---------------------------------------------------------
void BCModelManager::SetNChains(unsigned int n)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->SetNChains(n);
}

// ---------------------------------------------------------
void BCModelManager::Integrate()
{
    // initialize likelihood norm
    double normalization = 0.0;

    BCLog::OutSummary("Running normalization of all models.");

    for (unsigned int i = 0; i < GetNModels(); ++i) {
        fModels[i]->Integrate();

        // add to total normalization
        normalization += fModels[i]->GetIntegral() * fAPrioriProbability[i];
    }

    // set model a posteriori probabilities
    for (unsigned int i = 0; i < GetNModels(); ++i)
        fAPosterioriProbability[i] = (fModels[i]->GetIntegral() * fAPrioriProbability[i]) / normalization;
}

// ---------------------------------------------------------
double BCModelManager::BayesFactor(unsigned imodel1, unsigned imodel2) const
{
    if (imodel1 >= fModels.size() || imodel2 >= fModels.size())
        return -1;

    // Bayes Factors are the likelihoods integrated over the parameters
    // Is this equal to the posteriors?
    //    NOOOO
    // But it is equal to normalization factors.

    // check model 1
    const double norm1 = fModels[imodel1]->GetIntegral();
    if (norm1 < 0.) {
        BCLOG_ERROR(Form(" Model %s (index %d) not normalized. Cannot calculate Bayes factor.", fModels[imodel1]->GetName().data(), imodel1));
        return -1.;
    }

    // check model 2
    const double norm2 = fModels[imodel2]->GetIntegral();
    if (norm2 < 0) {
        BCLOG_ERROR(Form(" Model %s (index %d) not normalized. Cannot calculate Bayes factor.", fModels[imodel2]->GetName().data(), imodel2));
        return -1.;
    }

    // denominator cannot be zero
    if (norm2 < std::numeric_limits<double>::min() && norm1 >= std::numeric_limits<double>::min()) {
        BCLOG_ERROR(Form(" Model %s (index %d) has ZERO probability. Bayes factor is infinite.", fModels[imodel2]->GetName().data(), imodel2));
        return -1.;
    }

    // denominator cannot be zero unless also numerator is zero
    if (norm2 < std::numeric_limits<double>::min() && norm1 < std::numeric_limits<double>::min()) {
        BCLOG_WARNING(Form("Models %s and %s have ZERO probability. Bayes factor is unknown. Returning 1.", fModels[imodel2]->GetName().data(), fModels[imodel1]->GetName().data()));
        return 1.;
    }

    // now calculate the factor
    return norm1 / norm2;
}

// ---------------------------------------------------------
void BCModelManager::FindMode()
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->FindMode();
}

// ---------------------------------------------------------
void BCModelManager::MarginalizeAll()
{
    // marginalizes all models registered
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->MarginalizeAll();
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChain(bool flag)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->WriteMarkovChain(flag);
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChainRun(bool flag)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->WriteMarkovChainRun(flag);
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChainPreRun(bool flag)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->WriteMarkovChainPreRun(flag);
}

// ---------------------------------------------------------
void BCModelManager::WriteMarkovChain(const std::string& prefix, const std::string& option, bool flag_run, bool flag_prerun)
{
    for (unsigned i = 0; i < GetNModels(); ++i)
        GetModel(i)->WriteMarkovChain(prefix + GetModel(i)->GetSafeName() + ".root", option, flag_run, flag_prerun);
}

// ---------------------------------------------------------
void BCModelManager::PrintSummary() const
{
    BCLog::OutSummary("");
    BCLog::OutSummary("======================================");
    BCLog::OutSummary(" Summary");
    BCLog::OutSummary("======================================");
    BCLog::OutSummary("");
    BCLog::OutSummary(Form(" Number of models               %lu: ", fModels.size()));
    BCLog::OutSummary("");

    BCLog::OutSummary(" - Models:");
    BCLog::OutSummary("");
    for (unsigned i = 0; i < fModels.size(); ++i)
        fModels[i]->PrintSummary();

    BCLog::OutSummary(" - Data:");
    BCLog::OutSummary("");
    BCLog::OutSummary(Form("     Number of entries: %u", fDataSet->GetNDataPoints()));
    BCLog::OutSummary("");

    PrintModelComparisonSummary();
}

// ---------------------------------------------------------
void BCModelManager::PrintModelComparisonSummary() const
{
    BCLog::OutSummary("======================================");
    BCLog::OutSummary(" Model comparison:");
    BCLog::OutSummary("");

    BCLog::OutSummary(" - A priori probabilities:");
    for (unsigned i = 0; i < fModels.size(); ++i)
        BCLog::OutSummary(Form(" (%u) p(%s) = %f", i, fModels[i]->GetName().data(), fAPrioriProbability[i]));

    BCLog::OutSummary(" - A posteriori probabilities:");
    for (unsigned i = 0; i < fModels.size(); ++i)
        BCLog::OutSummary(Form(" (%u) p(%s | data) = %f", i, fModels[i]->GetName().data(), fAPosterioriProbability[i]));

    // Bayes factors summary
    BCLog::OutSummary(" - Bayes factors:");
    for (unsigned i = 0; i < fModels.size(); ++i)
        for (unsigned j = i + 1; j < fModels.size(); ++j)
            BCLog::OutSummary(Form("     K := p(data | %s) / p(data | %s) = %f", fModels[i]->GetName().data(), fModels[j]->GetName().data(),
                                   BayesFactor(i, j)));
    BCLog::OutSummary("===========================================");
}
