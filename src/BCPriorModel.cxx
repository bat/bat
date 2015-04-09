/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPriorModel.h"

#include <TH1D.h>

// ---------------------------------------------------------
BCPriorModel::BCPriorModel(BCModel* model, bool call_likelihood)
    : BCModel("prior_model")
    , fModel(model)
    , fCallLikelihood(call_likelihood)
{
    if (fModel) {
        SetName((fModel->GetName() + "_prior").data());
        PreparePriorModel();
    }
}

// ---------------------------------------------------------
BCPriorModel::~BCPriorModel()
{
}

// ---------------------------------------------------------
bool BCPriorModel::PreparePriorModel()
{
    // Clear Parameters and User-defined Observables
    fParameters = BCParameterSet();
    fObservables = BCVariableSet();

    if (!fModel)
        return false;

    // COPY parameters from model
    for (unsigned i = 0; i < fModel->GetNParameters(); ++i) {
        AddParameter(new BCParameter(*(fModel->GetParameter(i))));
        if (fModel->MarginalizedHistogramExists(i))
            // Set binning from marginalization rather than parameter copy
            GetParameter(i)->SetNbins(fModel->GetMarginalizedHistogram(i)->GetNbinsX());
    }
    // COPY user-observables from model
    for (unsigned i = 0; i < fModel->GetNObservables(); ++i) {
        AddObservable(new BCObservable(*(fModel->GetObservable(i))));
        if (fModel->MarginalizedHistogramExists(i + fModel->GetNParameters()))
            // Set binning from marginalization rather than observable copy
            GetObservable(i)->SetNbins(fModel->GetMarginalizedHistogram(i + fModel->GetNParameters())->GetNbinsX());
    }

    // set default MCMC setup to the one of the original model
    MCMCSetPrecision(fModel);

    return true;
}

// ---------------------------------------------------------
void BCPriorModel::CalculateObservables(const std::vector<double>& parameters)
{
    if (!fModel)
        return;
    fModel->MCMCSetCurrentChain(fMCMCCurrentChain);
    if (fCallLikelihood)
        fModel->LogLikelihood(parameters);
    fModel->CalculateObservables(parameters);
    for (unsigned i = 0; i < GetNObservables(); ++i)
        GetObservable(i)->Value(fModel->GetObservable(i)->Value());
}
