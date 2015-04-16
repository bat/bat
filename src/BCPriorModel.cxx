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
    if (!fModel)
        return false;

    // copy parameters and observables
    fParameters = fModel->GetParameters();
    fObservables = fModel->GetObservables();


    // but use binning that was used by model
    for (unsigned i = 0; i < fModel->GetNVariables(); ++i)
        if (fModel->MarginalizedHistogramExists(i))
            // Set binning from marginalization rather than parameter copy
            GetVariable(i).SetNbins(fModel->GetMarginalizedHistogram(i)->GetNbinsX());

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
        GetObservable(i).Value(fModel->GetObservable(i).Value());
}
