/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCPriorModel.h"

#include <TH1.h>

// ---------------------------------------------------------
BCPriorModel::BCPriorModel(BCModel& model, bool call_likelihood) :
    BCModel(model.GetName() + "_prior"),
    fModel(model),
    fCallLikelihood(call_likelihood)
{
    PreparePriorModel();
}

// ---------------------------------------------------------
bool BCPriorModel::PreparePriorModel()
{
    // copy parameters and observables
    fParameters = fModel.GetParameters();
    fObservables = BCObservableSet(fModel.GetObservables(), true);

    // but use binning that was used by model
    for (unsigned i = 0; i < fModel.GetNVariables(); ++i)
        if (fModel.MarginalizedHistogramExists(i))
            // Set binning from marginalization rather than parameter copy
            GetVariable(i).SetNbins(fModel.GetMarginalizedHistogram(i)->GetNbinsX());

    // set default MCMC setup to the one of the original model
    SetPrecision(fModel);

    return true;
}

// ---------------------------------------------------------
void BCPriorModel::CalculateObservables(const std::vector<double>& parameters)
{
    // help with thread safety and correctness of sampling
    fModel.UpdateChainIndex(GetCurrentChain());

    if (fCallLikelihood)
        fModel.LogLikelihood(parameters);
    fModel.CalculateObservables(parameters);
    for (unsigned i = 0; i < GetNObservables(); ++i)
        GetObservable(i).Value(fModel.GetObservable(i).Value());
}
