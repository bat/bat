/*
 * Copyright (C) 2007-2013, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCSummaryPriorModel.h"

#include "BCH1D.h"
#include "BCLog.h"
#include "BCParameter.h"
#include "BCUserObservable.h"

#include <TH1D.h>

// ---------------------------------------------------------
BCSummaryPriorModel::BCSummaryPriorModel()
   : BCModel()
   , fTestModel(0)
{
}

// ---------------------------------------------------------
BCSummaryPriorModel::BCSummaryPriorModel(const char * name)
   : BCModel(name)
   , fTestModel(0)
{
}

// ---------------------------------------------------------
BCSummaryPriorModel::BCSummaryPriorModel(BCModel * model, const char * name)
	: BCModel(name),
		fTestModel(0)
{
	SetModel(model);
}

// ---------------------------------------------------------
BCSummaryPriorModel::~BCSummaryPriorModel()
{}

// ---------------------------------------------------------
void BCSummaryPriorModel::SetModel(BCModel * model)
{
   fTestModel = model;

   // copy parameters
   for (unsigned i = 0; i < fTestModel->GetNParameters(); ++i) {
		 BCParameter * par = const_cast<BCParameter *>(fTestModel->GetParameter(i));
		 if (fTestModel->MarginalizedHistogramExists(i))
			 // Set binning from marginalization rather than parameter copy
			 par -> SetNbins(fTestModel->GetMarginalizedHistogram(i)->GetNbinsX());
		 AddParameter(par);
   }	
	 // copy user-observables
	 for (unsigned i = 0; i < fTestModel->GetNUserObservables(); ++i) {
	 	 BCUserObservable * obs = const_cast<BCUserObservable *>(fTestModel->GetUserObservable(i));
		 if (fTestModel->MarginalizedHistogramExists(i+fTestModel->GetNParameters()))
			 // Set binning from marginalization rather than observable copy
			 obs -> SetNbins(fTestModel->GetMarginalizedHistogram(i+fTestModel->GetNParameters())->GetNbinsX());
	 	 AddUserObservable(obs);
	 }
	 
   // set default MCMC setup to the one of the original model
	 MCMCSetPrecision(fTestModel);
}

// ---------------------------------------------------------
double BCSummaryPriorModel::LogLikelihood(const std::vector<double> & parameters)
{
   return fTestModel->LogAPrioriProbability(parameters);
}

// ---------------------------------------------------------
double BCSummaryPriorModel::LogAPrioriProbability(const std::vector<double> & /*parameters*/)
{
   return 0;
}
