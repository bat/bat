/*
 * Copyright (C) 2007-2014, the BAT core developer team
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
BCSummaryPriorModel::~BCSummaryPriorModel()
{}

// ---------------------------------------------------------
void BCSummaryPriorModel::SetModel(BCModel * model)
{
   fTestModel = model;

   // copy parameters
   int npar = fTestModel->GetNParameters();
   for (int i = 0; i < npar; ++i) {
      BCParameter * par = const_cast<BCParameter *>(fTestModel->GetParameter(i));
      AddParameter(par);
   }

   // set default histogram binning to the one of the original model
   for (int i = 0; i < npar; ++i) {
      // this construct has to go here, because otherwise there is a
      // warning from BCEngineMCMC:: MCMCGetH1Marginalized
      if (fTestModel->GetBestFitParameters().size() > 0){
         BCH1D* hist = fTestModel->GetMarginalized( fTestModel->GetParameter(i) );
         if (hist) {
            int nbins = hist->GetHistogram()->GetNbinsX();
            fTestModel->GetParameter(i)->SetNbins(nbins);
         }
      }
   }

   // set default MCMC setup to the one of the original model
   MCMCSetNChains( fTestModel->MCMCGetNChains() );
   MCMCSetNLag( fTestModel->MCMCGetNLag() );
   MCMCSetNIterationsMax( fTestModel->MCMCGetNIterationsMax() );
   MCMCSetNIterationsRun( fTestModel->MCMCGetNIterationsRun() );
   MCMCSetNIterationsPreRunMin( fTestModel->MCMCGetNIterationsPreRunMin() );
   MCMCSetNIterationsUpdate( fTestModel->MCMCGetNIterationsUpdate() );
   MCMCSetNIterationsUpdateMax( fTestModel->MCMCGetNIterationsUpdateMax() );
   MCMCSetRValueCriterion( fTestModel->MCMCGetRValueCriterion() );
   MCMCSetRValueParametersCriterion( fTestModel->MCMCGetRValueParametersCriterion() );

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
