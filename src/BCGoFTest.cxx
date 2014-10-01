/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCGoFTest.h"

#include "BCDataPoint.h"
#include "BCDataSet.h"
#include "BCLog.h"
#include "BCParameter.h"

#include <TH1D.h>
#include <TString.h>

// ---------------------------------------------------------

BCGoFTest::BCGoFTest(const char* name) : BCModel(name)
{
   // set original data set to zero
   fTemporaryDataSet = 0;

   // set test mode to zero
   fTestModel = 0;

   // reset pvalue and counter
   fPValue = 0;
   fPValueAbove = 0;
   fPValueBelow = 0;

   // reset loglikelihood and range
   fLogLikelihood = 0;
   fLogLikelihoodMin = 1e99;
   fLogLikelihoodMax = -1e99;

   // define new histogram
   fHistogramLogProb = 0;

   // set defaults for the MCMC
   MCMCSetNChains(5);
   MCMCSetNIterationsMax(100000);
   MCMCSetNIterationsRun(2000);
}

// ---------------------------------------------------------

BCGoFTest::~BCGoFTest()
{
   // restore original data set

   // get number of data points and values
   int ndatapoints = fTemporaryDataSet->GetNDataPoints();
   int ndatavalues = fTemporaryDataSet->GetDataPoint(0)->GetNValues();

   for (int i = 0; i < ndatapoints; ++i)
      for (int j = 0; j < ndatavalues; ++j)
         fTestModel->GetDataSet()->GetDataPoint(i)->SetValue(j, fTemporaryDataSet->GetDataPoint(i)->GetValue(j));

   // restore data point limits
   for (unsigned int i = 0; i < GetNParameters(); ++i)
      fTestModel->SetDataBoundaries(
            fMapDataValue[i],
            GetParameter(i)->GetLowerLimit(),
            GetParameter(i)->GetUpperLimit());

   // delete temporary data set
   delete fTemporaryDataSet;
}

// ---------------------------------------------------------

double BCGoFTest::LogLikelihood(const std::vector<double> & parameters)
{
   // set the original data set to the new parameters
   for (int i = 0; i < int(parameters.size()); ++i)
      fTestModel->GetDataSet()->GetDataPoint(fMapDataPoint[i])->SetValue(fMapDataValue[i], parameters.at(i));

   // calculate likelihood at the point of the original parameters
   double loglikelihood = fTestModel->LogLikelihood(fDataSet->GetDataPoint(0)->GetValues());

   // return likelihood
   return loglikelihood;
}

// ---------------------------------------------------------

void BCGoFTest::MCMCUserIterationInterface()
{
   int nchains = MCMCGetNChains();

   for (int i = 0; i < nchains; ++i)
   {
      // get likelihood at the point of the original parameters
      double loglikelihood = MCMCGetLogProbx(i);

      // calculate pvalue
      if (loglikelihood < fLogLikelihood)
         fPValueBelow++;
      else
         fPValueAbove++;

      // if histogram exists already, then fill it ...
      if (fHistogramLogProb)
         fHistogramLogProb->Fill(loglikelihood);
      // ...otherwise find range
      else
      {
         if (loglikelihood > fLogLikelihoodMax)
            fLogLikelihoodMax = loglikelihood;
         else if (loglikelihood < fLogLikelihoodMin)
            fLogLikelihoodMin = loglikelihood;
      }
   }
}

// ---------------------------------------------------------

int BCGoFTest::SetTestPoint(std::vector<double> parameters)
{
   // check if the boundaries of the original data set exist.
   if (!fTestModel->GetFlagBoundaries())
   {
      BCLog::OutError("BCGoFTest::SetTestDataPoint : Boundaries of the original data set are not defined.");
      return 0;
   }

   // reset histogram
   if (fHistogramLogProb)
   {
      delete fHistogramLogProb;
      fHistogramLogProb = 0;
   }

   // reset variables
   fPValue = 0;
   fPValueAbove = 0;
   fPValueBelow = 0;

   // create temporary data set ...
   fTemporaryDataSet = new BCDataSet();

   // ... and fill with the original one

   // get number of data points and values
   int ndatapoints = fTestModel->GetDataSet()->GetNDataPoints();
   int ndatavalues = fTestModel->GetDataSet()->GetDataPoint(0)->GetNValues();

   for (int i = 0; i < ndatapoints; ++i)
   {
      BCDataPoint * dp = new BCDataPoint(fTestModel->GetDataSet()->GetDataPoint(i)->GetValues());
      fTemporaryDataSet->AddDataPoint(dp);
   }

   // clear maps
   fMapDataPoint.clear();
   fMapDataValue.clear();

   int counter = 0;

   // remove parameters, but doesn't clear up memory
   fParameters = BCParameterSet();

   // loop through data points and values
   for (int i = 0; i < ndatapoints; ++i)
      for (int j = 0; j < ndatavalues; ++j)
      {
        // debugKK
        // needs to be fixed
        //         if (fTestModel->GetFixedDataAxis(j))
        //            continue;

         // add parameter to this model
         std::string parName = Form("parameter_%i", counter);
         AddParameter(
               parName.c_str(),
               fTestModel->GetDataPointLowerBoundary(j),
               fTestModel->GetDataPointUpperBoundary(j));

         // add another element to the maps
         fMapDataPoint.push_back(i);
         fMapDataValue.push_back(j);

         // increase counter
         counter ++;
      }

   // check if there are any non-fixed data values left
   if (counter == 0)
   {
      BCLog::OutError("BCGoFTest::SetTestDataPoint : No non-fixed data values left.");
      return 0;
   }

   // create a new data set containing the vector of parameters which
   // are to be tested
   BCDataPoint * datapoint = new BCDataPoint(parameters);
   BCDataSet * dataset = new BCDataSet();
   dataset->AddDataPoint(datapoint);

   // calculate likelihood of the original data set
   fLogLikelihood = fTestModel->LogLikelihood(parameters);

   // if data set has been set before, delete
   if (fDataSet)
      delete fDataSet;

   // set data set of this model
   fDataSet = dataset;

   // put proper range to new data set
   for (int i = 0; i < int(parameters.size()); ++i)
      SetDataBoundaries(
            i,
            fTestModel->GetParameter(i)->GetLowerLimit(),
            fTestModel->GetParameter(i)->GetUpperLimit());

   return 1;
}

// ---------------------------------------------------------

double BCGoFTest::GetCalculatedPValue(bool flag_histogram)
{
   // set histogram point to null
   fHistogramLogProb = 0;

   if (flag_histogram)
   {
      // perform first run to obtain limits for the log(likelihood)
      MarginalizeAll();

      // create histogram
      double D = fLogLikelihoodMax - fLogLikelihoodMin;
      fHistogramLogProb = new TH1D(Form("hist_%s_logprob", GetName().data()), ";ln(prob);N", 100, fLogLikelihoodMin - 0.1*D, fLogLikelihoodMax + 0.1*D);
      fHistogramLogProb->SetStats(kFALSE);
   }
   else
   {
   }

   // run MCMC
   MarginalizeAll();

   // check for convergence
   if (MCMCGetNIterationsConvergenceGlobal() < 0.)
   {
      BCLog::OutDetail(" --> MCMC did not converge in evaluation of the p-value.");
      return -1;
   }

   // calculate p-value
   fPValue = double(fPValueBelow) / double(fPValueBelow + fPValueAbove);

   // return p-value
   return fPValue;
}

// ---------------------------------------------------------
