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

#include <TString.h>

#include <limits>

// ---------------------------------------------------------
BCGoFTest::BCGoFTest(std::string name)
    : BCModel(name)
    , fTestModel(NULL)
    , fOriginalLogLikelihood(std::numeric_limits<double>::infinity())
{
    // reset pvalue and counter
    fPValue = -1;

    // set defaults for the MCMC
    MCMCSetNChains(5);
    MCMCSetNIterationsPreRunCheck(500);
    MCMCSetNIterationsClearConvergenceStats(5000);
    MCMCSetNIterationsPreRunMax(100000);
    MCMCSetNIterationsPreRunMin(10000);
    MCMCSetNIterationsRun(2000);

    // add observable to hold distribution of log(prob)
    AddObservable("log_prob", -10, 10, "log(prob)");
}

// ---------------------------------------------------------
BCGoFTest::BCGoFTest(BCModel* model, std::string name)
    : BCModel(name)
    , fTestModel(model)
    , fOriginalLogLikelihood(std::numeric_limits<double>::infinity())
{
    if (fName.empty() and fTestModel)
        SetName(fTestModel->GetName() + "_GoFTest");

    fPValue = -1;

    // set defaults for the MCMC
    MCMCSetNChains(5);
    MCMCSetNIterationsPreRunCheck(500);
    MCMCSetNIterationsClearConvergenceStats(5000);
    MCMCSetNIterationsPreRunMax(100000);
    MCMCSetNIterationsPreRunMin(10000);
    MCMCSetNIterationsRun(2000);

    // add observable to hold distribution of log(prob)
    AddObservable("log_prob", -10, 10, "log(prob)");
}

// ---------------------------------------------------------
BCGoFTest::~BCGoFTest()
{
}

// ---------------------------------------------------------
double BCGoFTest::LogLikelihood(const std::vector<double>& parameters)
{
    // Update data set from parameter set
    for (unsigned i = 0; i < parameters.size(); ++i)
        fTestModel->GetDataSet()[fDataMap[i].first][fDataMap[i].second] = parameters[i];

    // calculate likelihood at the point of the original parameters with new data set
    double loglikelihood = fTestModel->LogLikelihood(fOriginalParameters);

    if (loglikelihood < fOriginalLogLikelihood)
        ++fPValueBelow;
    else
        ++fPValueAbove;

    // return likelihood
    return loglikelihood;
}

// ---------------------------------------------------------
void BCGoFTest::CalculateObservables(const std::vector<double>& /*pars*/)
{
    // Set log(prob) observable to log(prob)
    GetObservable(0).Value(MCMCGetLogProbx(fMCMCCurrentChain));
}

// ---------------------------------------------------------
void BCGoFTest::MCMCMetropolisPreRun()
{
    // run usual pre run
    BCEngineMCMC::MCMCMetropolisPreRun();

    // set histogramming range of log(prob) observable from pre-run results
    // size of likelihood range:
    double LL_min = MCMCGetStatistics().minimum[GetNParameters()];
    double LL_max = MCMCGetStatistics().maximum[GetNParameters()];

    // set proper range for observable: log(prob)
    GetObservable("log_prob").SetLimits(LL_min - 0.1 * (LL_max - LL_min), LL_max + 0.1 * (LL_max - LL_min));
    GetObservable("log_prob").FillHistograms(true, false);
}

// ---------------------------------------------------------
double BCGoFTest::CalculatePValue(const std::vector<double>& parameters)
{
    // check if the boundaries of the original data set exist.
    if (!fTestModel->GetDataSet()->BoundsExist()) {
        BCLog::OutError("BCGoFTest::SetTestDataPoint : Boundaries of the original data set are not defined.");
        return -1;
    }

    // reset variables
    fPValue = 0;
    fPValueAbove = 0;
    fPValueBelow = 0;

    // store test model's data set in this object
    SetDataSet(fTestModel->GetDataSet());
    // create cop of data set and place into test model
    fTestModel->SetDataSet(new BCDataSet(*fDataSet));

    // clear map & reserce space
    fDataMap.clear();
    fDataMap.reserve(fTestModel->GetDataSet()->GetNDataPoints() * fTestModel->GetDataSet()->GetNValuesPerPoint());

    // remove existing parameters
    fParameters = BCParameterSet();

    // loop through data points and values
    int counter = 0;
    for (unsigned i = 0; i < fTestModel->GetDataSet()->GetNDataPoints(); ++i)
        for (unsigned j = 0; j < fTestModel->GetDataSet()->GetNValuesPerPoint(); ++j) {
            // skip fixed data axes
            if (fTestModel->GetDataSet()->IsFixed(j))
                continue;

            // add parameter to this model corresponding to data value j of data point i
            AddParameter(Form("parameter_%i", counter++), fTestModel->GetDataSet()->GetLowerBound(j), fTestModel->GetDataSet()->GetUpperBound(j));

            // add element to the map
            fDataMap.push_back(std::make_pair(i, j));
        }

    // do not use histograms to record MCMC for parameters created above
    fParameters.FillHistograms(false);

    // check if there are any non-fixed data values left
    if (counter == 0) {
        BCLog::OutError("BCGoFTest::SetTestDataPoint : No non-fixed data values left.");
        return -1;
    }

    // store parameters to be tested
    fOriginalParameters = parameters;

    // calculate likelihood of the original data set
    fOriginalLogLikelihood = fTestModel->LogLikelihood(fOriginalParameters);

    // run MCMC
    MarginalizeAll(BCIntegrate::kMargMetropolis);

    // check for convergence
    if (MCMCGetNIterationsConvergenceGlobal() < 0.) {
        BCLog::OutDetail(" --> MCMC did not converge in evaluation of the p-value.");
        return -1;
    }

    // calculate p-value
    fPValue = 1.* fPValueBelow / (fPValueBelow + fPValueAbove);

    // restore original data set
    fTestModel->SetDataSet(fDataSet);
    // remove this object's access to data set
    SetDataSet(0);

    // return p-value
    return fPValue;
}
