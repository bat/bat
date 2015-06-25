/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "test.h"
#include "GaussModel.h"

using namespace test;

namespace
{
GaussModel* gauss_check(bool multivariate_proposal)
{
    /* set up model */
    std::string name("BCEngineMCMC_TEST-model");
    if (multivariate_proposal)
        name += "multivariate";
    static const unsigned dim = 5;
    GaussModel& m = *new GaussModel(name.c_str(), dim);

    m.MCMCSetNChains(2);
    m.MCMCSetNIterationsPreRunCheck(500);
    m.MCMCSetNIterationsClearConvergenceStats(20000);
    m.MCMCSetNIterationsPreRunMax(1000000);
    m.MCMCSetNIterationsPreRunMin(1000);
    m.MCMCSetNIterationsRun(100000);
    m.MCMCSetMinimumEfficiency(0.15);
    m.MCMCSetMaximumEfficiency(0.35);

    m.MCMCSetCorrectRValueForSamplingVariability(true);
    m.MCMCSetMultivariateProposalFunction(multivariate_proposal);

    m.MCMCSetRandomSeed(23062015);

    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    /* analyze statistics */

    // no accidental calls to likelihood
    unsigned ncalls = m.MCMCGetNChains() * (m.MCMCGetNIterationsConvergenceGlobal() + m.MCMCGetNIterationsRun());
    ncalls *= multivariate_proposal ? 1 : dim;
    TEST_CHECK_EQUAL(m.Calls(), ncalls);
    //            m.PrintAllMarginalized(std::string(__FILE__) + "_gauss.pdf");
    const BCEngineMCMC::MCMCStatistics& s = m.MCMCGetStatistics();

    // samples are counted for the total of chains
    TEST_CHECK_EQUAL(s.n_samples, m.MCMCGetNIterationsRun() * m.MCMCGetNChains());
    // mean etc. are per parameter
    TEST_CHECK_EQUAL(s.mean.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.variance.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.covariance.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.minimum.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.maximum.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.mode.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(s.efficiency.size(), m.GetNParameters());
    for (unsigned i = 0; i < m.GetNParameters(); ++i) {
        /* compare to values set inside gauss model' likelihood */

        // error on the mean from independent samples
        static const double sigma = 2;
        static const double var = sigma * sigma;
        const double bestError = sqrt(var * m.MCMCGetNIterationsRun()) / m.MCMCGetNIterationsRun();

        // we don't know the effective sample size but should use it here
        // instead of an arbitrary factor
        static const double ess = 10;
        TEST_CHECK_NEARLY_EQUAL(s.mean[i], 0, ess * bestError);

        // variance estimate much worse than mean estimate
        TEST_CHECK_RELATIVE_ERROR(s.variance[i], var, 0.3);

        TEST_CHECK_EQUAL(s.covariance[i].size(), m.GetNParameters());

        TEST_CHECK_EQUAL(s.covariance[i].size(), m.GetNParameters());
        for (unsigned j = 0; j < m.GetNParameters(); ++j) {
            TEST_CHECK_NEARLY_EQUAL(s.covariance[i][j], i == j ? var : 0.0, 0.85);
        }

        // expect to be within 3.5-5 sigma for 1e5 samples
        static const double xmin = 3.5 * sigma;
        static const double xmax = 5 * sigma;
        TEST_CHECK(s.maximum[i] > xmin);
        TEST_CHECK(s.maximum[i] < xmax);
        TEST_CHECK(s.minimum[i] < -xmin);
        TEST_CHECK(s.minimum[i] > -xmax);

        // mcmc not a great mode finder => large uncertainty
        TEST_CHECK_NEARLY_EQUAL(s.mode[i], 0, 0.4);

        // 23.8% is "optimal" acceptance rate for Gaussian target in high dimensions
        // and Gaussian proposal
        TEST_CHECK_NEARLY_EQUAL(s.efficiency[i], 0.238, 0.1);
    }
    return &m;
}
}

class ConvergenceTest :
    public TestCase
{
public:
    ConvergenceTest() :
        TestCase("Convergence test")
    {
        BCLog::SetLogLevelScreen(BCLog::debug);
    }

    virtual void run() const
    {
        TEST_SECTION("factorized",
                     GaussModel* m = ::gauss_check(false);
                     delete m;
                    );
        TEST_SECTION("multivariate",
                     GaussModel* m = ::gauss_check(true);
                     delete m;
                    );
    }
} convergenceTest;

// Local Variables:
// compile-command: "make -C .. check TESTS= && (./BCEngineMCMC_TEST || cat test-suite.log)"
// End:
