/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "test.h"
#include "GaussModel.h"

#include <limits>

using namespace test;

namespace
{
void check_efficiency(const GaussModel& m, double efficiency)
{
    if (efficiency < m.GetMinimumEfficiency())
        TEST_CHECK_FAILED(stringify(efficiency) + " less than required " + stringify(m.GetMinimumEfficiency()));
    if (efficiency > m.GetMaximumEfficiency())
        TEST_CHECK_FAILED(stringify(efficiency) + " larger than required " + stringify(m.GetMaximumEfficiency()));
}

GaussModel* gauss_check(bool multivariate_proposal, double fix = std::numeric_limits<double>::infinity())
{
    /* set up model */
    const bool fixLast = (fix != std::numeric_limits<double>::infinity());
    std::string name("BCEngineMCMC_TEST-model");
    if (multivariate_proposal)
        name += ", multivariate";
    if (fixLast)
        name += ", fix last";
    static const unsigned ntotal = 5;
    GaussModel& m = *new GaussModel(name.c_str(), ntotal);

    // Fix the last parameter so comparison below is simpler to code
    // keep track of free and total number of parameters
    unsigned nfree = ntotal;
    if (fixLast) {
        m.GetParameter(ntotal - 1).Fix(fix);
        --nfree;
    }
    TEST_CHECK_EQUAL(m.GetNParameters(), ntotal);
    TEST_CHECK_EQUAL(m.GetNFreeParameters(), nfree);
    TEST_CHECK_EQUAL(m.GetNFixedParameters(), ntotal - nfree);

    m.SetNChains(2);
    m.SetNIterationsPreRunCheck(500);
    m.SetPreRunCheckClear(40);
    m.SetNIterationsPreRunMax(20000);
    m.SetNIterationsPreRunMin(1000);
    m.SetNIterationsRun(100000);

    m.SetMinimumEfficiency(0.10);
    m.SetMaximumEfficiency(0.30);

    m.SetCorrectRValueForSamplingVariability(true);
    m.SetProposeMultivariate(multivariate_proposal);

    m.SetRandomSeed(23062015);

    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    /* convergence */

    // the chains should converge before hitting the maximum prerun length for this simple problem
    TEST_CHECK(m.GetNIterationsConvergenceGlobal() > 0);
    TEST_CHECK(unsigned(m.GetNIterationsConvergenceGlobal()) >= m.GetNIterationsPreRunMin());
    TEST_CHECK(unsigned(m.GetNIterationsConvergenceGlobal()) < m.GetNIterationsPreRunMax());
    TEST_CHECK(m.GetNIterationsPreRun() >= m.GetNIterationsPreRunMin());
    TEST_CHECK(m.GetNIterationsPreRun() <= m.GetNIterationsPreRunMax());

    /* analyze statistics */

    // no accidental calls to likelihood
    unsigned ncalls = m.GetNChains() * (m.GetNIterationsConvergenceGlobal() + m.GetNIterationsRun());
    ncalls *= multivariate_proposal ? 1 : nfree;
    TEST_CHECK_EQUAL(m.Calls(), ncalls);

    // target variance
    const double var = m.sigma() * m.sigma();

    // check each chain
    const std::vector<BCEngineMCMC::Statistics>& s = m.GetStatisticsVector();
    for (unsigned j = 0; j < s.size(); ++j) {
        BCLog::OutSummary(Form("Checking chain %d -->", j));

        TEST_CHECK_EQUAL(s[j].n_samples, m.GetNIterationsRun());
        // mean etc. are per parameter
        TEST_CHECK_EQUAL(s[j].mean.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].variance.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].covariance.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].minimum.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].maximum.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].mode.size(), m.GetNParameters());
        TEST_CHECK_EQUAL(s[j].efficiency.size(), m.GetNParameters());
        for (unsigned i = 0; i < m.GetNFreeParameters(); ++i) {
            /* compare to values set inside gauss model' likelihood */
            // error on the mean from independent samples
            const double bestError = sqrt(var * m.GetNIterationsRun()) / m.GetNIterationsRun();

            // we don't know the effective sample size but should use it here
            // instead of an arbitrary factor
            static const double ess = 15;
            TEST_CHECK_NEARLY_EQUAL(s[j].mean[i], 0, ess * bestError);

            // variance estimate much worse than mean estimate
            TEST_CHECK_RELATIVE_ERROR(s[j].variance[i], var, 0.75);
            TEST_CHECK_RELATIVE_ERROR(s[j].covariance[i][i], var, 0.75);

            TEST_CHECK_EQUAL(s[j].covariance[i].size(), m.GetNParameters());

            // check consistency of covariance[i==k] and variance[i]
            TEST_CHECK_NEARLY_EQUAL(s[j].covariance[i][i], s[j].variance[i], 0.01);

            TEST_CHECK_EQUAL(s[j].covariance[i].size(), m.GetNParameters());

            // compare to fixed parameter, too
            for (unsigned k = 0; k < m.GetNParameters(); ++k)
                if (i != k)
                    TEST_CHECK_NEARLY_EQUAL(s[j].covariance[i][k], 0.0, 2);

            // mcmc not a great mode finder => large uncertainty
            TEST_CHECK_NEARLY_EQUAL(s[j].mode[i], 0, 1);

            if (!multivariate_proposal)
                check_efficiency(m, s[j].efficiency[i]);
        }

        // check fixed parameter
        if (fixLast) {
            TEST_CHECK_EQUAL(s[j].mean.back(), fix);
            TEST_CHECK_EQUAL(s[j].variance.back(), 0.0);
            TEST_CHECK_EQUAL(s[j].minimum.back(), fix);
            TEST_CHECK_EQUAL(s[j].maximum.back(), fix);
            TEST_CHECK_EQUAL(s[j].mode.back(), fix);
            TEST_CHECK_EQUAL(s[j].efficiency.back(), 0.0);

            continue;
        }

        if (multivariate_proposal)
            check_efficiency(m, s[j].efficiency.front());
        BCLog::OutSummary(Form("<-- Checked chain %d", j));
    }

    // check all-chains statistics
    const BCEngineMCMC::Statistics& S = m.GetStatistics();
    TEST_CHECK_EQUAL(S.n_samples, m.GetNIterationsRun() * m.GetNChains());
    // mean etc. are per parameter
    TEST_CHECK_EQUAL(S.mean.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.variance.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.covariance.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.minimum.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.maximum.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.mode.size(), m.GetNParameters());
    TEST_CHECK_EQUAL(S.efficiency.size(), m.GetNParameters());
    for (unsigned i = 0; i < m.GetNFreeParameters(); ++i) {
        /* compare to values set inside gauss model' likelihood */

        // error on the mean from independent samples
        const double bestError = sqrt(var * m.GetNIterationsRun()) / m.GetNIterationsRun();

        // we don't know the effective sample size but should use it here
        // instead of an arbitrary factor
        static const double ess = 15;
        TEST_CHECK_NEARLY_EQUAL(S.mean[i], 0, ess * bestError);

        // variance estimate much worse than mean estimate
        TEST_CHECK_RELATIVE_ERROR(S.variance[i], var, 0.3);

        TEST_CHECK_EQUAL(S.covariance[i].size(), m.GetNParameters());

        // check consistency of covariance[i==k] and variance[i]
        TEST_CHECK_NEARLY_EQUAL(S.covariance[i][i], S.variance[i], 0.01);

        TEST_CHECK_EQUAL(S.covariance[i].size(), m.GetNParameters());
        for (unsigned j = 0; j < m.GetNParameters(); ++j) {
            TEST_CHECK_NEARLY_EQUAL(S.covariance[i][j], i == j ? var : 0.0, 0.85);
        }

        // don't expect to fall outside 6 sigma
        static const double cut = 6;
        TEST_CHECK(S.maximum[i] > 0);
        TEST_CHECK(S.maximum[i] < cut);
        TEST_CHECK(S.minimum[i] < 0);
        TEST_CHECK(S.minimum[i] > -cut);

        // mcmc not a great mode finder => large uncertainty
        TEST_CHECK_NEARLY_EQUAL(S.mode[i], 0, 0.4);

        // 23.8% is "optimal" acceptance rate for Gaussian target in high dimensions
        // and Gaussian proposal
        if (!multivariate_proposal)
            check_efficiency(m, S.efficiency[i]);
    }
    if (multivariate_proposal)
        check_efficiency(m, S.efficiency.front());

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
        TEST_SECTION("multivariate, fix one",
                     GaussModel* m = ::gauss_check(true, 1.1);
                     delete m;
                    );
    }
} convergenceTest;

#if 0
class RValueTest :
    public TestCase
{
public:
    RValueTest() :
        TestCase("rvalue_test")
    {
    }

    virtual void run() const
    {
        static const double eps = 1e-14;
        static const bool strict = true;
        static const bool relaxed = false;

        // R-value calculation checked against implementation in EOS
        {
            std::vector<double> chain_means(3);
            chain_means[0] = 4.2;
            chain_means[1] = 4.25;
            chain_means[2] = 4.22;
            std::vector<double> chain_variances(3);
            chain_variances[0] = 0.1;
            chain_variances[1] = 0.15;
            chain_variances[2] = 0.19;

            unsigned points = 500;

            // strict always larger than relaxed
            TEST_CHECK_RELATIVE_ERROR(BCEngineMCMC::RValue(chain_means, chain_variances, points, relaxed), 1.0011584199407115, eps);
            TEST_CHECK_RELATIVE_ERROR(BCEngineMCMC::RValue(chain_means, chain_variances, points, strict), 1.0176292831481546, eps);

            // for more points visited, R-value increases
            points *= 3;

            TEST_CHECK_RELATIVE_ERROR(BCEngineMCMC::RValue(chain_means, chain_variances, points, relaxed), 1.0018240939164496, eps);
            TEST_CHECK_RELATIVE_ERROR(BCEngineMCMC::RValue(chain_means, chain_variances, points, strict), 1.0183054631320092, eps);
        }
    }
} rvalue_test;
#endif

// Local Variables:
// compile-command: "make check TESTS= && (./BCEngineMCMC.TEST || cat test-suite.log)"
// End:
