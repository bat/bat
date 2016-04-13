/*
 * Copyright (C) 2013, Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <config.h>
#include "GaussModel.h"
#include "test.h"

#include <TH1.h>

#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCParameter.h>

using namespace test;

class BCModelTest :
    public TestCase
{
public:
    BCModelTest():
        TestCase("BCModel")
    {
    }

    /*!
     * check that the right marginal distributions are defined
     * @par fixed index of the fixed parameter
     */
    static void count_marginals(BCModel& m, unsigned fixed)
    {
        m.MarginalizeAll();
        const unsigned nplots = m.PrintAllMarginalized(BAT_TESTDIR "BCModel_TEST.pdf");
        // 1D + 2D
        TEST_CHECK_EQUAL(nplots, (m.GetNParameters() - 1) + (m.GetNParameters() - 2) * (m.GetNParameters() - 1) / 2);
        for (unsigned i = 0 ; i < m.GetNParameters() ; ++i)
            TEST_CHECK(m.GetMarginalized(i).Valid() xor (i == fixed));
        for (unsigned i = 0 ; i < m.GetNParameters() ; ++i)
            for (unsigned j = i + 1 ; j < m.GetNParameters() ; ++j)
                TEST_CHECK(m.GetMarginalized(i, j).Valid() xor ((i == fixed) or (j == fixed)));
    }

    // turn on/off parameter storing
    void storing() const
    {
        const unsigned fixed = 1;
        GaussModel m("Turn off filling for par 1", 4);
        m.SetPrecision(BCEngineMCMC::kLow);
        m.GetParameter(fixed).FillHistograms(false);
        count_marginals(m, fixed);
    }

    // check parameter fixing with mcmc
    // 1) rerun one model with fixed and unfixed parameters
    // 2) count stored distributions
    // 3) marginal mean and variance should be right
    // 4) fixed parameters require less likelihood calls
    void fixing(bool multivariate) const
    {
        GaussModel m("fix first", 4);
        m.SetPrecision(BCEngineMCMC::kMedium);
        m.SetNIterationsRun(30000);
        m.SetNIterationsPreRunMax(6000);
        m.SetNIterationsPreRunMin(6000);
        m.SetRandomSeed(235);
        m.SetProposeMultivariate(multivariate);

        static const double eps = 5e-2;

        /* fix first parameter */
        int fixed = 0;
        m.GetParameter(fixed).Fix(0.);
        count_marginals(m, fixed);

        // gaussian with mean and sigma
        TEST_CHECK_NEARLY_EQUAL(m.GetMarginalizedHistogram(2)->GetMean(), m.mean(), eps);
        TEST_CHECK_RELATIVE_ERROR(m.GetMarginalizedHistogram(2)->GetRMS(), m.sigma(), eps);

        // make sure that fixed parameters really are skipped
        const unsigned long nCalls = m.Calls();

        /* so run again without fixing */
        m.SetName("all free");
        m.GetParameter(0).Unfix();
        TEST_CHECK_EQUAL(m.GetNFreeParameters(), 4);
        m.MarginalizeAll();

        if (!multivariate) {
            // should be extra 33% calls unless multivariate
            TEST_CHECK_RELATIVE_ERROR(double(nCalls) / (m.Calls() - nCalls), 0.75, 1e-2);
        }
        /* then fix last parameter */
        m.SetName("fix last");
        fixed = 3;
        m.GetParameter(fixed).Fix(0.23);
        count_marginals(m, fixed);

        // gaussian around zero with width two
        TEST_CHECK_NEARLY_EQUAL(m.GetMarginalizedHistogram(0u)->GetMean(), m.mean(), eps);
        TEST_CHECK_RELATIVE_ERROR(m.GetMarginalizedHistogram(0u)->GetRMS(), m.sigma(), eps);
    }

    void deltaPrior() const
    {
        GaussModel m("set delta prior for par 1", 4);
        m.SetPrecision(BCEngineMCMC::kMedium);
        const int fixed = 0;
        m.GetParameter(fixed).Fix(0.);
        count_marginals(m, fixed);
    }

    void copy() const
    {
        GaussModel m("copy", 1);

        // use non-default values
        m.SetNChains(m.GetNChains() + 1);

        m.MarginalizeAll();

        // basic test: no segfault
        GaussModel m2 = m;

        // non-default values should be taken over
        TEST_CHECK_EQUAL(m2.GetNChains(), m.GetNChains());
    }

    virtual void run() const
    {
        storing();
        fixing(true);
        fixing(false);
        deltaPrior();
        copy();
    }
} bcmodel_test;


// Local Variables:
// compile-command: "make -C ../ && make BCModel.TEST && ./BCModel.TEST"
// End:
