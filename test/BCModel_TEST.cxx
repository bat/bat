/*
* Copyright (C) 2007-2018, the BAT core developer team
* All rights reserved.
*
* For the licensing terms see doc/COPYING.
* For documentation see http://mpp.mpg.de/bat
*/

#include <config.h>
#include "GaussModel.h"
#include "test.h"

#include <TH1.h>
#include <TH2.h>

#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCParameter.h>

#include <stdexcept>

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

    // check that
    // 1) one can call MarginalizeAll() several times
    // 2) one can force the pre-run to be rerun
    // 3) MarginalizeAll() fails if the model is changed
    void multiple_runs() const
    {
        GaussModel m("mult_run", 4);
        unsigned N = 10000;
        m.SetPrecision(BCEngineMCMC::kMedium);
        m.SetNIterationsRun(N);

        TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

        for (unsigned c = 0; c < m.GetNChains(); ++c)
            TEST_CHECK_EQUAL(m.GetChainState(c).iteration, 2 * N);

        // Force pre-run to be rerun
        m.SetFlagPreRun(true);
        m.MarginalizeAll(BCIntegrate::kMargMetropolis);
        for (unsigned c = 0; c < m.GetNChains(); ++c)
            TEST_CHECK_EQUAL(m.GetChainState(c).iteration, N);

        // fix a parameter
        GaussModel m1("mult_run_1", 4);
        m1.MarginalizeAll(BCIntegrate::kMargMetropolis);
        m1.GetParameter(0).Fix(0);
        TEST_CHECK_THROWS(std::runtime_error, m1.MarginalizeAll(BCIntegrate::kMargMetropolis));

        // add a parameter
        GaussModel m2("mult_run_2", 4);
        m2.MarginalizeAll(BCIntegrate::kMargMetropolis);
        m2.AddParameter("new_parameter", 0, 1);
        TEST_CHECK_THROWS(std::runtime_error, m2.MarginalizeAll(BCIntegrate::kMargMetropolis));

        // change the number of chains
        GaussModel m3("mult_run_3", 4);
        m3.MarginalizeAll(BCIntegrate::kMargMetropolis);
        m3.SetNChains(2 * m.GetNChains());
        TEST_CHECK_THROWS(std::runtime_error, m3.MarginalizeAll(BCIntegrate::kMargMetropolis));

        // change from multivariate to factorized proposal
        GaussModel m4("mult_run_4", 4);
        m4.SetProposeMultivariate(true);
        m4.MarginalizeAll(BCIntegrate::kMargMetropolis);
        m4.SetProposeMultivariate(false);
        TEST_CHECK_THROWS(std::runtime_error, m4.MarginalizeAll(BCIntegrate::kMargMetropolis));

        // change from factorized to multivariate proposal
        GaussModel m5("mult_run_5", 4);
        m5.SetProposeMultivariate(false);
        m5.MarginalizeAll(BCIntegrate::kMargMetropolis);
        m5.SetProposeMultivariate(true);
        TEST_CHECK_THROWS(std::runtime_error, m5.MarginalizeAll(BCIntegrate::kMargMetropolis));
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
        m.SetFlagPreRun(true);
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
        m.SetFlagPreRun(true);
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

    void compare_hist(const TH1* ref, const TH1* target) const
    {
        TEST_CHECK(target != ref);
        TEST_CHECK(target != NULL);
        TEST_CHECK_EQUAL(target->GetNbinsX(), ref->GetNbinsX());
        for (int i = 0; i <= ref->GetNbinsX(); ++i) {
            TEST_CHECK_EQUAL(target->GetBinContent(i),
                             ref->GetBinContent(i));
        }
    }

    void compare_hist(const TH2* ref, const TH2* target) const
    {
        TEST_CHECK(target != ref);
        TEST_CHECK(target != NULL);
        TEST_CHECK_EQUAL(target->GetNbinsX(), ref->GetNbinsX());
        TEST_CHECK_EQUAL(target->GetNbinsY(), ref->GetNbinsY());
        for (int i = 0; i <= ref->GetNbinsX(); ++i) {
            for (int j = 0; j <= ref->GetNbinsY(); ++j) {
                TEST_CHECK_EQUAL(target->GetBinContent(i, j),
                                 ref->GetBinContent(i, j));
            }
        }
    }

    void copy() const
    {
        GaussModel m("copy", 3);

        // add a data set (not copied!) even though it is not used
        BCDataSet data(1);
        m.SetDataSet(&data);
        TEST_CHECK_EQUAL(m.GetDataSet(), &data);

        // use non-default values
        m.SetNChains(m.GetNChains() + 1);

        m.MarginalizeAll();

        // histograms are dynamically created, they are the criticial part
        // during copying.
        TH1* const h0 = m.GetMarginalizedHistogram(0);
        TH1* const h1 = m.GetMarginalizedHistogram(1);
        TH2* const h01 = m.GetMarginalizedHistogram(0, 1);

        TEST_SECTION("copy ctor",

                     GaussModel m2(m);

                     // non-default values should be taken over
                     TEST_CHECK_EQUAL(m2.GetNChains(), m.GetNChains());

                     // data set shared
                     TEST_CHECK_EQUAL(m2.GetDataSet(), &data);

                     // Original values should be untouched
                     TEST_CHECK_EQUAL(h0, m.GetMarginalizedHistogram(0));
                     TEST_CHECK_EQUAL(h1, m.GetMarginalizedHistogram(1));
                     TEST_CHECK_EQUAL(h01, m.GetMarginalizedHistogram(0, 1));

                     compare_hist(h0, m2.GetMarginalizedHistogram(0));
                     compare_hist(h1, m2.GetMarginalizedHistogram(1));
                     compare_hist(h01, m2.GetMarginalizedHistogram(0, 1));
                    );

        TEST_SECTION("assignment operator",

                     GaussModel m3("assignment", 1);
                     TEST_CHECK_EQUAL(m3.GetMarginalizedHistogram(0), NULL);
                     TEST_CHECK_EQUAL(m3.GetDataSet(), NULL);

                     m3 = m;

                     TEST_CHECK_EQUAL(m3.GetName(), "copy");
                     TEST_CHECK_EQUAL(m3.GetDataSet(), &data);
                     TEST_CHECK_EQUAL(m3.GetNChains(), m.GetNChains());

                     compare_hist(h0, m3.GetMarginalizedHistogram(0));
                     compare_hist(h1, m3.GetMarginalizedHistogram(1));
                     compare_hist(h01, m3.GetMarginalizedHistogram(0, 1));
                    );
    }

    virtual void run() const
    {
        storing();
        fixing(true);
        fixing(false);
        deltaPrior();
        copy();
        TEST_SECTION("multiple marginalize", {
            multiple_runs();
        });
    }
} bcmodel_test;
