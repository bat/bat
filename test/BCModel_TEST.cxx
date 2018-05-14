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
    // mult := flag for multivariate running
    void multiple_runs(bool mult) const
    {
        unsigned N = 10000;

        unsigned n_par = 4;

        std::string prefix = mult ? "mvar_" : "fact_";

        {
            GaussModel m(prefix + "mult_run", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));
            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            for (unsigned c = 0; c < m.GetNChains(); ++c)
                TEST_CHECK_EQUAL(m.GetChainState(c).iteration / (mult ? 1 : m.GetNChains()), 2 * N);

            // Force pre-run to be rerun
            m.SetFlagPreRun(true);
            m.MarginalizeAll(BCIntegrate::kMargMetropolis);
            for (unsigned c = 0; c < m.GetNChains(); ++c)
                TEST_CHECK_EQUAL(m.GetChainState(c).iteration / (mult ? 1 : m.GetNChains()), N);
        }

        // fix a parameter
        {
            GaussModel m(prefix + "mult_run_fix_par", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            m.GetParameter(1).Fix(0);

            TEST_CHECK_THROWS(std::runtime_error, m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // check works with a fixed parameter and fails when unfixing it
        {
            unsigned n_fix = 1;

            GaussModel m(prefix + "mult_run_fixed_par", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);
            m.GetParameter(n_fix).Fix(0);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));
            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            m.GetParameter(n_fix).Unfix();

            TEST_CHECK_THROWS(std::runtime_error, m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }


        // add a parameter
        {
            GaussModel m(prefix + "mult_run_add_par", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            m.AddParameter("new_parameter", 0, 1);

            TEST_CHECK_THROWS(std::runtime_error, m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // change the number of chains
        {
            GaussModel m(prefix + "mult_run_change_nchains", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            m.SetNChains(2 * m.GetNChains());

            TEST_CHECK_THROWS(std::runtime_error, m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // change from proposal function
        {
            GaussModel m(prefix + "mult_run_change_proposal", n_par);
            m.SetPrecision(BCEngineMCMC::kMedium);
            m.SetNIterationsRun(N);
            m.SetProposeMultivariate(mult);

            TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

            m.SetProposeMultivariate(!m.GetProposeMultivariate());

            TEST_CHECK_THROWS(std::runtime_error, m.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

    }


    // check that
    // 1) one can continue running a marginalization
    // 2) continued marginalization fails if the model is changed
    // mult := flag for multivariate running
    void continued_running(bool mult) const
    {
        unsigned N = 10000;

        unsigned n_par = 4;
        unsigned n_fixed = 2;
        unsigned n_fix = 3;

        std::string prefix = mult ? "mvar_" : "fact_";

        GaussModel m(prefix + "first_run", n_par);
        m.SetPrecision(BCEngineMCMC::kMedium);
        m.SetNIterationsRun(N);
        m.SetProposeMultivariate(mult);
        m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");
        m.GetParameter(n_fixed).Fix(0.15);

        TEST_CHECK_NO_THROW(m.MarginalizeAll(BCIntegrate::kMargMetropolis));

        // test it works
        {
            GaussModel mm(prefix + "cont_run", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            // check that the proposal function is correctly chosen
            TEST_CHECK_EQUAL(mm.GetProposeMultivariate(), m.GetProposeMultivariate());
            for (unsigned p = 0; p < m.GetNParameters(); ++p) {
                TEST_CHECK_EQUAL(mm.GetParameter(p).Fixed(), m.GetParameter(p).Fixed());
                if (m.GetParameter(p).Fixed())
                    TEST_CHECK_EQUAL(mm.GetParameter(p).GetFixedValue(), m.GetParameter(p).GetFixedValue());
            }
            TEST_CHECK_NO_THROW(mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // add a parameter
        {
            GaussModel mm(prefix + "cont_run_add_par", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            mm.AddParameter("new_parameter", 0, 1);
            TEST_CHECK_THROWS(std::runtime_error, mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // change the number of chains
        {
            GaussModel mm(prefix + "cont_run_change_nchains", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            mm.SetNChains(2 * m.GetNChains());
            TEST_CHECK_THROWS(std::runtime_error, mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // change from proposal function
        // switching from factorized -> multivariate does not break, but is not recommended
        {
            GaussModel mm(prefix + "cont_run_change_proposal", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            mm.SetProposeMultivariate(!m.GetProposeMultivariate());
            TEST_CHECK_THROWS(std::runtime_error, mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // fix a parameter
        {
            GaussModel mm(prefix + "cont_run_fix_par", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            mm.GetParameter(n_fix).Fix(0.);
            TEST_CHECK_THROWS(std::runtime_error, mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }

        // unfix a parameter
        {
            GaussModel mm(prefix + "cont_run_unfix_par", n_par);
            mm.PrepareToContinueMarginalization(m.GetSafeName() + "_mcmc.root");
            mm.GetParameter(n_fixed).Unfix();
            TEST_CHECK_THROWS(std::runtime_error, mm.MarginalizeAll(BCIntegrate::kMargMetropolis));
        }
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
        TEST_SECTION("multiple marginalization", {
            multiple_runs(true);
            multiple_runs(false);
        });
        TEST_SECTION("continued marginalization", {
            continued_running(true);
            continued_running(false);
        });
    }
} bcmodel_test;
