/*
 * Copyright (C) 2013, Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GaussModel.h"
#include "test.h"

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

    static void count_marginals(BCModel & m)
    {
          m.MarginalizeAll();
          const unsigned nplots = m.PrintAllMarginalized("BCModel_TEST.pdf");
          TEST_CHECK_EQUAL(nplots, 6);
          TEST_CHECK(m.GetMarginalized(0u));
          TEST_CHECK( !m.GetMarginalized(1));
          TEST_CHECK(m.GetMarginalized(2));
          TEST_CHECK( !m.GetMarginalized(0, 1));
          TEST_CHECK(m.GetMarginalized(0, 2));
          TEST_CHECK( !m.GetMarginalized(1, 2));
    }

    virtual void run() const
    {
       // turn on/off parameter storing
       {
          GaussModel m("model", 4);
          m.MCMCSetPrecision(BCEngineMCMC::kLow);
          m.GetParameter(1)->FillHistograms(false);
          count_marginals(m);
       }

       // check parameter fixing with mcmc
       {
          GaussModel m("model", 4);
          m.MCMCSetPrecision(BCEngineMCMC::kLow);
          m.GetParameter(1)->Fix(0.32);
          count_marginals(m);
       }

       // delta prior
       {
          GaussModel m("model", 4);
          m.MCMCSetPrecision(BCEngineMCMC::kLow);
          m.SetPriorDelta(1, 0.32);
          count_marginals(m);
       }
    }
} bcmodel_test;
