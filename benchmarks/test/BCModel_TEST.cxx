/*
 * Copyright (C) 2013, Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GaussModel.h"
#include "test.h"

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

   // turn on/off parameter storing
   void storing() const
   {
      GaussModel m("Turn off filling for par 1", 4);
      m.MCMCSetPrecision(BCEngineMCMC::kLow);
      m.GetParameter(1)->FillHistograms(false);
      count_marginals(m);
   }

   // check parameter fixing with mcmc
   // 1) rerun one model with fixed and unfixed parameters
   // 2) count stored distributions
   // 3) marginal mean and variance should be right
   // 4) fixed parameters require less likelihood calls
   void fixing() const
   {
      GaussModel m("Fix par 1", 4);
      m.MCMCSetRandomSeed(235);
      m.MCMCSetPrecision(BCEngineMCMC::kLow);
      m.MCMCSetNIterationsRun(30000);
      m.GetParameter(1)->Fix(0.32);
      count_marginals(m);

      // gaussian around zero with width two
      TEST_CHECK_NEARLY_EQUAL(0, m.GetMarginalized(2)->GetMean(), 5e-2);
      TEST_CHECK_RELATIVE_ERROR(2, m.GetMarginalized(2)->GetSTD(), 5e-2);

      // make sure that fixed parameters really are skipped
      const unsigned long nCalls = m.Calls();

      // so run again without fixing
      m.GetParameter(1)->Unfix();
      m.MarginalizeAll();

      // should be extra 33% calls
      const unsigned long unfixedCalls = m.Calls() - nCalls;

      TEST_CHECK_RELATIVE_ERROR(double(nCalls) / unfixedCalls, 0.75, 1e-2);
   }

   void deltaPrior() const
   {
      GaussModel m("set delta prior for par 1", 4);
      m.MCMCSetPrecision(BCEngineMCMC::kLow);
      m.SetPriorDelta(1, 0.32);
      count_marginals(m);
   }

   virtual void run() const
   {
      storing();
      fixing();
      deltaPrior();
   }
} bcmodel_test;
