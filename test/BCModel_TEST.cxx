/*
 * Copyright (C) 2013, Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <config.h>
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

   /*!
    * check that the right marginal distributions are defined
    * @par fixed index of the fixed parameter
    */
   static void count_marginals(BCModel & m, unsigned fixed)
   {
      m.MarginalizeAll();
      const unsigned nplots = m.PrintAllMarginalized(BAT_TESTDIR "BCModel_TEST.pdf");
      // 1D + 2D
      TEST_CHECK_EQUAL(nplots, (m.GetNParameters() - 1) +
                       (m.GetNParameters() - 2) * (m.GetNParameters() - 1) / 2);
      for (unsigned i = 0 ; i < m.GetNParameters() ; ++i)
         TEST_CHECK(bool(m.GetMarginalized(i)) xor (i == fixed));
      for (unsigned i = 0 ; i < m.GetNParameters() ; ++i)
         for (unsigned j = i + 1 ; j < m.GetNParameters() ; ++j)
            TEST_CHECK(bool(m.GetMarginalized(i, j)) xor ( (i == fixed) xor (j == fixed)));
   }

   // turn on/off parameter storing
   void storing() const
   {
      const unsigned fixed = 1;
      GaussModel m("Turn off filling for par 1", 4);
      m.MCMCSetPrecision(BCEngineMCMC::kLow);
      m.GetParameter(fixed)->FillHistograms(false);
      count_marginals(m, fixed);
   }

   // check parameter fixing with mcmc
   // 1) rerun one model with fixed and unfixed parameters
   // 2) count stored distributions
   // 3) marginal mean and variance should be right
   // 4) fixed parameters require less likelihood calls
   void fixing() const
   {
      GaussModel m("Fix", 4);
      m.MCMCSetPrecision(BCEngineMCMC::kMedium);
      m.MCMCSetNIterationsRun(30000);
      m.MCMCSetRandomSeed(235);

      static const double eps = 5e-2;

      /* fix first parameter */
      int fixed = 0;
      m.GetParameter(fixed)->Fix(0.);
      count_marginals(m, fixed);

      // gaussian around zero with width two
      TEST_CHECK_NEARLY_EQUAL(0, m.GetMarginalized(2)->GetMean(), eps);
      TEST_CHECK_RELATIVE_ERROR(2, m.GetMarginalized(2)->GetSTD(), eps);

      // make sure that fixed parameters really are skipped
      const unsigned long nCalls = m.Calls();

      /* so run again without fixing */
      m.GetParameter(0)->Unfix();
      m.MarginalizeAll();

      // should be extra 33% calls
      const unsigned long unfixedCalls = m.Calls() - nCalls;

      TEST_CHECK_RELATIVE_ERROR(double(nCalls) / unfixedCalls, 0.75, 1e-2);

      /* then fix last parameter */
      fixed = 3;
      m.GetParameter(fixed)->Fix(0.23);
      count_marginals(m, fixed);

      // gaussian around zero with width two
      TEST_CHECK_NEARLY_EQUAL(0, m.GetMarginalized(0u)->GetMean(), eps);
      TEST_CHECK_RELATIVE_ERROR(2, m.GetMarginalized(0u)->GetSTD(), eps);
   }

   void deltaPrior() const
   {
      GaussModel m("set delta prior for par 1", 4);
      m.MCMCSetPrecision(BCEngineMCMC::kMedium);
      const int fixed = 0;
      m.SetPriorDelta(fixed, 0.);
      count_marginals(m, fixed);
   }

   virtual void run() const
   {
      storing();
      fixing();
      deltaPrior();
   }
} bcmodel_test;
