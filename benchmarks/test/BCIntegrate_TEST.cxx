/*
 * Copyright (C) 2013, Daniel Greenwald and Frederik Beaujean
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "IntegrationModel.h"
#include "test.h"
#include <config.h>

using namespace test;

class BCIntegrateTest :
public TestCase
{
public:

   struct Precision
   {
      double monteCarlo, importance, vegas, suave, divonne, cuhre;

      // Turn off checks by default
      Precision(double value = -1)
      {
         Set(value);
      }

      void Set(double value)
      {
         monteCarlo = importance = vegas = suave = divonne = cuhre = value;
      }

   };

   BCIntegrateTest() :
      TestCase("BCIntegrate test")
   {
   }

   static IntegrationModel Factory(unsigned dim, unsigned modality, unsigned complexity)
   {
      IntegrationModel m;
      m.SetDimensionality(dim);
      m.SetModality(modality);
      m.SetComplexity(complexity);
      m.PopulatePolynomialDegrees();
      m.MCMCSetRandomSeed(1346);

      return m;
   }

   static void SingleCase(unsigned dim, unsigned modality, unsigned complexity, const Precision & p)
   {
      IntegrationModel m = Factory(dim, modality, complexity);
      std::cout << "Computing single integration test in " << dim
            << " dimensions with modality = " << modality
            << " and complexity = " << complexity << ".\n";

      if (p.monteCarlo > 0)
         TEST_CHECK_RELATIVE_ERROR(m.Integrate(BCIntegrate::kIntMonteCarlo), m.Integral(), p.monteCarlo);
      if (p.importance > 0)
         TEST_CHECK_RELATIVE_ERROR(m.Integrate(BCIntegrate::kIntImportance), m.Integral(), p.importance);
#if HAVE_CUBA_H
      if (p.vegas > 0)
         TEST_CHECK_RELATIVE_ERROR(m.IntegrateCuba(BCIntegrate::kCubaVegas), m.Integral(), p.vegas);
      if (p.suave > 0)
         TEST_CHECK_RELATIVE_ERROR(m.IntegrateCuba(BCIntegrate::kCubaSuave), m.Integral(), p.suave);
      if (p.divonne> 0)
         TEST_CHECK_RELATIVE_ERROR(m.IntegrateCuba(BCIntegrate::kCubaDivonne), m.Integral(), p.divonne);
      if (p.cuhre > 0)
         TEST_CHECK_RELATIVE_ERROR(m.IntegrateCuba(BCIntegrate::kCubaCuhre), m.Integral(), p.cuhre);
#endif
   }

   virtual void run() const
   {

      {
         Precision p;
         p.monteCarlo = 1e-2;
         p.vegas = 2e-2;
         p.suave = 1e-2;
         p.divonne = 2; // fails in 1D
         p.cuhre = 2; // fails in 1D
         SingleCase(1, 1, 1, p);
      }

      {
         Precision p;
         p.divonne = p.cuhre = 2e-3;
         SingleCase(2, 1, 1, p);
      }
      // todo reactivate divonne when cuba 3.1 is fixed
      {
         Precision p(5e-3);
         p.importance = -1;
         p.divonne = p.cuhre = -1;
         SingleCase(5, 1, 1, p);
      }

      {
         Precision p(5e-3);
         p.importance = -1;
         p.divonne = -1e-2;
         p.cuhre = -1e-2;
         SingleCase(5, 3, 1, p);
      }

      {
         Precision p(5e-3);
         p.importance = -1;
         p.divonne = p.cuhre = -1;
         SingleCase(5, 3, 3, p);
      }

      // tough case: suave off by a factor of three,
      // divonne by 35, cuhre by 11
      {
         Precision p(5e-2);
         p.importance = -1;
         p.suave = -3.5;
         p.divonne = -35;
         p.cuhre = -11;
         SingleCase(15, 5, 3, p);
      }

      // test optimization
      {
         // two modes at 0 and 1
         IntegrationModel m = Factory(1,1,1);
         m.MCMCSetRandomSeed(1346);
         TEST_CHECK_NEARLY_EQUAL(m.FindMode(std::vector<double>(1, 0.1)).front(), 0, 1e-3);
         TEST_CHECK_NEARLY_EQUAL(m.FindModeSA(std::vector<double>(1, 0.1)).front(), 1, 1e-3);
         TEST_CHECK_NEARLY_EQUAL(m.FindModeSA(std::vector<double>(1, 0.9)).front(), 1, 1e-3);
         TEST_CHECK_NEARLY_EQUAL(m.FindModeMinuit(std::vector<double>(1, 0.1)).front(), 0, 1e-10);
         TEST_CHECK_NEARLY_EQUAL(m.FindModeMinuit(std::vector<double>(1, 0.9)).front(), 1, 1e-10);
         TEST_CHECK_NEARLY_EQUAL(m.FindModeMCMC().front(), 1, 1e-3);

         m.PrintSummary();
         m.PrintShortFitSummary();
      }
   }
} bcIntegrateTest;