/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>
#include <BAT/BCMath.h>
#include <cmath>

using namespace test;
using namespace BCMath;

class BCPValueTest :
public TestCase
{
public:
   BCPValueTest() :
      TestCase("BCPValue test")
   {
   }

   virtual void run() const
   {
      static const double eps = 1e-13;
      // CorrectPValue
      {
         // values taken from histogram fitter example
         TEST_CHECK_RELATIVE_ERROR(CorrectPValue(0.17166, 4, 20), 0.05653668738694474, eps);
         TEST_CHECK_RELATIVE_ERROR(CorrectPValue(0.67686, 4, 20), 0.4099294170252848, eps);
         TEST_CHECK_NEARLY_EQUAL(CorrectPValue(0, 4, 20), 0, eps);
         TEST_CHECK_NEARLY_EQUAL(CorrectPValue(1e-38, 4, 20), 0, eps);
         TEST_CHECK_NEARLY_EQUAL(CorrectPValue(1, 4, 20), 1, eps);

         TEST_CHECK_THROWS(std::domain_error, CorrectPValue(-0.5, 4, 20));
         TEST_CHECK_THROWS(std::domain_error, CorrectPValue(-1.5, 4, 20));
         TEST_CHECK_THROWS(std::domain_error, CorrectPValue(0.5, 20, 4));
      }

      // FastPValue
      {
         size_t nbins = 10;

         /* most likely case, cannot do better */
         std::vector<unsigned> observed(nbins, 5);
         std::vector<double> expected(nbins, 5);

         TEST_CHECK_RELATIVE_ERROR(FastPValue(observed, expected, 123), 1, eps);

         /* bad fit */
         expected = std::vector<double>(nbins, 0.1);

         TEST_CHECK_NEARLY_EQUAL(FastPValue(observed, expected, 1e5, 123), 0, eps);

         /* simple one bin */
         observed = std::vector<unsigned> (1, 1);
         expected = std::vector<double>(1, 1.3);

         TEST_CHECK_RELATIVE_ERROR(FastPValue(observed, expected, 1e6, 1235), 1, eps);

         /* p = 1 - P(1|1.3) */
         observed = std::vector<unsigned> (1, 0);
         TEST_CHECK_RELATIVE_ERROR(FastPValue(observed, expected, 1e6, 123), 1 - 0.35429133094421639, 5e-3);
      }
   }
} bcPValueTest;

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
         TEST_CHECK_RELATIVE_ERROR(Rvalue(chain_means, chain_variances, points, relaxed), 1.0011584199407115, eps);
         TEST_CHECK_RELATIVE_ERROR(Rvalue(chain_means, chain_variances, points, strict), 1.0176292831481546, eps);

         // for more points visited, R-value increases
         points *= 3;

         TEST_CHECK_RELATIVE_ERROR(Rvalue(chain_means, chain_variances, points, relaxed), 1.0018240939164496, eps);
         TEST_CHECK_RELATIVE_ERROR(Rvalue(chain_means, chain_variances, points, strict), 1.0183054631320092,eps);
      }
   }
} rvalue_test;
