/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>
#include <BAT/BCMath.h>
#include <BAT/BCEngineMCMC.h>
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
        BCMath::CacheFactorials(100);

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

class GammaTest :
    public TestCase
{
public:
    GammaTest() :
        TestCase("Gamma test")
    {
    }

    virtual void run() const
    {
        // check mean and variance
        {
            TRandom3 rng;
            rng.SetSeed(235);

            static const unsigned N = 10000;
            static const double a = 1.5;
            static const double b = 1.0 / 25.0;

            // abuse stats object for one parameter, no observables, dummy probability
            BCEngineMCMC::Statistics s(1, 0);
            BCEngineMCMC::ChainState cs(0);
            cs.log_probability = 0.;
            cs.parameters.assign(1, 0);
            for (unsigned i = 0; i < N; ++i) {
                cs.parameters.front() = BCMath::Random::Gamma(&rng, a, b);
                s.Update(cs);
            }
            // 1st moment: O(1/sqrt(N)), 2nd moment: order of magnitude larger
            TEST_CHECK_RELATIVE_ERROR(s.mean.front(), a * b, 1.0 / std::sqrt(N));
            TEST_CHECK_RELATIVE_ERROR(s.variance.front(), a * b * b, 10.0 / std::sqrt(N));
        }
    }
} gammaTest;
