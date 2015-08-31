/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>

#include <limits>
#include <vector>

#include <BAT/BCPrior.h>

#include <BAT/BCPriorCauchy.h>
#include <BAT/BCPriorConstant.h>
#include <BAT/BCPriorGaussian.h>
#include <BAT/BCPriorLogTF1.h>
#include <BAT/BCPriorSplitGaussian.h>
#include <BAT/BCPriorTF1.h>
#include <BAT/BCPriorTH1.h>
#include <TF1.h>
#include <TH1D.h>

using namespace test;

class BCPriorTest :
    public TestCase
{
public:
    BCPriorTest() :
        TestCase("BCPrior test")
    {
    }

    // check the mode and integral of a prior distribution
    void TestPriorDistribution(BCPrior* prior, double xmin, double xmax, double mode, double integral, double eps = 1.e-5) const
    {
        // mode
        if (std::isfinite(mode))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMode(xmin, xmax), mode, eps );
        else
            TEST_CHECK_EQUAL( prior->GetMode(xmin, xmax), mode );

        if (std::isfinite(integral))
            TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(xmin, xmax), integral, eps);
        else
            TEST_CHECK_EQUAL( prior->GetIntegral(xmin, xmax), integral);
    }

    // check the mode, integral, mean, and variance of a prior distribution
    // ignores variance check if negative
    void TestPriorDistribution(BCPrior* prior, double xmin, double xmax, double mode, double integral, double mean, double variance, double eps = 1.e-5) const
    {
        // check mode & integral
        TestPriorDistribution(prior, xmin, xmax, mode, integral);

        // mean
        if (std::isfinite(mean))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMean(xmin, xmax), mean, eps );
        else
            TEST_CHECK_EQUAL( prior->GetMean(xmin, xmax), mean );

        // variance
        if (std::isfinite(variance))
            TEST_CHECK_NEARLY_EQUAL( prior->GetVariance(xmin, xmax), variance, eps );
        else
            TEST_CHECK_EQUAL( prior->GetVariance(xmin, xmax), variance );
    }


    // Test the implementation of the prior against the implementation
    // of the mode, integral, and moment (upto maxn'th moment)
    void TestPriorImplementation(BCPrior* prior, double xmin, double xmax, unsigned maxn = 2, double eps = 1.e-5) const
    {
        // BCPrior::Functions use ROOT routines to calculate from implementation of prior

        // mode
        TEST_CHECK_NEARLY_EQUAL( prior->GetMode(xmin, xmax),
                                 prior->BCPrior::GetMode(xmin, xmax),
                                 eps );

        // integral
        TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(xmin, xmax),
                                 prior->BCPrior::GetIntegral(xmin, xmax),
                                 eps );

        // raw moments
        for (unsigned n = 1; n <= maxn; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetRawMoment(n, xmin, xmax),
                                     prior->BCPrior::GetRawMoment(n, xmin, xmax),
                                     eps );

        // central moments
        for (unsigned n = 1; n <= maxn; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetCentralMoment(n, xmin, xmax),
                                     prior->BCPrior::GetCentralMoment(n, xmin, xmax),
                                     eps );
    }

    virtual void run() const
    {

        double neg_inf = -std::numeric_limits<double>::infinity();
        double pos_inf = +std::numeric_limits<double>::infinity();

        BCPrior* prior = NULL;

        // Cauchy Prior
        std::cout << "Testing BCPriorCauchy ... " << std::flush;
        prior = new BCPriorCauchy(1.5, 3);
        TestPriorImplementation( prior, -10, 10, 1 );
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.24334 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -2.93649 , 1.e-5); // at mean+scale
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -3.85278 , 1.e-5); // at mean-2*scale
        //                            xmin,    xmax,    mode, integral, mean,    var
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,  1,        1.5,     pos_inf);  // [ -inf , +inf ]
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,  0.5,      1.5,     2.45916);  // [  fin , fin  ]
        TestPriorDistribution( prior, 3,       pos_inf, 3,    0.35241,  pos_inf, pos_inf);  // [  fin , +inf ]
        TestPriorDistribution( prior, neg_inf, 3,       1.5,  0.64758,  neg_inf, pos_inf);  // [ -inf ,  fin ]
        delete prior;
        std::cout << "PASS" << std::endl;

        // Constant Prior
        std::cout << "Testing BCPriorConstant ... " << std::flush;
        prior = new BCPriorConstant();
        // check moments
        for (unsigned n = 1; n <= 4; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetRawMoment(n, -10, 10),
                                     prior->BCPrior::GetRawMoment(n, -10, 10),
                                     1.e-5 );
        TEST_CHECK_EQUAL( prior->GetLogPrior(0.4) , 0 );
        //                            xmin,    xmax,    mode,    I, mean,    var
        TestPriorDistribution( prior, neg_inf, pos_inf, 0,       1, 0,       pos_inf); // [ -inf , +inf ]
        TestPriorDistribution( prior, 0,       1,       0.5,     1, 0.5,     1. / 12); // [  fin ,  fin ]
        TestPriorDistribution( prior, 10,      pos_inf, pos_inf, 1, pos_inf, pos_inf); // [  fin , +inf ]
        TestPriorDistribution( prior, neg_inf, 10,      neg_inf, 1, neg_inf, pos_inf); // [ -inf ,  fin ]
        delete prior;
        std::cout << "PASS" << std::endl;

        // Gaussian Prior
        std::cout << "Testing BCPriorGaussian ... " << std::flush;
        prior = new BCPriorGaussian(1.5, 3);
        TestPriorImplementation( prior, -10, 10, 1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.01755 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -2.51755 , 1.e-5); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.01755 , 1.e-5); // at mean-2*sigma
        //                            xmin,    xmax,    mode, integral, mean,    var,
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,  1,        1.5,     9);        // [ -inf , +inf ]
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,  0.68269,  1.5,     2.62013);  // [  fin , fin  ]
        TestPriorDistribution( prior, 4.5,     pos_inf, 4.5,  0.15866,  6.07541, 1.79188);  // [  fin , +inf ] above mean
        delete prior;
        std::cout << "PASS" << std::endl;

        // Split Gaussian Prior
        std::cout << "Testing BCPriorSplitGaussian ... " << std::flush;
        prior = new BCPriorSplitGaussian(1.5, 3, 5);
        TestPriorImplementation( prior, -10, 10, 1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.30523 , 1.e-5); // at mode
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(6.5),  -2.80523 , 1.e-5); // at mode+sigma_above
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.30523 , 1.e-5); // at mode-2*sigma_below
        //                            xmin,    xmax,    mode, I,       mean,    var
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,  1,       3.09577, 16.45352); // [ -inf , +inf ]
        TestPriorDistribution( prior, -4.5,    6.5,     1.5,  0.78462, 1.76119,  7.06646); // [  fin, fin  ]
        TestPriorDistribution( prior, 6.5,     pos_inf, 6.5,  0.19832, 9.12568,  4.97744); // [  fin, +inf ] above mean
        delete prior;
        std::cout << "PASS" << std::endl;

        // TF1LogPrior
        std::cout << "Testing BCPriorLogTF1 ... " << std::flush;
        // TF1 of ln(Normal)
        prior = new BCPriorLogTF1("-0.5*((x-[0])/[1])^2 - log([1]) - 0.5*log(2*pi)", -10, 10);
        dynamic_cast<BCPriorLogTF1*>(prior)->GetLogFunction().SetParameters(1.5, 3);
        TestPriorImplementation( prior, -10, 10, 2);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.01755 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -2.51755 , 1.e-5); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.01755 , 1.e-5); // at mean-2*sigma
        //                            xmin,    xmax,    mode, integral, mean,    var,
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,  0.68269,  1.5,     2.62013);  // [  fin , fin  ]
        delete prior;
        std::cout << "PASS" << std::endl;


        // TF1 Prior
        std::cout << "Testing BCPriorTF1 ... " << std::flush;
        // TF1 of Normal
        prior = new BCPriorTF1("exp(-0.5*((x-[0])/[1])^2)/sqrt(2*pi)/[1]", -10, 10);
        prior->GetFunction().SetParameters(1.5, 3);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.01755 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -2.51755 , 1.e-5); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.01755 , 1.e-5); // at mean-2*sigma
        //                            xmin,    xmax,    mode, integral, mean,    var,
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,  0.68269,  1.5,     2.62013);  // [  fin , fin  ]
        delete prior;
        std::cout << "PASS" << std::endl;

        // TH1Prior
        std::cout << "Testing BCPriorTH1 with interpolation ... " << std::flush;
        TH1D h1_prior("h1_prior", "", 100, -5, 5);
        h1_prior.FillRandom("gaus", 1000000);
        prior = new BCPriorTH1(h1_prior, true);
        TestPriorImplementation( prior, -10, 10, 2, 1.e-2);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(0),  -0.91894 , 5.e-2); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1),  -1.41894 , 5.e-2); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-2), -2.91894 , 5.e-2); // at mean-2*sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetMode(-5, 5), 0, 1.e-1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(-5, 5), 1, 1.e-4);
        TEST_CHECK_NEARLY_EQUAL( prior->GetMean(-5, 5), 0, 1.e-2);
        TEST_CHECK_NEARLY_EQUAL( prior->GetVariance(-5, 5), 1, 1.e-2);
        delete prior;
        std::cout << "PASS" << std::endl;

        std::cout << "Testing BCPriorTH1 without interpolation ... " << std::flush;
        prior = new BCPriorTH1(h1_prior, false);
        TestPriorImplementation( prior, -10, 10, 2, 5.e-2);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(0),  -0.91894 , 0.15); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1),  -1.41894 , 0.15); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-2), -2.91894 , 0.15); // at mean-2*sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetMode(-5, 5), 0, 0.1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(-5, 5), 1, 1.e-4 );
        TEST_CHECK_NEARLY_EQUAL( prior->GetMean(-5, 5), 0, 0.1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetVariance(-5, 5), 1, 0.1);
        delete prior;
        std::cout << "PASS" << std::endl;

    }

} bcaux_Test;
