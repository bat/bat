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

#include <BAT/BCCauchyPrior.h>
#include <BAT/BCConstantPrior.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCSplitGaussianPrior.h>
#include <BAT/BCTF1LogPrior.h>
#include <BAT/BCTF1Prior.h>
#include <BAT/BCTH1Prior.h>

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
    void TestPriorDistribution(BCPrior * prior, double xmin, double xmax, double mode, double integral) const
    {
        // mode
        if (std::isfinite(mode))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMode(xmin,xmax), mode, 1.e-5 );
        else
            TEST_CHECK_EQUAL( prior->GetMode(xmin,xmax), mode );
        
        if (std::isfinite(integral))
            TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(xmin,xmax), integral, 1.e-5);
        else
            TEST_CHECK_EQUAL( prior->GetIntegral(xmin,xmax), integral);
    }

    // check the mode, integral, mean, and variance of a prior distribution
    // ignores variance check if negative
    void TestPriorDistribution(BCPrior * prior, double xmin, double xmax, double mode, double integral, double mean, double variance=-1) const
    {
        // check mode & integral
        TestPriorDistribution(prior,xmin,xmax,mode,integral);

        // mean
        if (std::isfinite(mean))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMean(xmin,xmax), mean, 1.e-5 );
        else
            TEST_CHECK_EQUAL( prior->GetMean(xmin,xmax), mean );
        
        // variance
        if (variance >= 0) {
            if (std::isfinite(variance))
                TEST_CHECK_NEARLY_EQUAL( prior->GetVariance(xmin,xmax), variance, 1.e-5 );
            else
                TEST_CHECK_EQUAL( prior->GetVariance(xmin,xmax), variance );
        }
    }


    // Test the implementation of the prior against the implementation
    // of the mode, integral, and moment (upto maxn'th moment)
    void TestPriorImplementation(BCPrior * prior, double xmin, double xmax, unsigned maxn) const
    {
        // BCPrior::Functions use ROOT routines to calculate from implementation of prior

        // mode
        TEST_CHECK_NEARLY_EQUAL( prior->GetMode(xmin,xmax),
                                 prior->BCPrior::GetMode(xmin,xmax),
                                 1.e-5 ); 
        
        // integral
        TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(xmin,xmax),
                                 prior->BCPrior::GetIntegral(xmin,xmax),
                                 1.e-5 ); 
        
        // raw moments
        for (unsigned n=1; n<=maxn; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetRawMoment(n,xmin,xmax),
                                     prior->BCPrior::GetRawMoment(n,xmin,xmax),
                                     1.e-5 );

        // central moments
        for (unsigned n=1; n<=maxn; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetCentralMoment(n,xmin,xmax),
                                     prior->BCPrior::GetCentralMoment(n,xmin,xmax),
                                     1.e-5 );
    }

    virtual void run() const
    {

        double neg_inf = -std::numeric_limits<double>::infinity();
        double pos_inf = +std::numeric_limits<double>::infinity();

        BCPrior * prior = NULL;

        // Cauchy Prior
        std::cout << "Testing BCCauchyPrior ... " << std::flush;
        prior = new BCCauchyPrior(1.5,3);
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
        std::cout << "Testing BCConstantPrior ... " << std::flush;
        prior = new BCConstantPrior();
        // check moments
        for (unsigned n=1; n<=4; ++n)
            TEST_CHECK_NEARLY_EQUAL( prior->GetRawMoment(n,-10,10),
                                     prior->BCPrior::GetRawMoment(n,-10,10),
                                     1.e-5 );
        TEST_CHECK_EQUAL( prior->GetLogPrior(0.4) , 0 );
        //                            xmin,    xmax,    mode,    I, mean,    var
        TestPriorDistribution( prior, neg_inf, pos_inf, 0,       1, 0,       pos_inf); // [ -inf , +inf ]
        TestPriorDistribution( prior, 0,       1,       0.5,     1, 0.5,     1./12);   // [  fin ,  fin ]
        TestPriorDistribution( prior, 10,      pos_inf, pos_inf, 1, pos_inf, pos_inf); // [  fin , +inf ]
        TestPriorDistribution( prior, neg_inf, 10,      neg_inf, 1, neg_inf, pos_inf); // [ -inf ,  fin ]
        delete prior;
        std::cout << "PASS" << std::endl;

        // Gaussian Prior
        std::cout << "Testing BCGaussianPrior ... " << std::flush;
        prior = new BCGaussianPrior(1.5,3);
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
        std::cout << "Testing BCSplitGaussianPrior ... " << std::flush;
        prior = new BCSplitGaussianPrior(1.5,3,5);
        TestPriorImplementation( prior, -10, 10, 1);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.30523 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(6.5),  -2.80523 , 1.e-5); // at mean+sigma_above
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.30523 , 1.e-5); // at mean-2*sigma_below
        //                            xmin,    xmax,    mode, I,       mean, var
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,  1,       3.09577); // [ -inf , +inf ]
        TestPriorDistribution( prior, -4.5,    6.5,     1.5,  0.78462, 1.76119); // [  fin, fin  ]
        TestPriorDistribution( prior, 6.5,     pos_inf, 6.5,  0.19832, 9.12568); // [  fin, +inf ] above mean
        delete prior;
        std::cout << "PASS" << std::endl;
        
        
        // TF1LogPrior

        // TF1 Prior

        // TH1Prior

    }

} bcaux_Test;
