/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <test.h>

#include <vector>
#include <limits>

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

    // struct for testing evaluation of a prior
    struct EvaluationTestVals {
        EvaluationTestVals(double X, double L) : x(X), log_val(L) {}
        double x;               ///< point to evaluate at
        double log_val;         ///< true value of log(prior) at point
    };

    // struct for testing a prior's distribution
    struct DistributionTestVals {
        DistributionTestVals(double x0, double x1, double mn, double md, double var, double I)
            : xmin(x0), xmax(x1), mean(mn), mode(md), variance(var), integral(I) {}
        double xmin;            ///< lower limit of range to evaluate over
        double xmax;            ///< upper limit of range to evaluate over
        double mean;            ///< true value of mean in range
        double mode;            ///< true value of mode in range
        double variance;        ///< true value of variance in range
        double integral;        ///< true value of integral in range
    };

    struct TestVals {
        std::vector<EvaluationTestVals> evaluation_tests;
        std::vector<DistributionTestVals> distribution_tests;
    };

    // check a prior
    void TestPriorDistribution(BCPrior * prior, double xmin, double xmax, double mean, double mode, double variance, double integral) const {
        // mean
        if (std::isfinite(mean))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMean(xmin,xmax), mean, 1.e-5 );
        else
            TEST_CHECK_EQUAL( prior->GetMean(xmin,xmax), mean );
        
        // mode
        if (std::isfinite(mode))
            TEST_CHECK_NEARLY_EQUAL( prior->GetMode(xmin,xmax), mode, 1.e-5 );
        else
            TEST_CHECK_EQUAL( prior->GetMode(xmin,xmax), mode );
        
        // variance
        if (std::isfinite(variance))
            TEST_CHECK_NEARLY_EQUAL( prior->GetVariance(xmin,xmax), variance, 1.e-5 );
        else
            TEST_CHECK_EQUAL( prior->GetVariance(xmin,xmax), variance );
        
        if (std::isfinite(integral))
            TEST_CHECK_NEARLY_EQUAL( prior->GetIntegral(xmin,xmax), integral, 1.e-5);
        else
            TEST_CHECK_EQUAL( prior->GetIntegral(xmin,xmax), integral);
    }

    virtual void run() const
    {

        double neg_inf = -std::numeric_limits<double>::infinity();
        double pos_inf = +std::numeric_limits<double>::infinity();

        BCPrior * prior = NULL;

        // Cauchy Prior
        prior = new BCCauchyPrior(1.5,3);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.24334 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -2.93649 , 1.e-5); // at mean+scale
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -3.85278 , 1.e-5); // at mean-2*scale
        //                            xmin,    xmax,    mean,    mode,    var,     integral
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,     1.5,     pos_inf, 1);       // [ -inf , +inf ]
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,     1.5,     2.45916, 0.5);     // [  fin , fin  ]
        TestPriorDistribution( prior, 3,       pos_inf, pos_inf, 3,       pos_inf, 0.35241); // [  fin , +inf ]
        TestPriorDistribution( prior, neg_inf, 3,       neg_inf, 1.5,     pos_inf, 0.64758); // [ -inf ,  fin ]
        delete prior;

        // Constant Prior
        prior = new BCConstantPrior();
        TEST_CHECK_EQUAL( prior->GetLogPrior(0.4) , 0 );
        //                            xmin,    xmax,    mean,    mode,    var,     integral
        TestPriorDistribution( prior, neg_inf, pos_inf, 0,       0,       pos_inf, 1); // [ -inf , +inf ]
        TestPriorDistribution( prior, 0,       1,       0.5,     0.5,     1./12.,  1); // [  fin ,  fin ]
        TestPriorDistribution( prior, 10,      pos_inf, pos_inf, pos_inf, pos_inf, 1); // [  fin , +inf ]
        TestPriorDistribution( prior, neg_inf, 10,      neg_inf, neg_inf, pos_inf, 1); // [ -inf ,  fin ]
        delete prior;

        // Gaussian Prior
        prior = new BCGaussianPrior(1.5,3);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.93649 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(4.5),  -3.43649 , 1.e-5); // at mean+sigma
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.93649 , 1.e-5); // at mean-2*sigma
        //                            xmin,    xmax,    mean,    mode,    var,     integral
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,     1.5,     9,       1);       // [ -inf , +inf ]
        TestPriorDistribution( prior, -1.5,    4.5,     1.5,     1.5,     2.62013, 0.68269); // [  fin , fin  ]

        // Split Gaussian Prior
        prior = new BCSplitGaussianPrior(1.5,3,5);
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(1.5),  -2.53102 , 1.e-5); // at mean
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(6.5),  -3.03102 , 1.e-5); // at mean+sigma_above
        TEST_CHECK_NEARLY_EQUAL( prior->GetLogPrior(-4.5), -4.53102 , 1.e-5); // at mean-2*sigma_below
        //                            xmin,    xmax,    mean,    mode,    var,     integral
        TestPriorDistribution( prior, neg_inf, pos_inf, 1.5,     1.5,     9,       1);       // [ -inf , +inf ]
        // TestPriorDistribution( prior, -1.5,    4.5,     1.5,     1.5,     2.62013, 0.68269); // [  fin , fin  ]
        
        
        // TF1LogPrior

        // TF1 Prior

        // TH1Prior

    }

} bcaux_Test;
