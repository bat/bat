/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include "test.h"

#include <models/base/BCEfficiencyFitter.h>
#include <models/base/BCGraphFitter.h>
#include <models/base/BCHistogramFitter.h>

#include <TCanvas.h>

using namespace test;

class FitterTest :
    public TestCase
{
public:
    FitterTest() :
        TestCase("fitter_test"),
        fRandom(1234)
    {
    }

    template <class ConcreteFitter>
    void TestCopy(ConcreteFitter& m) const
    {
        // basic test: copy and assignment don't lead to segfaults
        ConcreteFitter copy(m);
        ConcreteFitter assigment(m);
        assigment = m;
    }

    void test_hist() const
    {
        TF1 f("f", "200/sqrt(2*pi)/[1] * exp(-0.5*((x-[0])/[1])^2)", -5, 5);
        f.SetParNames("mu", "sigma");
        f.SetParLimits(0, -2, 2);
        f.SetParLimits(1, 0, 5);

        // create data from standard Gaussian
        static const double mu = 0;
        static const double sigma = 1;
        TH1D hist("data", ";x;N", 20, -5, 5);
        for (int i = 0; i < 200; ++i)
            hist.Fill(fRandom.Gaus(mu, sigma));

        BCHistogramFitter hf(hist, f);
        hf.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        hf.SetPrecision(BCEngineMCMC::kQuick);
        hf.SetFlagIntegration(false);
        hf.Fit();

        TEST_CHECK_NEARLY_EQUAL(hf.GetBestFitParameters()[0], mu, 0.05);
        TEST_CHECK_NEARLY_EQUAL(hf.GetBestFitParameters()[1], sigma, 0.05);

        TEST_CHECK_NEARLY_EQUAL(hf.GetStatistics().mean[0], mu, 0.05);
        TEST_CHECK_NEARLY_EQUAL(hf.GetStatistics().mean[1], sigma, 0.05);
#if 0
        TCanvas c1("c1");
        hf.DrawFit("", true); // draw with a legend
        c1.Print("histfit.pdf");
#endif

        TestCopy(hf);
    }

    virtual void run() const
    {
        test_hist();
    }

private:
    mutable TRandom3 fRandom;
} fitter_test;
