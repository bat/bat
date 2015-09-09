/*
 * Copyright (C) 2009, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#include "include/ReleaseTestSuite.h"
#include <include/PerfTestVarPar.h>
#include <include/PerfTest1DFunction.h>
#include <include/PerfTest2DFunction.h>

#include <BAT/BCParameter.h>
#include <iostream>

//______________________________________________________________________________
ReleaseTestSuite::ReleaseTestSuite(bool multivariate, double dof):
    TestSuite(multivariate, dof)
{
}

//______________________________________________________________________________
int ReleaseTestSuite::PrepareTests()
{
    /* 1D */

    // these functions are covered by binomial
#if 1
    // 1D slope
    TF1* testfunc_1d_slope = new TF1("Slope", "x", 0., 10.);
    PerfTest1DFunction* perftest_1d_slope = new PerfTest1DFunction("1d_slope", testfunc_1d_slope);
    AddTest(perftest_1d_slope);

    // 1D squared
    TF1* testfunc_1d_squared = new TF1("Squared", "400.-x*x", -20., 20.);
    PerfTest1DFunction* perftest_1d_squared = new PerfTest1DFunction("1d_squared", testfunc_1d_squared);
    AddTest(perftest_1d_squared);
#endif

#if 1
    // 1D Gaussian
    TF1* testfunc_1d_gaus = new TF1("Gaus", "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2/[1]/[1])", -25., 25.);
    testfunc_1d_gaus->FixParameter(0, 0.0);
    testfunc_1d_gaus->FixParameter(1, 5.0);
    PerfTest1DFunction* perftest_1d_gaus = new PerfTest1DFunction("1d_gaus", testfunc_1d_gaus);
    AddTest(perftest_1d_gaus);

    // 1D Poissons
    for (int i = 0; i <= 2; i++) {
        double xmax = 15;
        if (i > 3)
            xmax = 10.0 * sqrt(double(i));
        TF1* testfunc_1d_poisson = new TF1("Poisson", "TMath::PoissonI([0], x)", 0., xmax);
        testfunc_1d_poisson->FixParameter(0, double(i));
        PerfTest1DFunction* perftest = new PerfTest1DFunction(Form("1d_poisson_%i", i), testfunc_1d_poisson);
        AddTest(perftest);
    }

    // 1D Binomials
    for (int N = 1; N < 3; N++) {
        for (int k = 0; k <= N; ++k) {
            TF1* testfunc = new TF1("Binomial", "([0]+1)*TMath::Binomial([0], [1]) * TMath::Power(x, [1]) * TMath::Power(1-x, [0]-[1])", 0., 1.);
            testfunc->FixParameter(0, N);
            testfunc->FixParameter(1, k);
            PerfTest1DFunction* perftest = new PerfTest1DFunction(Form("1d_binomial_%i_%i", k, N), testfunc);
            AddTest(perftest);
        }
    }

    // 1D exponential
    TF1* testfunc_1d_exponential = new TF1("Exponential", "1/[0]*exp(-x/[0])", 0., 100.);
    testfunc_1d_exponential->FixParameter(0, 5);
    PerfTest1DFunction* perftest_1d_exponential = new PerfTest1DFunction("1d_exponential", testfunc_1d_exponential);
    AddTest(perftest_1d_exponential);

    // 1D Cauchy
    TF1* testfunc_1d_cauchy = new TF1("Cauchy", "[1] / (3.14159 * ( (x-[0])**2 +[1]**2))", -50., 50.);
    testfunc_1d_cauchy->FixParameter(0, 0.);
    testfunc_1d_cauchy->FixParameter(1, 5.);
    PerfTest1DFunction* perftest_1d_cauchy = new PerfTest1DFunction("1d_cauchy", testfunc_1d_cauchy);
    AddTest(perftest_1d_cauchy);

    // 1D Lognormal
    TF1* testfunc_1d_lognormal = new TF1("Lognormal", "1./sqrt(2*TMath::Pi()*[1])*1/x*exp(-(log(x)-[0])*(log(x)-[0])/2/[1]/[1])", 0., 10.);
    testfunc_1d_lognormal->FixParameter(0, 0.);
    testfunc_1d_lognormal->FixParameter(1, 1.);
    PerfTest1DFunction* perftest_1d_lognormal = new PerfTest1DFunction("1d_lognormal", testfunc_1d_lognormal);
    AddTest(perftest_1d_lognormal);

    // 1D x^4 sin^2(x)
    TF1* testfunc_1d_sin2 = new TF1("x4Sin2", "x*x*x*x*sin(x)*sin(x)", 2., 25.);
    PerfTest1DFunction* perftest_1d_sin2 = new PerfTest1DFunction("1d_sin2", testfunc_1d_sin2);
    AddTest(perftest_1d_sin2);

    // 1D 2 Gaussians
    TF1* testfunc_1d_2gaus = new TF1("2gaus1d", "gaus + gaus(3)", -25., 50.);
    testfunc_1d_2gaus->FixParameter(0,  1.0);
    testfunc_1d_2gaus->FixParameter(1, -10.0);
    testfunc_1d_2gaus->FixParameter(2,  2.0);
    testfunc_1d_2gaus->FixParameter(3,  2.0);
    testfunc_1d_2gaus->FixParameter(4, 30.0);
    testfunc_1d_2gaus->FixParameter(5,  1.0);
    PerfTest1DFunction* perftest_1d_2gaus = new PerfTest1DFunction("1d_2gaus", testfunc_1d_2gaus);
    perftest_1d_2gaus->GetParameter("x").SetNbins(200);
    perftest_1d_2gaus->GetSubtest("mode")->SetStatusOff(true);
    AddTest(perftest_1d_2gaus);
#endif

    /* 2D */
#if 1

    // 2D flat
//   TF2* testfunc_2d_flat = new TF2("Flat", "1", -5., 5., -5., 5.);
// the above definition of a constant 2d function doesn't work in Root (not sure why)
// so we use a workaround
// using TF2("Flat", "1", 1.) would also work but would not have the range defined correctly
    TF2* testfunc_2d_flat = new TF2("Flat", "y*0. + 1", -5., 5., -5., 5.);
    PerfTest2DFunction*   perftest_2d_flat = new PerfTest2DFunction("2d_flat", testfunc_2d_flat);
    AddTest(perftest_2d_flat);

    // 2D Gaussian product Gaussian
    TF2* testfunc_2d_gaus = new TF2("Gaus", "xygaus", -3., 3., -5., 7.);
    // const,meanx,sigmax,meany,sigmay
    testfunc_2d_gaus->SetParameters(1, 0, 1, 1, 2);
    PerfTest2DFunction*   perftest_2d_gaus = new PerfTest2DFunction("2d_gaus", testfunc_2d_gaus);
    AddTest(perftest_2d_gaus);

    // 2D 2 Gaussians
    TF2* testfunc_2d_2gaus = new TF2("2twoGaus2d",
                                     "[0] * ( [1]*exp(-0.5*((x-[2])/[3])**2)*exp(-0.5*((y-[4])/[5])**2) + [6]*exp(-0.5*((x-[7])/[8])**2)*exp(-0.5*((y-[9])/[10])**2))",
                                     -20., 20., -20., 20);
    testfunc_2d_2gaus->SetParameters(1.,   10., 0., 1.0,  5., 1.0,    10., 5., 1.0,  10., 1.0);
    PerfTest2DFunction*   perftest_2d_2gaus = new PerfTest2DFunction("2d_2gaus", testfunc_2d_2gaus);
    AddTest(perftest_2d_2gaus);
#endif

    /* variable parameters */
#if 1
    std::vector<double> values_lag;
    for (int i = 1; i <= 10; ++i)
        values_lag.push_back(i);

    // 1D Gauss with varying lag
    TF1* testfunc_1d_gaus_lag = new TF1("1dGaus_lag", "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2/[1]/[1])", -25., 25.);
    testfunc_1d_gaus_lag->FixParameter(0, 0.0);
    testfunc_1d_gaus_lag->FixParameter(1, 5.0);
    PerfTest1DFunction*   perftest_1d_gaus_lag = new PerfTest1DFunction("1d_gaus_lag", testfunc_1d_gaus_lag);
    PerfTestVarPar* varpar_gaus_lag = new PerfTestVarPar("1d_gaus_lag", perftest_1d_gaus_lag);
    varpar_gaus_lag->AddVarPar(values_lag, "lag");
    AddTest(varpar_gaus_lag);

    // 2D Gaussian with varying lag
    TF2* testfunc_2d_gaus_lag = new TF2("2dGaus_lag", "xygaus", -3., 3., -5., 7.);
    testfunc_2d_gaus_lag->SetParameters(1, 0, 1, 1, 2);
    PerfTest2DFunction*   perftest_2d_gaus_lag = new PerfTest2DFunction("2d_gaus_lag", testfunc_2d_gaus_lag);
    PerfTestVarPar* varpar_2dgaus_lag = new PerfTestVarPar("2d_gaus_lag", perftest_2d_gaus_lag);
    varpar_2dgaus_lag->AddVarPar(values_lag, "lag");
    AddTest(varpar_2dgaus_lag);

    std::vector<double> values_iter;
    for (int i = 1; i <= 10; ++i)
        values_iter.push_back(i * 100000);

    // 1D Gauss with varying number of iterations
    TF1* testfunc_1d_gaus_iter = new TF1("Gaus", "1.0/sqrt(2.0*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2/[1]/[1])", -25., 25.);
    testfunc_1d_gaus_iter->FixParameter(0, 0.0);
    testfunc_1d_gaus_iter->FixParameter(1, 5.0);
    PerfTest1DFunction*   perftest_1d_gaus_iter = new PerfTest1DFunction("1d_gaus_iter", testfunc_1d_gaus_iter);
    PerfTestVarPar* varpar_gaus_iter = new PerfTestVarPar("1d_gaus_iter", perftest_1d_gaus_iter);
    varpar_gaus_iter->AddVarPar(values_iter, "iterations");
    AddTest(varpar_gaus_iter);

    // 2D Gaussian with varying iter
    TF2* testfunc_2d_gaus_iter = new TF2("2dGaus_iter", "xygaus", -3., 3., -5., 7.);
    testfunc_2d_gaus_iter->SetParameters(1, 0, 1, 1, 2);
    PerfTest2DFunction*   perftest_2d_gaus_iter = new PerfTest2DFunction("2d_gaus_iter", testfunc_2d_gaus_iter);
    PerfTestVarPar* varpar_2dgaus_iter = new PerfTestVarPar("2d_gaus_iter", perftest_2d_gaus_iter);
    varpar_2dgaus_iter->AddVarPar(values_iter, "iteration");
    AddTest(varpar_2dgaus_iter);
#endif

    // proposal for all MCMC  tests
    for (unsigned i = 0; i < GetNTests(); ++i)
        GetTest(i)->SetProposal(fMultivariate, fDof);

    // no error
    return 1;
}

//______________________________________________________________________________
void ReleaseTestSuite::WebpageSetup()
{
    IncludeHtmlHeader(false);
    IncludeHtmlFooter(false);
    SetLinkPrefix("<?php echo $linkPrefix; ?>");
    SetFileLinkPrefix("<?php echo $fileLinkPrefix; ?>");
    SetHtmlFileExtension(".php");
}

//______________________________________________________________________________
