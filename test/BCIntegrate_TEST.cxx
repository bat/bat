/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

#include <config.h>

#include "IntegrationModel.h"
#include "GaussModel.h"
#include "test.h"

#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

using namespace test;

class BCIntegrateTest :
    public TestCase
{
public:

    struct Precision {
        double monteCarlo, grid, vegas, suave, divonne, cuhre, cubaDefault;

        // Turn off checks by default
        Precision(double value = -1)
        {
            Set(value);
        }

        void Set(double value)
        {
            monteCarlo = grid = vegas = suave = divonne = cuhre = cubaDefault = value;
        }

    };

    BCIntegrateTest() :
        TestCase("BCIntegrate test")
    {
    }

    static IntegrationModel Factory(unsigned dim, unsigned modality, unsigned complexity)
    {
        IntegrationModel m;
        // set seed immediately before polynomial degrees are
        // populated through random number generator
        m.SetRandomSeed(1346);
        m.SetDimensionality(dim);
        m.SetModality(modality);
        m.SetComplexity(complexity);
        m.PopulatePolynomialDegrees();

        return m;
    }

#if HAVE_CUBA_H
    static void IntegrateCuba(IntegrationModel& m, BCIntegrate::BCCubaMethod method, const double& eps)
    {
        // aim for higher precision, as naively 1/3 of tests would lie outside interval
        // if eps were used
        m.SetRelativePrecision(eps / 3.);
        m.SetAbsolutePrecision(1e-20);
        m.SetCubaIntegrationMethod(method);
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(), m.Integral(), eps);
    }
#endif
    static void IntegrateBat(IntegrationModel& m, BCIntegrate::BCIntegrationMethod method, const double& eps)
    {
        // aim for higher precision, as naively 1/3 of tests would lie outside interval
        // if eps were used
        m.SetRelativePrecision(eps / 2.);
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(method), m.Integral(), eps);
    }

    static void IntegrateAllMethods(unsigned dim, unsigned modality, unsigned complexity, const Precision& p)
    {
        IntegrationModel m = Factory(dim, modality, complexity);

        std::cout << "Computing single integration test in " << dim
                  << " dimensions with modality = " << modality
                  << " and complexity = " << complexity
                  << " with target value " << m.Integral() << std::endl;

        if (p.monteCarlo > 0)
            IntegrateBat(m, BCIntegrate::kIntMonteCarlo, p.monteCarlo);
        if (p.grid > 0)
            IntegrateBat(m, BCIntegrate::kIntGrid, p.grid);

        // don't use Laplace: most test functions don't resemble a Gaussian
#if HAVE_CUBA_H
        if (p.vegas > 0)
            IntegrateCuba(m, BCIntegrate::kCubaVegas, p.vegas);
        if (p.suave > 0)
            IntegrateCuba(m, BCIntegrate::kCubaSuave, p.suave);
        if (p.divonne > 0)
            IntegrateCuba(m, BCIntegrate::kCubaDivonne, p.divonne);
        if (p.cuhre > 0)
            IntegrateCuba(m, BCIntegrate::kCubaCuhre, p.cuhre);
        if (p.cubaDefault > 0)
            IntegrateCuba(m, BCIntegrate::kCubaDefault, p.cubaDefault);
#endif
    }

    void Optimization() const
    {
        // two modes at 0 and 1
        IntegrationModel m = Factory(1, 1, 1);

        // can't be sure in which order the best value is found. If the
        // 2nd optimization doesn't improve on previous call, the old mode
        // is returned.
        m.SetFlagIgnorePrevOptimization(true);
        TEST_CHECK_NEARLY_EQUAL(m.FindMode(std::vector<double>(1, 0.1)).front(), 0, 1e-3);

        // minuit
        {
            static const double error = 1.19878e-01;

            // run optimization
            std::vector<double> mode = m.FindMode(BCIntegrate::kOptMinuit, std::vector<double>(1, 0.1));
            TMinuitMinimizer& minuit = m.GetMinuit();
            TEST_CHECK(minuit.Hesse());

            TEST_CHECK_NEARLY_EQUAL(mode.front(), 0, 1e-10);
            TEST_CHECK_NEARLY_EQUAL(m.GetBestFitParameterErrors().front(), error, 1e-5);
            TEST_CHECK_NEARLY_EQUAL(minuit.Errors()[0], error, 1e-5);

            /* mode at the left boundary, so covariance matrix is meaningless */
            double errhi, errlo;
            minuit.GetMinosError(0, errlo, errhi);

            /* but error and errhi only agree to 10%. I don't understand how `error` is computed */
            TEST_CHECK_RELATIVE_ERROR(errhi, error, 0.1);

            TEST_CHECK_NEARLY_EQUAL(m.FindMode(BCIntegrate::kOptMinuit, std::vector<double>(1, 0.85)).front(), 1, 1e-10);
            TEST_CHECK_NEARLY_EQUAL(m.GetBestFitParameterErrors().front(), error, 1e-5);
        }

        TEST_CHECK_NEARLY_EQUAL(m.FindMode(BCIntegrate::kOptSimAnn, std::vector<double>(1, 0.1)).front(), 1, 1e-3);
        // for different seeds, find either 0 or 1. Doesn't depend (much?) on starting point
        m.SetRandomSeed(6414);
        TEST_CHECK_NEARLY_EQUAL(m.FindMode(BCIntegrate::kOptSimAnn, std::vector<double>(1, 0.9)).front(), 1, 1e-3);

        // chain can find either mode
        m.SetNIterationsPreRunMin(5000);
        m.SetRandomSeed(1346);
        const double& mode = m.FindMode(BCIntegrate::kOptMetropolis).front();
        const double target = mode > 0.5 ? 1.0 : 0.0;
        TEST_CHECK_NEARLY_EQUAL(mode, target, 1e-3);

        m.PrintSummary();
        m.PrintShortFitSummary();
    }

    void FixedParameters(const unsigned ndim) const
    {
        GaussModel m("Fixed parameter example", ndim);
        m.SetRandomSeed(613);
        if (ndim >= 4)
            m.GetParameter(3).Fix(0.5);
        // set to boundary of parameter range
        // but the range should be small so the integral doesn't vanish
        if (ndim == 5) {
            BCParameter& p = m.GetParameter(4);
            p.SetLimits(-15, 3);
            p.Fix(3);
        }

        /* optimization */

        // perturb starting point
        std::vector<double> mode(ndim, -0.3);
        for (unsigned i = 0; i < ndim; ++i)
            mode.at(i) = mode.front() + i * 0.5351;

        mode = m.FindMode(BCIntegrate::kOptMinuit, mode);
        TMinuitMinimizer& minuit = m.GetMinuit();
        minuit.Hesse();

        std::vector<double> cov(ndim * ndim);
        minuit.GetCovMatrix(&cov[0]);

        std::vector<double> hessian(cov);
        minuit.GetHessianMatrix(&hessian[0]);

        double errlo, errhi;

        for (unsigned i = 0 ; i < ndim ; ++i) {
            const BCParameter& p = m.GetParameter(i);

            /* fixed parameter induces row/col. with zeros. */
            if (p.Fixed()) {
                TEST_CHECK_EQUAL(mode.at(i), p.GetFixedValue());
                TEST_CHECK_EQUAL(m.GetBestFitParameterErrors().at(i), 0.0);
                for (unsigned j = 0; j < ndim; ++j)
                    if (i != j) {
                        TEST_CHECK_EQUAL(cov.at(i * ndim + j), 0.0);
                        TEST_CHECK_EQUAL(cov.at(j * ndim + i), 0.0);
                    }
            }
            /* Expect result proportional to unit matrix for both
             * covariance and hessian, no correlation */
            else {
                TEST_CHECK_NEARLY_EQUAL(mode.at(i), m.mean(), 5e-4);

                const double eps = 2e-2;
                TEST_CHECK_RELATIVE_ERROR(m.GetBestFitParameterErrors().at(i), m.sigma(), eps);
                TEST_CHECK(minuit.GetMinosError(i, errlo, errhi));
                TEST_CHECK_RELATIVE_ERROR(-errlo, m.sigma(), eps);
                TEST_CHECK_RELATIVE_ERROR(errhi, m.sigma(), eps);
                TEST_CHECK_RELATIVE_ERROR(cov.at(i * ndim + i), m.sigma() * m.sigma(), eps);
                TEST_CHECK_RELATIVE_ERROR(hessian.at(i * ndim + i), 1.0 / (m.sigma() * m.sigma()), eps);

                for (unsigned j = 0; j < ndim; ++j) {
                    /* get NaN if comparing to fixed parameter */
                    if (i != j && !m.GetParameter(j).Fixed())
                        TEST_CHECK_NEARLY_EQUAL(minuit.Correlation(i, j), 0.0, 1e-3);
                }
            }
        }

        const double evidence = m.evidence();

        // Ask for higher precision than what we compare to. The error
        // is of a stochastic nature, and if we could interpret it as
        // a Gaussian standard deviation, then with 30% the true value
        // is outside the quoted interval.
        static const double eps = 5e-2;
        m.SetAbsolutePrecision(1e-12);

        // Laplace itself is naturally implemented on the log scale
        TEST_CHECK_RELATIVE_ERROR(m.IntegrateLaplace(), std::log(evidence), eps);

        // but the general interface is on the linear scale
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(BCIntegrate::kIntLaplace), evidence, eps);

        // sample mean needs huge number of evaluations
        m.SetNIterationsMax(3e6);
        m.SetRelativePrecision(eps / 3);
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(BCIntegrate::kIntMonteCarlo), evidence, eps);

#if HAVE_CUBA_H
        m.SetIntegrationMethod(BCIntegrate::kIntCuba);

        // cuba needs far less iterations
        m.SetNIterationsMax(1e5);
        m.SetRelativePrecision(eps / 5);

        m.SetCubaIntegrationMethod(BCIntegrate::kCubaVegas);
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(), evidence, eps);

        // suave systematically wrong by an order of magnitude even with huge number of evaluations.
        // I tried combinations of parameters with no success (cuba v4.2, Sep 25 2015)
#if 0
        m.SetCubaIntegrationMethod(BCIntegrate::kCubaSuave);
        BCCubaOptions::Suave o = m.GetCubaSuaveOptions();
        o.flatness = 5;
        o.nnew = 5000;
        o.nmin = 15;
        m.SetNIterationsMax(1e7);
        m.SetCubaOptions(o);
        TEST_CHECK_RELATIVE_ERROR(m.Integrate(), evidence, eps);
#endif
        if (ndim > 1) {
            m.SetCubaIntegrationMethod(BCIntegrate::kCubaDivonne);
            TEST_CHECK_RELATIVE_ERROR(m.Integrate(), evidence, eps);

            m.SetCubaIntegrationMethod(BCIntegrate::kCubaCuhre);
            TEST_CHECK_RELATIVE_ERROR(m.Integrate(), evidence, eps);
        }
#endif

    }

    void Integration() const
    {
        // top down: normalization should work automatically
        {
            IntegrationModel m = Factory(2, 1, 1);
            static const double eps = 1e-2;
            m.SetRelativePrecision(eps / 2);
            m.Normalize();
            TEST_CHECK_RELATIVE_ERROR(m.GetIntegral(), m.Integral(), eps * 1.3);
        }

        {
            Precision p;
            p.monteCarlo = 1e-2;
            p.grid = 1e-2;
            p.vegas = 2e-2;
            p.suave = 1e-2;
            p.divonne = -1; // fails in 1D
            p.cuhre = -1; // fails in 1D
            IntegrateAllMethods(1, 1, 1, p);
        }

        {
            Precision p(2e-2);
            p.divonne = p.cuhre = -2e-3;
            IntegrateAllMethods(2, 1, 1, p);
        }

        {
            Precision p(5e-3);
            p.grid = -1; // not available in 5D
            p.suave = -1;
            IntegrateAllMethods(5, 1, 1, p);
        }

        // tough case: low precision required
        {
            Precision p(0.1);
            p.grid = -1; // not available in 5D
            IntegrateAllMethods(5, 8, 4, p);
        }

        // run Laplace without running minuit before
        {
            GaussModel m("Gauss Laplace", 3);

            TEST_CHECK_RELATIVE_ERROR(m.Integrate(BCIntegrate::kIntLaplace), m.evidence(), 5e-8);
        }
    }

    void Slice() const
    {
        GaussModel m("slice", 1);

        // set bins of width 0.1
        // maximum should be at the center of one of the bins that have an edge at zero
        m.GetParameter(0).SetLimits(-3, 3);
        m.GetParameter(0).SetNbins(60);
        m.MarginalizeAll(BCIntegrate::kMargGrid);
        TEST_CHECK_RELATIVE_ERROR(0.05, std::abs(m.GetBestFitParameters()[0]), 1e-14);

        // mode finding should start from previous solution and converge to (0, 0)
        m.FindMode();
        TEST_CHECK_NEARLY_EQUAL(0, m.GetBestFitParameters()[0], 5e-5);
    }

    void Observables() const
    {
        static const unsigned npar = 3;
        GaussModel m("observables", npar);
        m.AddObservable("obs1", 0, 20, "#obs_1");
        TEST_CHECK_EQUAL(m.GetNObservables(), 1);

        {
            m.MarginalizeAll();

            /* make sure observables don't show up in best-fit parameters */
            TEST_CHECK_EQUAL(m.GetNParameters(), npar);
            TEST_CHECK_EQUAL(m.GetBestFitParameters().size(), npar);
        }

        {
            m.FindMode();

            TEST_CHECK_EQUAL(m.GetNParameters(), npar);
            TEST_CHECK_EQUAL(m.GetBestFitParameters().size(), npar);
            TEST_CHECK_EQUAL(m.GetBestFitParameterErrors().size(), npar);
        }
    }

    virtual void run() const
    {
        Observables();
        Integration();
        Optimization();
        FixedParameters(2);
        FixedParameters(5);
        Slice();
    }
} bcIntegrateTest;
