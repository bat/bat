/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMath.h"

#include "BCLog.h"

#include <Math/PdfFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>
#include <TMath.h>
#include <TRandom3.h>

#include <limits>
#include <math.h>

// ---------------------------------------------------------
double BCMath::LogGaus(double x, double mean, double sigma, bool norm)
{
    if (sigma == 0.)							// delta function !
        return (fabs(x - mean) < std::numeric_limits<double>::epsilon()) ? std::numeric_limits<double>::infinity() : 0;

    if (sigma < 0.)						// if sigma is negative use absolute value
        sigma *= -1.;

    double arg = (x - mean) / sigma;
    double result = -.5 * arg * arg;

    // check if we should add the normalization constant
    if (!norm)
        return result;

    // return result with subtraction of log of normalization constant
    return result - 0.5 * log(2 * M_PI) - log(sigma);
}

// ---------------------------------------------------------
double BCMath::LogSplitGaus(double x, double mode, double sigma_below, double sigma_above, bool norm)
{
    double norm_const = (norm) ? 0.5 * log(2 / M_PI) - log(sigma_below + sigma_above) : 0;
    if (x > mode)
        return LogGaus(x, mode, sigma_above, false) + norm_const;
    return LogGaus(x, mode, sigma_below, false) + norm_const;
}

// ---------------------------------------------------------
double BCMath::LogPoisson(double x, double lambda)
{
    if (x < 0)
        return -std::numeric_limits<double>::infinity();

    if (lambda == 0) {
        BCLog::OutWarning("BCMath::LogPoisson : expectation value (lambda) cannot be zero.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (lambda < 0) {
        BCLog::OutWarning("BCMath::LogPoisson : expectation value (lambda) cannot be negative.");
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (lambda > 899)
        return LogGaus(x, lambda, sqrt(lambda), true);

    if (x == 0.)
        return -lambda;

    return x * log(lambda) - lambda - ApproxLogFact(x);
}

// ---------------------------------------------------------
double BCMath::LogApproxBinomial(unsigned n, unsigned k, double p)
{
    if (k > n)										// impose k < n requirement
        return std::numeric_limits<double>::quiet_NaN();

    // p in [0,1]
    if (p < 0 or p > 1)
        return std::numeric_limits<double>::quiet_NaN();

    if (p == 0)										// no chance of success
        return (k == 0) ? 0 : -std::numeric_limits<double>::infinity();
    if (p == 1)										// no chance of failure
        return (k == n) ? 0 : -std::numeric_limits<double>::infinity();

    return LogBinomFactor(n, k) + k * log(p) + (n - k) * log(1 - p);
}

// ---------------------------------------------------------
double BCMath::ApproxBinomial(unsigned n, unsigned k, double p)
{
    return exp(LogApproxBinomial(n, k, p));
}

namespace
{
/**
 * Vector of cached log factorials.
 * Hidden in anonymous namespace. */
static std::vector<double> CachedLogFactorials;
}

// ---------------------------------------------------------
double BCMath::LogBinomFactor(unsigned n, unsigned k)
{
    if (n < k)										// impose k < n requirement
        return std::numeric_limits<double>::quiet_NaN();

    if (k == 0 or k == n)
        return 0.;
    if (k == 1 or k == n - 1)
        return log((double)n);

    // if no approximation needed
    if ( n < CachedLogFactorials.size() and (n - k) < 10 )
        return LogBinomFactorExact(n, k);

    // calculate final log(n over k) using approximations if necessary
    return ApproxLogFact((double)n) - ApproxLogFact((double)k) - ApproxLogFact((double)(n - k));
}

// ---------------------------------------------------------
double BCMath::ApproxLogFact(double x)
{
    if (x < 0)
        return std::numeric_limits<double>::quiet_NaN();
    unsigned n = static_cast<unsigned>(floor(x + 0.5));
    if (n < CachedLogFactorials.size())
        return CachedLogFactorials[n];
    return x * log(x) - x + log(x * (1 + 4 * x * (1 + 2 * x))) / 6 + log(M_PI) / 2;
}

// ---------------------------------------------------------
double BCMath::LogFact(unsigned n)
{
    // return cached value if available
    if (n < CachedLogFactorials.size())
        return CachedLogFactorials[n];

    // calculate factorial starting from the highest cached value
    double logfact = (CachedLogFactorials.empty()) ? 0 : CachedLogFactorials.back();
    for (unsigned i = CachedLogFactorials.size(); i <= n; ++i)
        logfact += log((double)i);

    return logfact;
}

// ---------------------------------------------------------
unsigned BCMath::CacheFactorials(unsigned n)
{
    if (n < CachedLogFactorials.size())
        // log factorials have already been cached up to n
        return CachedLogFactorials.size();

    // reserve memory
    CachedLogFactorials.reserve(n);

    // add log(0!) if not there
    if (CachedLogFactorials.empty())
        CachedLogFactorials.push_back(0);

    // calculate new log factorials
    for (unsigned i = CachedLogFactorials.size(); i <= n; ++i)
        CachedLogFactorials.push_back(CachedLogFactorials.back() + log(static_cast<double>(i)));

    return CachedLogFactorials.size();
}

namespace
{
/**
 * Caches first factorials statically.
 * Hide in anonymous namespace. */
static unsigned cached_factorials = BCMath::CacheFactorials(899);
}

// ---------------------------------------------------------
int BCMath::Nint(double x)
{
    return static_cast<int>(((x > 0) ? 1 : -1) * floor(fabs(x) + 0.5));
}

// ---------------------------------------------------------
double BCMath::LogBinomFactorExact(unsigned n, unsigned k)
{
    if (n < k)										// impose k<n requirement
        return std::numeric_limits<double>::quiet_NaN();

    if ( k == 0 or k == n )
        return 0;
    if ( k == 1 or k == n - 1)
        return log((double)n);

    int lmax = std::max(k, n - k);
    int lmin = std::min(k, n - k);

    double log_bf = 0;

    for (int i = n; i > lmax; --i)
        log_bf += log((double)i);
    log_bf -= LogFact(lmin);

    return log_bf;
}

// ---------------------------------------------------------
double BCMath::LogBreitWignerNonRel(double x, double mean, double Gamma, bool norm)
{
    double bw = log(Gamma) - log((x - mean) * (x - mean) + Gamma * Gamma / 4.);

    if (norm)
        bw -= log(2. * M_PI);

    return bw;
}

// ---------------------------------------------------------

double BCMath::LogBreitWignerRel(double x, double mean, double Gamma)
{
    return -log((x * x - mean * mean) * (x * x - mean * mean) + mean * mean * Gamma * Gamma);
}

// ---------------------------------------------------------

double BCMath::LogChi2(double x, int n)
{
    if (x < 0) {
        BCLog::OutWarning("BCMath::LogChi2 : parameter cannot be negative!");
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (x == 0 && n == 1) {
        BCLog::OutWarning("BCMath::LogChi2 : returned value is infinity!");
        return std::numeric_limits<double>::infinity();
    }

    double nOver2 = ((double) n) / 2.;

    return (nOver2 - 1.) * log(x) - x / 2. - nOver2 * log(2) - log(TMath::Gamma(nOver2));
}

// ---------------------------------------------------------
double BCMath::LogVoigtian(double x, double sigma, double gamma)
{
    if (sigma <= 0 || gamma <= 0) {
        BCLog::OutWarning("BCMath::LogVoigtian : widths are negative or zero!");
        return std::numeric_limits<double>::quiet_NaN();
    }

    return log(TMath::Voigt(x, sigma, gamma));
}

// ---------------------------------------------------------
double BCMath::LogGammaPDF(double x, double alpha, double beta)
{
    return alpha * log(beta) - TMath::LnGamma(alpha) + (alpha - 1) * log(x) - beta * x;
}

// ---------------------------------------------------------
double BCMath::LogLogNormal(double x, double mean, double sigma)
{
    // if we have a delta function, return fixed value
    if (sigma == 0.)
        return 0;

    // if sigma is negative use absolute value
    if (sigma < 0.)
        sigma *= -1.;

    double arg = (log(x) - mean) / sigma;
    double result = -.5 * arg * arg;

    return result - log(x * sqrt(2. * M_PI) * sigma);
}

// ---------------------------------------------------------
double BCMath::CorrectPValue(const double& pvalue, const unsigned& npar, const unsigned& nobservations) throw (std::domain_error)
{
    // avoid pathologies
    if (pvalue < 0 or pvalue > 1)
        throw std::domain_error(Form("BCMath::CorrectPValue: pvalue (%g) out of range", pvalue));

    if (pvalue < std::numeric_limits<double>::epsilon())
        return 0;

    // dof = nobs - npar
    if (npar >= nobservations)
        throw std::domain_error(Form("BCMath::CorrectPValue: "
                                     "npar exceeds nobservations, %u vs %u", npar, nobservations));

    // convert upper-tail p value to a chi squared
    const double chi2 = ROOT::Math::chisquared_quantile_c(pvalue, nobservations);

    // corrected degrees of freedom
    const unsigned dof = nobservations - npar;

    // transform back to p value
    return TMath::Prob(chi2, dof);
}

// ---------------------------------------------------------
double BCMath::FastPValue(const std::vector<unsigned>& observed, const std::vector<double>& expected,
                          unsigned nIterations, unsigned seed) throw (std::invalid_argument)
{
    size_t nbins = observed.size();
    if (nbins != expected.size()) {
        throw std::invalid_argument(Form("BCMath::FastPValue: "
                                         "size of expected and observed do not match, %u vs %u", unsigned(expected.size()), unsigned(nbins)));
    }

    // temporary histogram to be modified in each MCMC step
    std::vector<unsigned> histogram(nbins, 0);

    // fix seed to iterations for reproducible results
    TRandom3 rng(seed);

    // keep track of log of probability and count data sets with larger value
    double logp = 0;
    double logp_start = 0;
    int counter_pvalue = 0;

    // define starting distribution as histogram with most likely entries
    for (size_t ibin = 0; ibin < nbins; ++ibin) {

        // get the number of expected events
        double yexp = expected[ibin];

        //most likely observed value = int(expected value)
        histogram[ibin] = size_t(yexp);

        // calculate test statistic (= likelihood of entire histogram) for starting distr.
        logp += LogPoisson(size_t(yexp), yexp);

        //statistic of the observed data set
        logp_start += LogPoisson(observed[ibin], yexp);
    }

    // loop over iterations
    for (unsigned iiter = 0; iiter < nIterations; ++iiter) {
        // loop over bins
        for (size_t ibin = 0; ibin < nbins; ++ibin) {
            // random step up or down in statistics for this bin
            double ptest = rng.Rndm() - 0.5;

            // increase statistics by 1
            if (ptest > 0) {
                // calculate factor of probability
                double r = expected[ibin] / double(histogram[ibin] + 1);

                // walk, or don't (this is the Metropolis part)
                if (rng.Rndm() < r) {
                    ++histogram[ibin];
                    logp += log(r);
                }
            }

            // decrease statistics by 1
            else if (ptest <= 0 && histogram[ibin] > 0) {
                // calculate factor of probability
                double r = double(histogram[ibin]) / expected[ibin];

                // walk, or don't (this is the Metropolis part)
                if (rng.Rndm() < r) {
                    --histogram[ibin];
                    logp += log(r);
                }
            }
        } // end of looping over bins

        // increase counter
        if (logp <= logp_start)
            ++counter_pvalue;
        // handle the case where a -b +b > a because of precision loss
        else if (logp_start && logp != logp_start && fabs((logp - logp_start) / logp_start) < 1e-15)
            ++counter_pvalue;

    } // end of looping over iterations

    // calculate p-value
    return double(counter_pvalue) / nIterations;
}

double BCMath::Random::Chi2(TRandom* rng, double dof)
{
    return 2.0 * Gamma(rng, dof / 2.0, 1.0);
}

double BCMath::Random::Gamma(TRandom* rng, double a, double b)
{
    /* assume a > 0 */

    if (a < 1) {
        // u in ]0,1[
        double u = rng->Uniform(1);
        return Gamma(rng, 1.0 + a, b) * std::pow(u, 1.0 / a);
    }

    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    while (1) {
        do {
            // unit Gauss
            x = rng->Gaus();
            v = 1.0 + c * x;
        } while (v <= 0);

        v = v * v * v;
        u = rng->Uniform(1);

        if (u < 1 - 0.0331 * x * x * x * x)
            break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
            break;
    }

    return b * d * v;
}
