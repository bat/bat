/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCMath.h"
#include "BCLog.h"

#include <math.h>
#include <limits>

#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>

#include <Math/PdfFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>

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
    return result - log(sqrt(2 * M_PI) * sigma);
}

// ---------------------------------------------------------
double BCMath::LogPoisson(double x, double lambda)
{
    if (x < 0)
        return -std::numeric_limits<double>::infinity();

    if (lambda <= 0)
        return -std::numeric_limits<double>::infinity();

    if (lambda > 899)
        return LogGaus(x, lambda, sqrt(lambda), true);

    if (x == 0.)
        return -lambda;

    return x * log(lambda) - lambda - ApproxLogFact(x);
}

// ---------------------------------------------------------
double BCMath::LogApproxBinomial(unsigned n, unsigned k, double p)
{
    if (k < n)										// impose k < n requirement
        return LogApproxBinomial(k, n, p);

    // p in [0,1]
    if (p < 0 or p > 1)
        return -std::numeric_limits<double>::infinity();

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
        return LogBinomFactor(k, n);

    if (k == 0 or k == n)
        return 0.;
    if (k == 1 or k == n - 1)
        return log((double)n);

    // if no approximation needed
    if ( n < BCMATH_NFACT_ALIMIT or (n < CachedLogFactorials.size() and  (n - k) < 10) )
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
        return LogBinomFactorExact(k, n);

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
        return -std::numeric_limits<double>::infinity();
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
        return -std::numeric_limits<double>::infinity();
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
double BCMath::chi2(double* x, double* par)
{
    return ROOT::Math::chisquared_pdf(x[0], par[0]);
}

// ---------------------------------------------------------
void BCMath::RandomChi2(std::vector<double>& randoms, int K)
{
    // fixed upper cutoff to 1000, might be too small
    TF1* f = new TF1("chi2", chi2, 0.0, 1000, 1);
    f->SetParameter(0, K);
    f->SetNpx(500);
    // uses inverse-transform method
    // fortunately CDF only built once
    for (unsigned int i = 0; i < randoms.size(); i++)
        randoms.at(i) = f->GetRandom();

    delete f;
}

// ---------------------------------------------------------
TH1D* BCMath::ECDF(const std::vector<double>& data)
{

    // copy data to new vector to eventually contain bin lower edges:
    std::vector<double> bins = data;
    // sort to increasing order
    std::stable_sort(bins.begin(), bins.end(), std::less<double>());
    // remove duplicates
    std::vector<double>::iterator last = std::unique(bins.begin(), bins.end());
    bins.erase(last, bins.end());

    // create histogram where
    // lower edge of first bin = min. data
    // upper edge of last bin = max. data
    TH1D* ECDF = new TH1D("ECDF", "Empirical cumulative distribution function", bins.size() - 1, &bins[0]);

    // fill the data in to find multiplicities
    for (unsigned i = 0; i < data.size(); ++i)
        ECDF->Fill(data[i]);

    // divide first bin by number of data points
    ECDF->SetBinContent(1, ECDF->GetBinContent(1) / data.size());
    ECDF->SetBinError(1, 0.0);

    // construct the ecdf
    for (int nBin = 2; nBin <= ECDF->GetNbinsX(); nBin++) {
        // new content = prior bin content (already summed up and divided by n) + this bin content / n
        ECDF->SetBinContent(nBin, ECDF->GetBinContent(nBin - 1) + ECDF->GetBinContent(nBin) / data.size());
        ECDF->SetBinError(nBin, 0.0);
    }

    // adjust for nice plotting
    ECDF->SetMinimum(0.);
    ECDF->SetMaximum(1.);

    return ECDF;
}

// ---------------------------------------------------------

std::vector<int> BCMath::longestRuns(const std::vector<bool>& bitStream)
{
    // initialize counter variables
    unsigned int maxRunAbove, maxRunBelow, currRun;
    maxRunAbove = 0;
    maxRunBelow = 0;
    currRun = 1;
    // set both entries to zero
    std::vector<int> runs(2, 0);

    if (bitStream.empty())
        return runs;

    // flag about kind of the currently considered run
    bool aboveRun = bitStream.at(0);

    // start at second variable
    std::vector<bool>::const_iterator iter = bitStream.begin();
    ++iter;
    while (iter != bitStream.end()) {

        // increase counter if run continues
        if (*(iter - 1) == *iter)
            currRun++;
        else {
            // compare terminated run to maximum
            if (aboveRun)
                maxRunAbove = std::max(maxRunAbove, currRun);
            else
                maxRunBelow = std::max(maxRunBelow, currRun);
            // set flag to run of opposite kind
            aboveRun = !aboveRun;
            // restart at length one
            currRun = 1;
        }
        // move to next bit
        ++iter;
    }

    // check last run
    if (aboveRun)
        maxRunAbove = std::max(maxRunAbove, currRun);
    else
        maxRunBelow = std::max(maxRunBelow, currRun);

    // save the longest runs
    runs.at(0) = maxRunBelow;
    runs.at(1) = maxRunAbove;

    return runs;
}
// ---------------------------------------------------------

std::vector<double> BCMath::longestRunsChi2(
    const std::vector<double>& yMeasured,
    const std::vector<double>& yExpected, const std::vector<double>& sigma)
{
    //initialize counter variables
    double maxRunAbove, maxRunBelow, currRun;
    maxRunAbove = 0;
    maxRunBelow = 0;
    currRun = 0;
    //set both entries to zero
    std::vector<double> runs(2, 0);

    //check input size
    if (yMeasured.size() != yExpected.size() || yMeasured.size() != sigma.size()
            || yExpected.size() != sigma.size()) {
        //should throw exception
        return runs;
    }

    //exclude zero uncertainty
    //...

    int N = yMeasured.size();
    if ( N <= 0)
        return runs;
    //BCLog::OutDebug(Form("N = %d", N));


    //flag about kind of the currently considered run
    double residue = (yMeasured.at(0) - yExpected.at(0)) / sigma.at(0);
    bool aboveRun = residue >= 0 ? true : false;
    currRun = residue * residue;

    //start at second variable
    for (int i = 1; i < N; i++) {
        residue = (yMeasured.at(i) - yExpected.at(i)) / sigma.at(i);
        //run continues
        if ((residue >= 0) == aboveRun) {
            currRun += residue * residue;
        } else {
            //compare terminated run to maximum
            if (aboveRun)
                maxRunAbove = std::max(maxRunAbove, currRun);
            else
                maxRunBelow = std::max(maxRunBelow, currRun);
            //set flag to run of opposite kind
            aboveRun = !aboveRun;
            //restart at current residual
            currRun = residue * residue;
        }
    }

    //check last run
    if (aboveRun)
        maxRunAbove = std::max(maxRunAbove, currRun);
    else
        maxRunBelow = std::max(maxRunBelow, currRun);

    //save the longest runs
    runs.at(0) = maxRunBelow;
    runs.at(1) = maxRunAbove;

    return runs;
}
// ---------------------------------------------------------
double BCMath::longestRunFrequency(unsigned longestObserved, unsigned int nTrials)
{
    // can't observe run that's longer than the whole sequence
    if (longestObserved >= nTrials)
        return 0.;

    // return value
    double prob = 0.;

    // short cuts
    typedef unsigned int uint;
    uint Lobs = longestObserved;
    uint n = nTrials;

    // the result of the inner loop is the cond. P given r successes
    double conditionalProb;

    /* first method: use the gamma function for the factorials: bit slower and more inaccurate
     * in fact may return NaN for n >= 1000.
     * alternative using log factorial approximations, is faster and more accurate
     */

    double tempLog = 0;
    for (uint r = 0; r <= n; r++) {
        conditionalProb = 0.0;

        for (uint i = 1; (i <= n - r + 1) && (i <= uint(r / double(Lobs + 1))); i++) {
            tempLog = ApproxLogFact(n - i * (Lobs + 1)) - ApproxLogFact(i)
                      - ApproxLogFact(n - r - i + 1)
                      - ApproxLogFact(r - i * (Lobs + 1));
            if (i % 2)
                conditionalProb += exp(tempLog);
            else
                conditionalProb -= exp(tempLog);
        }
        prob += (1 + n - r) * conditionalProb;
    }

    // Bernoulli probability of each permutation
    prob *= pow(2., -double(n));

    return prob;
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
