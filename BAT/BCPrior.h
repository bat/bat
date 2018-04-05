#ifndef __BCPRIOR__H
#define __BCPRIOR__H

/*!
 * \class BCPrior
 * \brief A class to represent the prior of a parameter
 * \author Daniel Greenwald
 * \version 1.0
 * \date 01.2015
 *
 * \defgroup Priors Classes for defining priors
 *
 * If the prior factorizes into a product of 1D priors, members of
 * this group can be used to model the individual 1D
 * distributions. For example, a standard Gaussian for the first parameter is set by
 *
 * ~~~{.cpp}
 * model.GetParameter(0).SetPrior(new BCPriorGaussian(1, 0))
 * ~~~
 */


/*
 * Copyright (C) 2007-2018, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------


#include "BCH1D.h"
#include "BCH2D.h"

#include <TF1.h>

#include <cstddef>
#include <cmath>
#include <limits>

class TH1;
class TH2;
class TRandom;

// ---------------------------------------------------------

class BCPrior
{

public:

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * Empty constructor */
    BCPrior();

    /**
     * Copy constructor */
    BCPrior(const BCPrior& other);

    /**
     * Destructor */
    virtual ~BCPrior();

    /** @} */

    /** \name operators and swap */
    /** @{ */

    /** swap */
    friend void swap(BCPrior& A, BCPrior& B);

    /** @} */

    /** \name methods that must be overloaded in derived classes */
    /** @{ */

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x) = 0;

    /**
     * Clone function. [Copy constructor must also be provided.] */
    virtual BCPrior* Clone() const = 0;
    // { return new [Derived Class](*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const = 0;

    /** @} */

    /** \name Getters */
    /** @{ **/

    /**
     * Get log of normalized prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogNormalizedPrior(double x)
    { return GetLogPrior(x) - fLogIntegral;}

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @param normalize Whether to normalize prior with stored integral
     * @return prior */
    virtual double GetPrior(double x, bool normalize = false);

    /**
     * Return back ROOT TF1 evaluating BCPrior::GetPrior */
    virtual TF1& GetFunction()
    { return fPriorFunction; }

    /**
     * Return back ROOT TF1 evaluating BCPrior::GetPrior */
    virtual const TF1& GetFunction() const
    { return fPriorFunction; }

    /**
     * Set range of ROOT TF1 function. */
    virtual void SetFunctionRange(double xmin, double xmax);

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return raw moment of prior distribution */
    virtual double GetRawMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get integral of prior.
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return integral of prior */
    virtual double GetIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get central moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return central of prior distribution */
    virtual double GetCentralMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get standardised moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standardised moment of prior distribution */
    virtual double GetStandardizedMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity());

    /**
     * Get mean of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mean of prior distribution */
    virtual double GetMean(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return GetRawMoment(1, xmin, xmax); }

    /**
     * Get variance of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return variance of prior distribution */
    virtual double GetVariance(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return GetCentralMoment(2, xmin, xmax); }

    /**
     * Get standard deviation of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standard deviation of prior distribution */
    virtual double GetStandardDeviation(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return sqrt(fabs(GetVariance(xmin, xmax))); }

    /**
     * Get skewness of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return skewness of prior distribution */
    virtual double GetSkewness(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return GetStandardizedMoment(3, xmin, xmax); }

    /**
     * Get kurtosis of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return kurtosis of prior distribution */
    virtual double GetKurtosis(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { return GetStandardizedMoment(4, xmin, xmax); }

    /**
     * @return a random value distributed according to the prior.
     * @param xmin lower limit of range to generate value in
     * @param xmax upper limit of range to generate value in
     * @param R Pointer to the random generator to be used, if needed.
     * @return random value. */
    virtual double GetRandomValue(double xmin, double xmax, TRandom* const R = NULL);

    /** @} **/
    /** \name Functions related to ROOT access. */
    /** @{ **/

    /**
     * For accessing prior as ROOT TF1 */
    virtual double GetPriorForROOT(double* x, double* /*p*/)
    { return GetPrior(x[0]); }

    /**
     * For accessing normalized prior as ROOT TF1 */
    virtual double GetNormalizedPriorForROOT(double* x, double* /*p*/)
    { return GetPrior(x[0]); }

    /**
     * For accessing log(prior) as ROOT TF1 */
    virtual double GetLogPriorForROOT(double* x, double* /*p*/)
    { return GetLogPrior(x[0]); }

    /**
     * For accessing normalized log(prior) as ROOT TF1 */
    virtual double GetNormalizedLogPriorForROOT(double* x, double* /*p*/)
    { return GetLogPrior(x[0]) - fLogIntegral; }

    /**
     * Calculate and store integral for use in normalized TF1s */
    virtual double CalculateAndStoreIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { fLogIntegral = log(GetIntegral(xmin, xmax)); return fLogIntegral;}

    /**
     * Store integral; */
    virtual void StoreIntegral(double I)
    { fLogIntegral = log(I); }

    /**
     * Store log(integral); */
    virtual void StoreLogIntegral(double logI)
    { fLogIntegral = logI; }

    /**
     * Get stored integral. */
    virtual double GetStoredIntegral() const
    { return exp(fLogIntegral); }

    /**
     * Get stored integral. */
    virtual double GetStoredLogIntegral() const
    { return fLogIntegral; }

    /**
     * Fill histogram by prior evaluated at bin center. */
    virtual void FillHistogramByCenterValue(TH1* h);

    /**
     * Fill histogram by integrating prior over bin and dividing by bin width. */
    virtual void FillHistogramByIntegral(TH1* h);

    /**
     * Get BCH1D object for prior.
     * @param bins pointer to TH1 object defining binning to use (and axis names)
     * @param name name to assigned to the histogram created for BCH1D object
     * @return BCH1D object for prior. */
    virtual BCH1D GetBCH1D(TH1* bins, const std::string& name = "prior");

    /**
     * Get BCH2D object for prior.
     * @param ordinate pointer to BCPrior object to be the ordinate axis
     * @param bins pointer to TH2 object defining binning to use (and axis names)
     * @param name name to give histogram created for BCH2D object
     * @return BCH2D object for prior. */
    virtual BCH2D GetBCH2D(BCPrior* ordinate, TH2* bins, const std::string& name = "prior");

    /** @} **/

protected:
    TF1 fPriorFunction; ///< TF1 for use in default raw moment calculation

    double fLogIntegral; ///< Log of integral of unnormalized pdf over the range.
};

#endif
