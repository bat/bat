#ifndef __BCPRIOR__H
#define __BCPRIOR__H

/*!
 * \class BCPrior
 * \brief A class to represent the prior of a parameter
 * \author Daniel Greenwald
 * \version 1.0
 * \date 01.2015
 */


/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <limits>
#include <cmath>
#include <cstddef>

class TF1;
class TH1;
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
     * Constructor explicitely setting fPriorFunction.
     * @param f TF1 to be copied into prior. */
    BCPrior(const TF1* const f);

    /**
     * Copy constructor */
    BCPrior(const BCPrior& other);

    /**
     * Destructor */
    virtual ~BCPrior();

    /** @} */

    /** \name methods that must be overloaded in derived classes */
    /** @{ */

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x) const = 0;

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
     * Get prior
     * @param x value to evaluate prior at
     * @return prior */
    virtual double GetPrior(double x) const;

    /**
     * Return back ROOT TF1 evaluating BCPrior::GetPrior */
    virtual TF1* GetFunction()
    { return fPriorFunction; }

    /**
     * Set range of ROOT TF1 function. */
    virtual void SetFunctionRange(double xmin, double xmax);

    /**
     * Return mode of prior (in range).
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mode of prior in range. */
    virtual double GetMode(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get raw moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return raw moment of prior distribution */
    virtual double GetRawMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get integral of prior.
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return integral of prior */
    virtual double GetIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get central moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return central of prior distribution */
    virtual double GetCentralMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get standardised moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standardised moment of prior distribution */
    virtual double GetStandardizedMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get mean of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return mean of prior distribution */
    virtual double GetMean(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { return GetRawMoment(1, xmin, xmax); }

    /**
     * Get variance of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return variance of prior distribution */
    virtual double GetVariance(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { return GetCentralMoment(2, xmin, xmax); }

    /**
     * Get standard deviation of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standard deviation of prior distribution */
    virtual double GetStandardDeviation(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { return sqrt(fabs(GetVariance(xmin, xmax))); }

    /**
     * Get skewness of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return skewness of prior distribution */
    virtual double GetSkewness(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { return GetStandardizedMoment(3, xmin, xmax); }

    /**
     * Get kurtosis of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return kurtosis of prior distribution */
    virtual double GetKurtosis(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { return GetStandardizedMoment(4, xmin, xmax); }

    /**
     * @return a random value distributed according to the prior.
     * @param xmin lower limit of range to generate value in
     * @param xmax upper limit of range to generate value in
     * @param R Pointer to the random generator to be used, if needed.
     * @return random value. */
    virtual double GetRandomValue(double xmin, double xmax, TRandom* const R = NULL) const;

    /** @} **/
    /** \name Functions related to ROOT access. */
    /** @{ **/

    /**
     * For accessing prior as ROOT TF1 */
    double GetPriorForROOT(double* x, double* /*p*/)
    { return GetPrior(x[0]); }

    /**
     * For accessing normalized prior as ROOT TF1 */
    double GetNormalizedPriorForROOT(double* x, double* /*p*/)
    { return GetPrior(x[0]); }

    /**
     * For accessing log(prior) as ROOT TF1 */
    double GetLogPriorForROOT(double* x, double* /*p*/)
    { return GetLogPrior(x[0]) / fIntegral; }

    /**
     * For accessing normalized log(prior) as ROOT TF1 */
    double GetNormalizedLogPriorForROOT(double* x, double* /*p*/)
    { return GetLogPrior(x[0]) - log(fIntegral); }

    /**
     * Calculate and store integral for use in normalized TF1s */
    double CalculateAndStoreIntegral(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity())
    { fIntegral = GetIntegral(xmin, xmax); return fIntegral;}

    /**
     * Store integral; */
    void StoreIntegral(double I)
    { fIntegral = I; }

    /**
     * Get stored integral. */
    double GetStoredIntegral() const
    { return fIntegral; }

    /**
     * Fill histogram by prior evaluated at bin center. */
    void FillHistogramByCenterValue(TH1* h);

    /**
     * Fill histogram by integrating prior over bin and dividing by bin width. */
    void FillHistogramByIntegral(TH1* h);

    /** @} **/

protected:
    TF1* fPriorFunction;					///< TF1 for use in default raw moment calculation

    double fIntegral;                  ///< Integral of unnormalized pdf over the range.
};

#endif
