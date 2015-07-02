#ifndef __BCTH1PRIOR__H
#define __BCTH1PRIOR__H

/*!
 * \class BCTH1Prior
 * \brief A class to represent the prior of a parameter by a TH1
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

#include "BCPrior.h"

#include <limits>

#include <TH1.h>

// ---------------------------------------------------------

class BCTH1Prior : public BCPrior
{
public:
    /** \name Constructor & Destructor */
    /** @{ **/

    /** Constructor */
    BCTH1Prior(TH1& h, bool interpolate = false);

    /** Constructor */
    BCTH1Prior(TH1* h, bool interpolate = false);

    /** Copy constructor */
    BCTH1Prior(const BCTH1Prior& other);

    /** Destrcutor */
    virtual ~BCTH1Prior();

    /** @} **/

    /** \name operator and swap */
    /** @{ **/

    /** assignment operator */
    BCTH1Prior& operator=(const BCTH1Prior& rhs);

    /** swap */
    friend void swap(BCTH1Prior& A, BCTH1Prior& B);

    /** @} **/

    /** \name Functions overloaded from BCPrior **/
    /** @{ **/

    /** Clone function */
    virtual BCPrior* Clone() const
    { return new BCTH1Prior(*this); }

    /**
     * @return Whether everything needed for prior is set and prior can be used. */
    virtual bool IsValid() const;

    /**
     * Get prior
     * @param x value to evaluate prior at
     * @param normalize Whether to normalize prior with stored integral
     * @return prior */
    virtual double GetPrior(double x, bool normalize = false) const
    { return ((fInterpolate) ? const_cast<TH1*>(&fPriorHistogram)->Interpolate(x) : const_cast<TH1*>(&fPriorHistogram)->GetBinContent(const_cast<TH1*>(&fPriorHistogram)->FindFixBin(x))) * ((normalize) ? exp(-fLogIntegral) : 1); }

    /**
     * Get log of prior
     * @param x value to evaluate log of prior at
     * @return log of prior */
    virtual double GetLogPrior(double x) const;

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
     * Get standardised moment of prior distrubion. If limits are infinite, use exact value from prior type.
     * @param n moment number
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standardised moment of prior distribution */
    virtual double GetStandardizedMoment(unsigned n, double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const;

    /**
     * Get variance of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return variance of prior distribution */
    virtual double GetVariance(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { double s = GetStandardDeviation(xmin, xmax); return s * s; }

    /**
     * Get standard deviation of prior. If limits are infinite, use exact value from prior type
     * @param xmin lower limit of range to evaluate over
     * @param xmax upper limit of range to evaluate over
     * @return standard deviation of prior distribution */
    virtual double GetStandardDeviation(double xmin = -std::numeric_limits<double>::infinity(), double xmax = std::numeric_limits<double>::infinity()) const
    { (void)xmin; (void)xmax; return fPriorHistogram.GetRMS(); }

    /**
     * Get BCH1D object for prior.
     * @param h pointer to TH1 object defining binning to use.
     * @param name name to give histogram created for BCH1D object
     * @return BCH1D object for prior. */
    virtual BCH1D GetBCH1D(TH1* bins, const char* name = "prior");

    /** @} **/

    /** \name Setters */
    /** @{ **/

    virtual void SetInterpolate(bool interpolate)
    { fInterpolate = interpolate; }

    /**
     * @return a random value distributed according to the prior.
     * @param xmin lower limit of range to generate value in
     * @param xmax upper limit of range to generate value in
     * @param R Pointer to the random generator to be used, if needed.
     * @return random value. */
    virtual double GetRandomValue(double /*xmin*/, double /*xmax*/, TRandom* const /*R*/ = NULL) const
    { return const_cast<TH1*>(&fPriorHistogram)->GetRandom(); }

    /** @} **/

    /** \name Getters */
    /** @{ **/

    virtual TH1& GetHistogram()
    { return fPriorHistogram; }

    virtual const TH1& GetHistogram() const
    { return fPriorHistogram; }

    virtual bool GetInterpolate()
    { return fInterpolate; }

    /** @} **/

    /** \name Misc */
    /** @{ **/

    /** Normalize the histogram holding the prior. */
    void NormalizeHistogram();

    /** @} **/

protected:

    TH1& fPriorHistogram;  //< TH1 holding prior

    bool fInterpolate; //< whether to interpolate values of hist for prior function
};

#endif
