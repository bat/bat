#ifndef __BCH1D__H
#define __BCH1D__H

/**
 * @class BCH1D
 * @brief A class for handling 1D distributions.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @author Daniel Greenwald
 * @version 1.0
 * @date 08.2008
 * @details This class contains a TH1 histogram and some additional
 * functions. It is used for marginalized distributions.
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include "BCHistogramBase.h"

#include <vector>

class TH1;

// ---------------------------------------------------------

class BCH1D : public BCHistogramBase
{

public:

    /** \name Enumerators */
    /** @{ */

    /**
     * Enum for type of bands to be drawn on plot. */
    enum BCH1DBandType {
        kNoBands          = -1,
        kCentralInterval  =  0,
        kSmallestInterval =  1,
        kUpperLimit       =  2,
        kLowerLimit       =  3,
        kUserSpecified    =  4
    };

    /** @} */

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * The default constructor. */
    BCH1D(const TH1* const hist = 0);

    /**
     * Copy constructor. */
    BCH1D(const BCH1D& other);

    /**
     * The default destructor. */
    virtual ~BCH1D() {};

    /** @} */

    /** \name Member functions (get)  */
    /** @{ */

    using BCHistogramBase::GetLocalMode;
    /**
     * @return local mode. */
    double GetLocalMode()
    { return GetLocalMode(0); }

    using BCHistogramBase::GetBestFitParameters;
    /**
     * @return global mode. */
    double GetBestFitParameters()
    { return GetBestFitParameters(0); }

    /**
     * @return The median of the distribution. */
    double GetMedian()
    { return this->GetQuantile(0.5); };

    /**
     * Returns the quantile of the distribution.
     * @param probability The probability.
     * @return The quantile of the distribution for the probability.
     * @see GetLimit(double probability) */
    double GetQuantile(double probability);

    /**
     * Return the quantile of the distribution
     * @param probability The probability.
     * @return The quantile of the distribution for the probability.
     * @see GetQuantile(double probablity) */
    double GetLimit(double probability)
    { return this->GetQuantile(probability); };

    /**
     * @return Band type. */
    BCH1DBandType GetBandType()
    { return fBandType; }

    /**
     * @return Number of quantiles to draw. */
    unsigned GetNQuantiles()
    { return fNQuantiles; }

    /**
     * @return Quantile line color. */
    int GetQuantileLineColor()
    { return fQuantileLineColor; }

    /**
     *@return whether to draw median. */
    bool GetDrawMedian()
    { return fDrawMedian; }

    /**
     * @return whether to draw central 68% interval. */
    bool GetDrawCentral68()
    { return fDrawCentral68; }

    /** @} */

    /** \name Member functions (set)  */
    /** @{ */

    using BCHistogramBase::CopyOptions;

    /**
     * Copy options from other. */
    void CopyOptions(const BCH1D& other);

    /**
     * Sets the color scheme. */
    void SetColorScheme(BCHColorScheme scheme);

    using BCHistogramBase::SetGlobalMode;

    /**
     * Set global mode */
    void SetGlobalMode(double mode)
    { SetGlobalMode(std::vector<double>(1, mode)); }

    using BCHistogramBase::SetLocalMode;

    /**
     * Set local mode */
    void SetLocalMode(double mode)
    { SetLocalMode(std::vector<double>(1, mode)); }

    /**
     * Set band type. */
    void SetBandType(BCH1DBandType bt)
    { fBandType = bt; }

    /**
     * Set draw quantiles.
     * @param n N divisions of quantiles to draw, set to zero or one to disable drawing of quantiles. */
    void SetDrawQuantiles(unsigned n)
    { fNQuantiles = n; }

    /**
     * Set quantile line color.
     * @param c Quantile line color. */
    void SetQuantileLineColor(int c)
    { fQuantileLineColor = c; }

    /**
     * Set drawing of median.
     * @param flag Toggles drawing of median.
     * @param central68 Toggles drawing of arrows for central 68% interval. (Automatically suppressed if median is suppressed.)*/
    void SetDrawMedian(bool flag = true, bool central68 = true)
    { fDrawMedian = flag; fDrawCentral68 = central68;}

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    using BCHistogramBase::CheckIntervals;

    /**
     * Check intervals: remove values below 0 or above 1,
     * and sort to proper order for band type. */
    virtual void CheckIntervals(std::vector<double>& intervals);

    /**
     * Return default intervals.
     * @param nbands nbands to give defaults for; if negative, use fNBands.
     * @return vector of default values for band intervals. */
    virtual std::vector<double> DefaultIntervals(int nbands = -1);

    /**
     * Draw bands. */
    virtual void DrawBands(const std::string& options = "same");

    /**
     * Draw markers: global mode, local mode, mean, quantiles, median. */
    virtual void DrawMarkers();

    /**
     * Draw quantiles. */
    virtual void DrawQuantiles(unsigned n);

    /**
     * Draw median & central 68% interval. */
    virtual void DrawMedian();

    /**
     * Print information to log
     * @param prefix String to be prepended to every line.
     * @param intervals Vector of intervals to print.
     * @param prec Precision of doubles to output. */
    void PrintSummary(const std::string& prefix = "", unsigned prec = 6, std::vector<double> intervals = std::vector<double>(0));

    /**
     * \struct BCH1DInterval
     * Contains information about an interval. */
    struct BCH1DInterval {
        double xmin;
        double xmax;
        double mode;
        double relative_height;
        double relative_mass;

        BCH1DInterval();
    };

    /**
     * \struct BCH1DIntervals
     * Vector of intervals with information about total mass. */
    struct BCH1DSmallestInterval {
        std::vector<BCH1D::BCH1DInterval> intervals;
        double total_mass;
        double mode;
        double max_val;

        BCH1DSmallestInterval();
    };

    /**
     * @param masses Masses to return smallest intervals of.
     * @return vector of BCH1DSmallestInterval's. */
    std::vector<BCH1D::BCH1DSmallestInterval> GetSmallestIntervals(std::vector<double> masses);

    /**
     * @param mass Probability mass to return smallest intervals of.
     * @return BCH1DSmallestInterval for probability mass. */
    BCH1D::BCH1DSmallestInterval GetSmallestIntervals(double mass)
    { return GetSmallestIntervals(std::vector<double>(1, mass)).at(0); }

    /**
     * Get histogram with bins outside min, max band being zero. The
     * new histogram can have 2 more bins than the original one as the
     * bins where min and max fall into will be split in two (except for the
     * case when min and/or max are equal to some of the original bin
     * boundaries.
     * @param min lower boundary of the non-zero interval
     * @param max upper boundary of the non-zero interval
     * @param name Name for new histogram; empty string (default) appends "subhist" to histogram name.
     * @param preserve_range If true, preserves original histograms range, setting bins outside subhistogram range to zero.
     * @return new histogram which is nonzero only between min and max */
    TH1* GetSubHistogram(double min, double max, const std::string& name = "", bool preserve_range = false);

    /** @} */

protected:

    /**
     * Band type */
    BCH1DBandType fBandType;

    /**
     * Number of quantiles to draw. */
    unsigned fNQuantiles;

    /**
     * Quantile line color. */
    int fQuantileLineColor;

    /**
     * flag for drawing median. */
    bool fDrawMedian;

    /**
     * flag for darwing central 68% interval arrows. */
    bool fDrawCentral68;
};

// ---------------------------------------------------------

#endif
