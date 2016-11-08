#ifndef __BCH2D__H
#define __BCH2D__H

/**
 * @class BCH2D
 * @brief  A class for handling 2D distributions.
 * @author Daniel Kollar
 * @author Kevin Kr&ouml;ninger
 * @author Daniel Greenwald
 * @version 1.0
 * @date 08.2008
 * @details This class contains a TH2 histogram and some additional
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
#include <string>

class TH2;
class TGraph;
class TObject;

// ---------------------------------------------------------

class BCH2D : public BCHistogramBase
{

public:

    /** \name Enumerators */
    /** @{ */

    /**
     * Enum for type of bands to draw on plot. */
    enum BCH2DBandType {
        kNoBands                  = -1, //!< no bands
        kSmallestInterval          = 0  //!< smallest intervals containing probability mass
                                     // kCentralIntervalOfYGivenX  = 1, //!< 1D central intervals in Y for each bin in X
                                     // kCentralIntervalOfXGivenY  = 2, //!< 1D central intervals in X for each bin in Y
                                     // kSmallestIntervalOfYGivenX = 3, //!< 1D smallest intervals in Y for each bin in X
                                     // kSmallestIntervalOfXGivenY = 4, //!< 1D smallest intervals in X for each bin in Y
                                     // kUpperLimitOfYGivenX       = 5, //!< 1D upper limits of Y for each bin in X
                                     // kUpperLimitOfXGivenY       = 6, //!< 1D upper limits of X for each bin in Y
                                     // kLowerLimitOfYGivenX       = 7, //!< 1D lower limits of Y for each bin in X
                                     // kLowerLimitOfXGivenY       = 8  //!< 1D lower limits of X for each bin in Y
    };

    /**
     * Enum for type of profile. */
    enum BCH2DProfileType {
        kProfileMean   = 0,
        kProfileMedian = 1,
        kProfileMode   = 2
    };

    /**
     * Enum for axis to profile. */
    enum BCH2DProfileAxis {
        kProfileX = 0,
        kProfileY = 1
    };

    /** @} */

    /** \name Constructor and destructors */
    /** @{ */

    /**
     * The complete constructor. */
    BCH2D(const TH2* const h = 0);

    /**
     * Copy constuctor. */
    BCH2D(const BCH2D& other);

    /**
     * The default destructor. */
    virtual ~BCH2D() {};

    /** @} */

    /** \name Member functions (get)  */
    /** @{ */

    /**
     * Return the TH2 histogram
     * @return The TH2 histogram. */
    TH2* GetHistogram()
    { return (TH2*)fHistogram; }

    /**
     * @return Band type. */
    BCH2DBandType GetBandType()
    { return fBandType; }

    /**
     * @return flag for plotting log on z axis. */
    bool GetLogz()
    { return fLogz; }

    /**
     * @return flag for drawing x profile. */
    bool GetDrawProfileX()
    { return fDrawProfileX; }

    /**
     * @return profile type for x profile. */
    BCH2DProfileType GetProfileXType()
    { return fProfileXType; }

    /**
     * @return line color for x profile. */
    int GetProfileXLineColor()
    { return fProfileXLineColor; }

    /**
     * @return line style for x profile. */
    int GetProfileXLineStyle()
    { return fProfileXLineStyle; }

    /**
     * @return flag for drawing y profile. */
    bool GetDrawProfileY()
    { return fDrawProfileY; }

    /**
     * @return profile type for y profile. */
    BCH2DProfileType GetProfileYType()
    { return fProfileYType; }

    /**
     * @return line color for y profile. */
    int GetProfileYLineColor()
    { return fProfileYLineColor; }

    /**
     * @return line style for y profile. */
    int GetProfileYLineStyle()
    { return fProfileYLineStyle; }

    /** @} */


    /** \name Member functions (set)  */
    /** @{ */

    using BCHistogramBase::CopyOptions;

    /**
     * copy options from */
    void CopyOptions(const BCH2D& other);

    using BCHistogramBase::SetGlobalMode;

    /**
     * Set global mode.
     * @param x Global mode in x.
     * @param y Global mode in y. */
    void SetGlobalMode(double x, double y)
    { std::vector<double> m(1, x); m.push_back(y); SetGlobalMode(m); }

    using BCHistogramBase::SetLocalMode;

    /**
     * Set local mode.
     * @param x Global mode in x.
     * @param y Global mode in y. */
    void SetLocalMode(double x, double y)
    { std::vector<double> m(1, x); m.push_back(y); SetLocalMode(m); }

    /** Set band type. */
    void SetBandType(BCH2DBandType bt)
    { fBandType = bt; }

    /**
     * Sets drawing of z axis in log. */
    void SetLogz(bool flag = true)
    { fLogz = flag; }

    /**
     * Set drawing of x profile. */
    void SetDrawProfileX(bool flag = true)
    { fDrawProfileX = flag; }

    /**
     * Set profile type of x profile. */
    void SetProfileXType(BCH2DProfileType pt)
    { fProfileXType = pt; }

    /**
     * Set line color of x profile. */
    void SetProfileXLineColor(int c)
    { fProfileXLineColor = c; }

    /**
     * Set line style of x profile. */
    void SetProfileXLineStyle(int s)
    { fProfileXLineStyle = s; }

    /**
     * Set drawing of y profile. */
    void SetDrawProfileY(bool flag = true)
    { fDrawProfileY = flag; }

    /**
     * Set profile type of y profile. */
    void SetProfileYType(BCH2DProfileType pt)
    { fProfileYType = pt; }

    /**
     * Set line color of y profile. */
    void SetProfileYLineColor(int c)
    { fProfileYLineColor = c; }

    /**
     * Set line style of y profile. */
    void SetProfileYLineStyle(int s)
    { fProfileYLineStyle = s; }

    /** @} */


    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    using BCHistogramBase::CheckIntervals;

    /**
     * Check intervals: remove values below 0 or above 1,
     * and sort to proper order band type. */
    virtual void CheckIntervals(std::vector<double>& intervals);

    /**
     * Return default intervals.
     * @param nbands nbands to give defaults for; if negative, use fNBands.
     * @return vector of default values for band intervals. */
    virtual std::vector<double> DefaultIntervals(int nbands = -1);

    /**
     * Draw band, or if band type set to no bands, histogram. */
    virtual void DrawBands(const std::string& options = "same");

    /**
     * Draw Markers: global mode, local mode, etc. */
    virtual void DrawMarkers();

    /**
     * Return a graph of the profile along x or y. The profile is
     * calculated by scanning through the one axis, e.g., y, and
     * finding the mode, mean, median, etc. with respect to the
     * other axis, e.g. x.
     * @param axis x-axis (0) or y-axis (1)
     * @param pt Type of profile to construct. */
    TGraph* CalculateProfileGraph(BCH2DProfileAxis axis, BCH2DProfileType pt = kProfileMean);

    /**
     * Draw the profiles along x and y */
    void DrawProfileGraphs();

    /** @} */

protected:

    /**
     * band type. */
    BCH2DBandType fBandType;

    /**
     * flag for drawing log z. */
    bool fLogz;

    /**
     * flag for drawing x profile. */
    bool fDrawProfileX;

    /**
     * profile type of x profile. */
    BCH2DProfileType fProfileXType;

    /**
     * color of x profile. */
    int fProfileXLineColor;

    /**
     * line style of x profile. */
    int fProfileXLineStyle;

    /**
     * flag for drawing y profile. */
    bool fDrawProfileY;

    /**
     * profile type of y profile. */
    BCH2DProfileType fProfileYType;

    /**
     * color of y profile. */
    int fProfileYLineColor;

    /**
     * line style of y profile. */
    int fProfileYLineStyle;

};

// ---------------------------------------------------------

#endif
