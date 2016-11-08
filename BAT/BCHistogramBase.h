#ifndef __BCHistogramBase__H
#define __BCHistogramBase__H

/*!
 * \class BCHistogramBase
 * \brief A base class for drawing histograms in BAT style
 * \author Daniel Greenwald
 * \version 1.0
 */

/*
 * Copyright (C) 2007-2015, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>
#include <string>
#include <utility>

#include <TLegend.h>

class TH1;
class TLegendEntry;
class TObject;

// ---------------------------------------------------------

class BCHistogramBase
{
public:

    /**	\name Enumerators */
    /** @{ */

    /**
     * An enumerator for the color schemes. */
    enum BCHColorScheme {
        kBlackWhite     = 0,
        kGreenYellowRed = 1,
        kBlueOrange     = 2,
        kRedGreen       = 3
    };

    /** } */

    /** \name Constructors and destructors */
    /** @{ */

    /**
     * The default constructor.
     * @param hist Histogram to build BCHistogram around.
     * @param dimension Dimension of the histogram (for checking purposes), if negative, no check is made. */
    BCHistogramBase(const TH1* const hist = 0, int dimension = 0);

    /**
     * Copy constructor. */
    BCHistogramBase(const BCHistogramBase& other);

    /**
     * The default destructor. */
    virtual ~BCHistogramBase();

    /** @} */

    /** \name Operators and swap*/
    /** @{ */

    /**
     * Assignment operator. */
    BCHistogramBase& operator=(BCHistogramBase other);

    /**
     * Swap function. */
    friend void swap(BCHistogramBase& first, BCHistogramBase& second);

    /**
     * Assign TH1 histogram with operator. */
    BCHistogramBase& operator=(const TH1* const hist)
    { SetHistogram(hist); return *this; }

    /**
     * Assign TH1 histogram with operator. */
    BCHistogramBase& operator=(const TH1& hist)
    { SetHistogram(&hist); return *this; }

    /** @} */

    /** \name Member functions (get)  */
    /** @{ */

    /**
     * @return The histogram. */
    TH1* GetHistogram()
    { return fHistogram; };

    /**
     * @return Dimensions of BC Histogram. */
    unsigned GetDimension() const
    { return fDimension; }

    /**
     * @return The legend. */
    TLegend& GetLegend()
    { return fLegend; }

    /**
     * @return the number of columns to be set into the legend. */
    unsigned GetNLegendColumns() const
    { return fNLegendColumns; }

    /**
     * @return The global mode. */
    std::vector<double>&  GetBestFitParameters()
    { return fGlobalMode; }

    /**
     * @return The global mode. */
    const std::vector<double>&  GetBestFitParameters() const
    { return fGlobalMode; }

    /**
     * @return i'th component of global mode.
     * @param i index of coordinate to return. */
    double GetBestFitParameters(unsigned i) const
    { return fGlobalMode.at(i); }

    /**
     * @return The local mode. */
    std::vector<double>& GetLocalMode()
    { return fLocalMode; }

    /**
     * @return The local mode. */
    const std::vector<double>& GetLocalMode() const
    { return fLocalMode; }

    /**
     * @return i'th component of local mode.
     * @param i index of coordinate to return. */
    double GetLocalMode(unsigned i)
    { return fLocalMode.at(i); }

    /**
     * Returns a band color of the current color scheme.
     * @param index the color index
     * @return the color number. */
    int GetBandColor(int index) const
    { return fBandColors.at(index); };

    /**
     * @return flag for whether bands should be drawn to overcover (true) or undercover (false). */
    bool GetBandOvercoverage() const
    { return fBandOvercoverage; }

    /**
     * @return histogram line color. */
    int GetLineColor() const
    { return fLineColor; }

    /**
     * Returns the marker colors (used for mean, median, and mode. */
    int GetMarkerColor() const
    { return fMarkerColor; }

    /**
     * @return scale of marker size. */
    double GetMarkerScale() const
    { return fMarkerScale; }

    /**
     * @return Band fill style. */
    short GetBandFillStyle() const
    { return fBandFillStyle; }

    /**
     * @return N bands to draw. */
    unsigned GetNBands() const
    { return fNBands; }

    /**
     * @return N times to smooth. */
    unsigned GetNSmooth() const
    { return fNSmooth; }

    /**
     * @return whether to draw global mode. */
    bool GetDrawGlobalMode() const
    { return fDrawGlobalMode; }

    /**
     * @return whether to draw global mode arrow. */
    bool GetDrawGlobalModeArrows() const
    { return fDrawGlobalMode and fDrawGlobalModeArrows; }

    /**
     * @return Global mode marker style. */
    int GetBestFitParametersMarkerStyle() const
    { return fGlobalModeMarkerStyle; }

    /**
     * @return whether to draw local mode. */
    bool GetDrawLocalMode() const
    { return fDrawLocalMode; }

    /**
     * @return whether to draw global mode arrow. */
    bool GetDrawLocalModeArrows() const
    { return fDrawLocalMode and fDrawLocalModeArrows; }

    /**
     * @return Local mode marker style. */
    int GetLocalModeMarkerStyle() const
    { return fLocalModeMarkerStyle; }

    /**
     * @return whether to draw mean. */
    bool GetDrawMean() const
    { return fDrawMean; }

    /**
     * @return Mean marker style. */
    int GetMeanMarkerStyle() const
    { return fMeanMarkerStyle; }

    /**
     * @return whether to draw standard deviation. */
    bool GetDrawStandardDeviation() const
    { return fDrawMean and fDrawStandardDeviation; }

    /**
     * @return flag for plotting log on x axis. */
    bool GetLogx() const
    { return fLogx; }

    /**
     * @return flag for plotting log on y axis. */
    bool GetLogy() const
    { return fLogy; }

    /**
     * @return flag for plotting log on z axis. */
    bool GetLogz() const
    { return fLogz; }

    /**
     * @return flag for drawing grid on x axis. */
    bool GetGridx() const
    { return fGridx; }

    /**
     * @return flag for drawing grid on y axis. */
    bool GetGridy() const
    { return fGridy; }

    /**
     * @return vector of intervals to draw. */
    std::vector<double>& GetIntervals()
    { return fIntervals; }

    /**
     * @return vector of intervals to draw. */
    const std::vector<double>& GetIntervals() const
    { return fIntervals; }

    /**
     * @return ROOT drawing options. */
    std::string& GetROOToptions()
    { return fROOToptions; }

    /**
     * @return ROOT drawing options. */
    const std::string& GetROOToptions() const
    { return fROOToptions; }

    /** @} */

    /** \name Member functions (set)  */
    /** @{ */

    /**
     * Copy options from. */
    virtual void CopyOptions(const BCHistogramBase& other);

    /**
     * Set global mode. */
    void SetGlobalMode(std::vector<double> gm)
    { fGlobalMode = gm; }

    /**
     * Unset global mode. */
    void UnsetGlobalMode()
    { fGlobalMode.clear(); }

    /**
     * Set global mode element. */
    void SetGlobalMode(unsigned i, double lm)
    { if (i < fGlobalMode.size()) fGlobalMode[i] = lm; }

    /**
     * Set local mode. */
    void SetLocalMode(std::vector<double> lm)
    { fLocalMode = lm; }

    /**
     * Unset local mode. */
    void UnsetLocalMode()
    { fLocalMode.clear(); }

    /**
     * Set local mode element. */
    void SetLocalMode(unsigned i, double lm)
    { if (i < fLocalMode.size()) fLocalMode[i] = lm; }

    /**
     * Sets the color scheme.
     * @param scheme the scheme index \n
     * 0 : black and white
     * 1 : yellow-green-red
     * 2 : blueish colors
     * 2 : redish colors
     * 2 : blueish colors
     */
    virtual void SetColorScheme(BCHColorScheme scheme);

    /**
     * Add band color.
     * @param c Color to add. */
    void AddBandColor(int c)
    { fBandColors.push_back(c); }

    /**
     * Set band color.
     * @param i Index of color to set.
     * @param c Color to set it to. */
    void SetBandColor(unsigned i, int c)
    { if (i < fBandColors.size()) fBandColors[i] = c; }

    /**
     * Set band coverage to be overcoverage (true) or undercoverage (false). */
    void SetBandOvercoverage(bool flag = true)
    { fBandOvercoverage = flag; }

    /**
     * Set histogram line color.
     * @param c Color for histogram line. */
    void SetLineColor(int c)
    { fLineColor = c; }

    /**
     * Set marker color (used for mean, median, and mode).
     * @param c Color for markers. */
    void SetMarkerColor(int c)
    { fMarkerColor = c; }

    /**
     * Set marker size scale. */
    void SetMarkerScale(double s)
    { fMarkerScale = s; }

    /**
     * Set band fill style.
     * @param f Fill style for bands. */
    void SetBandFillStyle(short f)
    { fBandFillStyle = f; }

    /**
     * Sets the histogram. */
    virtual void SetHistogram(const TH1* const hist);

    /**
     * Sets drawing of x axis in log. */
    void SetLogx(bool flag = true)
    { fLogx = flag; }

    /**
     * Sets drawing of y axis in log. */
    void SetLogy(bool flag = true)
    { fLogy = flag; }

    /**
     * Sets drawing of z axis in log. */
    void SetLogz(bool flag = true)
    { fLogz = flag; }

    /**
     * Sets drawing of grid on x axis. */
    void SetGridx(bool flag = true)
    { fGridx = flag; }

    /**
     * Sets drawing of grid on y axis. */
    void SetGridy(bool flag = true)
    { fGridy = flag; }

    /**
     * Sets number of credibility interval bands to draw. */
    void SetNBands(unsigned n)
    { fNBands = n; fIntervals = DefaultIntervals(); }

    /**
     * Sets number of times to smooth the histogram using ROOT's smoothing function. */
    void SetNSmooth(unsigned n)
    { fNSmooth = n; }

    /**
     * Set drawing of global mode. */
    void SetDrawGlobalMode(bool flag = true, bool arrows = true)
    { fDrawGlobalMode = flag; fDrawGlobalModeArrows = arrows;}

    /**
     * Set global mode marker style.
     * @param s Marker style. */
    void SetGlobalModeMarkerStyle(int s)
    { fGlobalModeMarkerStyle = s; }

    /**
     * Set drawing of global mode. */
    void SetDrawLocalMode(bool flag = true, bool arrows = true)
    { fDrawLocalMode = flag; fDrawLocalModeArrows = arrows;}

    /**
     * Set Local mode marker style.
     * @param s Marker style. */
    void SetLocalModeMarkerStyle(int s)
    { fLocalModeMarkerStyle = s; }

    /**
     * Set drawing of mean.
     * @param flag Toggles mean drawing.
     * @param stddev Toggles standard deviaton drawing. (Automatically suppressed if mean is suppressed.) */
    void SetDrawMean(bool flag = true, bool stddev = true)
    { fDrawMean = flag; fDrawStandardDeviation = stddev;}

    /**
     * Set mean marker style.
     * @param s Marker style. */
    void SetMeanMarkerStyle(int s)
    { fMeanMarkerStyle = s; }

    /**
     * Set drawing of legend. */
    void SetDrawLegend(bool flag = true)
    { fDrawLegend = flag; }

    /**
     * Set number of columns in legend. */
    void SetNLegendColumns(unsigned n)
    { fNLegendColumns = n; }

    /**
     * Set drawing of ROOT histogram stats box. */
    void SetStats(bool flag = true)
    { fDrawStats = flag; }

    /**
     * Set intervals to be drawn. */
    void SetIntervals(std::vector<double> intervals)
    { fIntervals = intervals; }

    /**
     * Set intervals to one single value. */
    void SetInterval(double interval)
    { fIntervals.clear(); fIntervals.push_back(interval); }

    /** Add interval value. */
    void AddInterval(double interval)
    { fIntervals.push_back(interval); }

    /** Set ROOT drawing options. */
    void SetROOToptions(const std::string& options)
    { fROOToptions = options; }

    /** @} */

    /** \name Member functions (miscellaneous methods) */
    /** @{ */

    /**
     * Whether histogram has been set and filled.
     * @return whether BC Histogram object is valid */
    virtual bool Valid() const;

    void ClearBandColors()
    { fBandColors.clear(); }

    void ClearIntervals()
    { fIntervals.clear(); }

    /**
     * Applying ROOT smoothing to histogram, and renormalize.
     * @param n Number of times to smooth; fNSmooth, if n is negative. */
    void Smooth(int n = -1);

    /**
     * Draw distribution into the active pad.
     */
    virtual void Draw();

    /**
     * Draw bands. */
    virtual void DrawBands(const std::string& /*options*/)
    { }

    /**
     * Draw markers (global mode, local mode, etc.). */
    virtual void DrawMarkers();

    /**
     * Draw global mode. */
    virtual void DrawGlobalMode();

    /**
     * Draw global mode. */
    virtual void DrawLocalMode();

    /**
     * Draw mean and standard deviation. */
    virtual void DrawMean();

    /**
     * Resize legend and set it for placement at the top of the pad.
     * @return new lower y coordinate of legend*/
    virtual double ResizeLegend();

    /**
     * Resize histogram and draw legend. */
    virtual void DrawLegend();

    /**
     * Fill vector with values and integrals of nonzero bins sorted by value.
     * @param bin_dens_mass vector of bin densities (values) and masses (integrals).
     * @param sort order to sort in: -1 (default) decreasing; +1 increasing; 0 unsorted. */
    void GetNonzeroBinDensityMassVector(std::vector<std::pair<double, double> >& bin_dens_mass, int sort = -1);

    /**
     * Get probability density levels bounding from below the
     * smallest-interval levels with probability mass near that provided
     * in the argument. CAUTION: This function will sort all bins of the
     * histogram; if the histogram has many bins, this can be an
     * expensive calculation.
     * @param masses Vector of probability masses to bound.
     * @param overcoverage Flag for providing values that overcover (if true), or undercover (if false).
     * @return Vector of pairs (probability density, probability mass lower-bounded by it). */
    std::vector<std::pair<double, double> > GetSmallestIntervalBounds(std::vector<double> masses, bool overcoverage = true);

    /**
     * Get smallest interval sizes in dimensions of histogram: length (1D),
     * area (2D), volume (3D).
     * @param masses vector of probability masses defining smallest intervals.
     * @param overcoverage Flag for providing values that overcover (if true), or undercover (if false).
     * @return vector of smallest interval sizes. */
    virtual std::vector<double> GetSmallestIntervalSize(std::vector<double> masses, bool overcoverage = true);

    /**
     * Get smallest interval size in dimensions of histogram: length (1D),
     * area (2D), volume (3D).
     * @param mass probability mass defining smallest interval.
     * @param overcoverage Flag for providing values that overcover (if true), or undercover (if false).
     * @return smallest interval size. */
    virtual double GetSmallestIntervalSize(double mass, bool overcoverage = true);

    /**
     * Check intervals: remove values below 0 or above 1.
     * @param sort if positive, sort by increasing order (default);
     * if negative, sort by decreasing order; if zero, don't sort.*/
    virtual void CheckIntervals(std::vector<double>& intervals, int sort);

    /**
     * Return default intervals.
     * @param nbands nbands to give defaults for; if negative, use fNBands.
     * @return vector of default values for band intervals. */
    virtual std::vector<double> DefaultIntervals(int nbands = -1);

    /**
     * Add legend entry, checking first for unused extra entries. */
    TLegendEntry* AddLegendEntry(TObject* obj, const std::string& label, const std::string& options);

    /**
     * Add band legend entry, creating unused extra entries if necessary. */
    TLegendEntry* AddBandLegendEntry(TObject* obj, const std::string& label, const std::string& options);

    /** @} */

protected:

    /**
     * The histogram */
    TH1* fHistogram;

    /**
     * Legend. */
    TLegend fLegend;

    /**
     * number of columns to be set for legend. */
    unsigned fNLegendColumns;

    /**
     * Global mode */
    std::vector<double> fGlobalMode;

    /**
     * Local mode */
    std::vector<double> fLocalMode;

    /**
     * The colors of the color scheme. */
    std::vector<int> fBandColors;

    /**
     * flag for whether bands should overcover. */
    bool fBandOvercoverage;

    /**
     * The fill style for the bands. */
    short fBandFillStyle;

    /**
     * histogram line color. */
    int fLineColor;

    /**
     * Marker color. */
    int fMarkerColor;

    /**
     * Marker size scale. */
    double fMarkerScale;

    /**
     * Flag for controlling log plotting of x axis. */
    bool fLogx;

    /**
     * Flag for controlling log plotting of y axis. */
    bool fLogy;

    /**
     * Flag for controlling log plotting of z axis. */
    bool fLogz;

    /**
     * Flag for drawing of grid on x axis. */
    bool fGridx;

    /**
     * Flag for drawing of grid on y axis. */
    bool fGridy;

    /**
     * Number of credibility interval bands to draw. */
    unsigned fNBands;

    /**
     * Number of times to smooth histogram using ROOT's smooth function. */
    unsigned fNSmooth;

    /**
     * Flag for drawing global mode. */
    bool fDrawGlobalMode;

    /**
     * Flag for drawing global mode arrows. */
    bool fDrawGlobalModeArrows;

    /**
     * Holds option for global-mode marker style. */
    int fGlobalModeMarkerStyle;

    /**
     * Flag for drawing local mode. */
    bool fDrawLocalMode;

    /**
     * Flag for drawing local mode arrows. */
    bool fDrawLocalModeArrows;

    /**
     * Holds option for local-mode marker style. */
    int fLocalModeMarkerStyle;

    /**
     * Flag for drawing mean. */
    bool fDrawMean;

    /**
     * Holds option for mean marker style. */
    int fMeanMarkerStyle;

    /**
     * Flag for drawing standard deviation. */
    bool fDrawStandardDeviation;

    /**
     * Flag for drawing legend. */
    bool fDrawLegend;

    /**
     * Flag for drawing stats box. */
    bool fDrawStats;

    /**
     * Dimension of histogram. */
    int fDimension;

    /**
     * intervals to draw. */
    std::vector<double> fIntervals;

    /**
     * ROOT drawing options. */
    std::string fROOToptions;

    /**
     * Storage for plot objects. */
    mutable std::vector<TObject*> fROOTObjects;

    /**
     * Storage for unused legend entries (for use with multicolumn legends). */
    std::vector<TLegendEntry*> fExtraLegendEntries;

};

// ---------------------------------------------------------

#endif
