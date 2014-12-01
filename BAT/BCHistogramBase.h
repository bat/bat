#ifndef __BCHistogramBase__H
#define __BCHistogramBase__H

/*!
 * \class BCHistogramBase
 * \brief A base class for drawing histograms in BAT style
 * \author Daniel Greenwald
 * \version 1.0
 */

/*
 * Copyright (C) 2007-2014, the BAT core developer team
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 * For documentation see http://mpp.mpg.de/bat
 */

// ---------------------------------------------------------

#include <vector>

#include <TH1.h>

class TLegend;

// ---------------------------------------------------------

class BCHistogramBase {
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
  BCHistogramBase(TH1 * hist = 0, int dimension=0);

	/**
	 * Copy constructor. */
	BCHistogramBase(const BCHistogramBase & other);

  /**
   * The default destructor. */
  virtual ~BCHistogramBase();

  /** @} */

  /** \name Member functions (get)  */
  /** @{ */

  /**
   * @return The histogram. */
  TH1 * GetHistogram()
	{ return fHistogram; };

	/**
	 * @return The legend. */
	TLegend * GetLegend()
	{ return fLegend; }

	/**
	 * @return The global mode. */
	std::vector<double> GetGlobalMode()
	{ return fGlobalMode; }

	/**
	 * @return i'th component of global mode.
	 * @param i index of coordinate to return. */
	double GetGlobalMode(unsigned i) 
	{ return (i<fGlobalMode.size()) ? fGlobalMode[i] : 0; }

	/**
	 * @return The local mode. */
	std::vector<double> GetLocalMode()
	{ return fLocalMode; }

	/**
	 * @return i'th component of local mode.
	 * @param i index of coordinate to return. */
	double GetLocalMode(unsigned i) 
	{ return (i<fLocalMode.size()) ? fLocalMode[i] : 0; }

  /**
   * Returns a band color of the current color scheme.
   * @param index the color index
   * @return the color number. */
  int GetBandColor(int index) const
  { return fBandColors.at(index); };

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
	 * @return whether to draw local mode. */
	bool GetDrawLocalMode() const
	{ return fDrawLocalMode; }

	/**
	 * @return whether to draw global mode arrow. */
	bool GetDrawLocalModeArrows() const
	{ return fDrawLocalMode and fDrawLocalModeArrows; }

	/**
	 * @return whether to draw mean. */
	bool GetDrawMean() const
	{ return fDrawMean; }
	
	/**
	 * @return whether to draw standard deviation. */
	bool GetDrawStandardDeviation() const
	{ return fDrawMean and fDrawStandardDeviation; }

	/**
	 * @return flag for plotting log on x axis. */
	bool GetLogx()
	{ return fLogx; }

	/**
	 * @return flag for plotting log on y axis. */
	bool GetLogy()
	{ return fLogy; }

	/**
	 * @return vector of intervals to draw. */
	std::vector<double> GetIntervals()
	{ return fIntervals; }

  /** @} */

  /** \name Member functions (set)  */
  /** @{ */

	/**
	 * Copy options from. */
	virtual void CopyOptions(const BCHistogramBase & other);

	/**
	 * Set global mode. */
	void SetGlobalMode(std::vector<double> gm)
	{ fGlobalMode = gm; }

	/**
	 * Unset global mode. */
	void UnsetGlobalMode()
	{ fGlobalMode.clear(); }

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
	{ if (i<fBandColors.size()) fBandColors[i] = c; }

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
  virtual void SetHistogram(TH1 * hist);

	/**
	 * Sets drawing of x axis in log. */
	void SetLogx(bool flag=true)
	{ fLogx = flag; }

	/**
	 * Sets drawing of x axis in log. */
	void SetLogy(bool flag=true)
	{ fLogy = flag; }
	
	/**
	 * Sets number of credibility interval bands to draw. */
	void SetNBands(unsigned n)
	{ fNBands = n; }

	/**
	 * Sets number of times to smooth the histogram using ROOT's smoothing function. */
	void SetNSmooth(unsigned n)
	{ fNSmooth = n; }

	/**
	 * Set drawing of global mode. */
	void SetDrawGlobalMode(bool flag=true, bool arrows=true)
	{ fDrawGlobalMode = flag; fDrawGlobalModeArrows=arrows;}

	/**
	 * Set drawing of global mode. */
	void SetDrawLocalMode(bool flag=true, bool arrows=true)
	{ fDrawLocalMode = flag; fDrawLocalModeArrows=arrows;}

	/**
	 * Set drawing of mean.
	 * @param flag Toggles mean drawing.
	 * @param stddev Toggles standard deviaton drawing. (Automatically suppressed if mean is suppressed.) */
	void SetDrawMean(bool flag=true, bool stddev=true)
	{ fDrawMean = flag; fDrawStandardDeviation = stddev;}

	/**
	 * Set drawing of legend. */
	void SetDrawLegend(bool flag=true)
	{ fDrawLegend = flag; }

	/**
	 * Set drawing of ROOT histogram stats box. */
	void SetStats(bool flag=true)
	{ fDrawStats = flag; }

	/**
	 * Set intervals to be drawn. */
	void SetIntervals(std::vector<double> intervals) 
	{ fIntervals = intervals; }

	/** Add interval value. */
	void AddInterval(double interval) 
	{ fIntervals.push_back(interval); }

  /** @} */

  /** \name Member functions (miscellaneous methods) */
  /** @{ */

	void ClearBandColors()
	{ fBandColors.clear(); }

	/**
	 * Applying ROOT smoothing to histogram, and renormalize.
	 * @param n Number of times to smooth. */
	void Smooth(unsigned n);

  /**
   * Draw distribution into the active pad.
   * @param options ROOT drawing options
   * @param intervals the intervals
   */
  virtual void Draw(std::string options="", std::vector<double> intervals=std::vector<double>(0))
	{ }

  /**
   * Draw distribution into the active pad.
   * @param options ROOT drawing options
   * @param interval a single interval.
   */
  virtual void Draw(std::string options, double interval)
	{ Draw(options,std::vector<double>(1,interval)); }

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
	 * Resize histogram and draw legend. */
	virtual void DrawLegend();

  /** @} */

protected:

  /**
   * The histogram */
  TH1 * fHistogram;

	/**
	 * Legend. */
	TLegend * fLegend;

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
	 * The fill style for the bands. */
	short fBandFillStyle;

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
	 * Flag for drawing local mode. */
	bool fDrawLocalMode;

	/**
	 * Flag for drawing local mode arrows. */
	bool fDrawLocalModeArrows;

	/**
	 * Flag for drawing mean. */
	bool fDrawMean;

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
   * Storage for plot objects. */
  mutable std::vector<TObject*> fROOTObjects;

};

// ---------------------------------------------------------

#endif
